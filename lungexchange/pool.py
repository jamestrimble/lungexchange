from itertools import combinations
from collections import namedtuple, defaultdict
import sys
from ortools.linear_solver import pywraplp

class DPMatch(object):
    def __init__(self, donor, patient):
        self.donor = donor
        self.patient = patient

    def __gt__(self, other):
        if self.donor.d_id > other.donor.d_id:
            return True
        elif self.donor.d_id==other.donor.d_id and self.patient.p_id>other.patient.p_id:
            return True
        else:
            return False

    def __eq__(self, other):
        return self.donor==other.donor and self.patient==other.patient

class Patient(object):
    def __init__(self, disease, gender, age, bt, weight, p_id):
        self.disease = disease
        self.gender = gender
        self.age = age
        self.bt = bt
        self.weight = weight
        self.p_id = p_id
        self.paired_donors = []

    @classmethod
    def from_string(cls, s, p_id):
        disease, gender, age, bt, weight = s.split()
        age = int(age)
        weight = int(weight)
        return cls(disease, gender, age, bt, weight, p_id)

    def __str__(self):
        return "{0} {1} {2} {3} {4}".format(
            self.disease, self.gender, self.age, self.bt, self.weight)

    def n_compatible_paired_donors(self):
        return sum(self in d.compat_patients for d in self.paired_donors)

    def draw_number_of_donors(self):
        # Generates the number of donors that should
        # be created for this patient
        return 2

class Donor(object):
    def __init__(self, gender, age, bt, weight, d_id, patient):
        self.gender = gender
        self.age = age
        self.bt = bt
        self.weight = weight
        self.d_id = d_id
        self.patient = patient # TODO: avoid keeping this reference
        self.compat_patients = []

    @classmethod
    def from_string(cls, s, d_id, patients):
        p_id, gender, age, bt, weight = s.split()
        p_id = int(p_id)
        age = int(age)
        weight = int(weight)
        patient = patients[p_id]
        return cls(gender, age, bt, weight, d_id, patient)

    def is_compatible(self, patient, use_weight=True):
        bt_compat = (self.bt=="O" or
                     patient.bt=="AB" or
                     self.bt==patient.bt)
        size_compat = self.weight >= patient.weight or not use_weight
        return bt_compat and size_compat

    def __str__(self):
        return "{0} {1} {2} {3} {4}".format(
            self.patient.p_id, self.gender, self.age, self.bt, self.weight)


class Exchange(object):
    # TODO: implement __hash__(self)

    def __init__(self, edges):
        self.edges = sorted(edges)
    
    def show(self):
        for edge in self.edges:
            print edge.donor.d_id, edge.patient.p_id
        print

    def __eq__(self, other):
        return self.edges==other.edges

    def __gt__(self, other):
        if len(self.edges) > len(other.edges):
            return True
        elif len(self.edges) < len(other.edges):
            return False
        else:
            for e0, e1 in zip(self.edges, other.edges):
                if e0 > e1:
                    return True
            return False

class OptimisationException(Exception):
    pass

class InvalidInputException(Exception):
    pass

def create_solver():
    return pywraplp.Solver('CoinsGridCLP',
            pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)

class Pool(object):
    def __init__(self):
        self.patients = []
        self.donors = []

    @classmethod
    def from_file(cls, filename):
        with open(filename) as f:
            n_patients, n_donors = f.readline().strip().split()
            n_patients = int(n_patients)
            n_donors = int(n_donors)

            pool = Pool()

            for i in range(n_patients):
                s = f.readline().strip()
                pool.patients.append(Patient.from_string(s, i))

            for i in range(n_donors):
                s = f.readline().strip()
                donor = Donor.from_string(s, i, pool.patients)
                pool.donors.append(donor)
                pool.patients[donor.patient.p_id].paired_donors.append(donor)

            return pool

    def n_patients_with_2_compat_donors(self):
        return sum(p.n_compatible_paired_donors()>=2 for p in self.patients)

    def create_donor_patient_arcs(self, use_weight=True):
        """
        For each donor, add to the donor's compat_patients list those
        patients to whom the donor can donate. An exception is made for patients
        who have two or more compatible paired donors. Such a patient's paired
        donors have no outgoing arcs to other patients. Furthermore, such a
        patient has no incoming arcs from other patients' paired donors.
        """
        
        # Create a set of patients who can receive transplants from
        # both of their paired donors
        direct_donation_patients = set()

        for patient in self.patients:
            if sum(donor.is_compatible(patient, use_weight)
                         for donor in patient.paired_donors) >= 2:
                direct_donation_patients.add(patient)
                for donor in patient.paired_donors:
                    if donor.is_compatible(patient):
                        donor.compat_patients.append(patient)

        for patient in self.patients:
            if patient not in direct_donation_patients:
                for donor in patient.paired_donors:
                    for patient2 in self.patients:
                        if (patient2 not in direct_donation_patients and
                                    donor.is_compatible(patient2, use_weight)):
                            donor.compat_patients.append(patient2)


    def find_exchanges(self, max_size):
        exchanges = []   # Return value: a list of feasible exchanges
        
        # adj_mat[i][j] is the edge if there is an edge from donor i
        # to patient j, and None otherwise
        adj_mat = [[None for patient in self.patients] for donor in self.donors]
        for donor in self.donors:
            for patient in donor.compat_patients:
                adj_mat[donor.d_id][patient.p_id] = DPMatch(donor, patient)

        def can_be_added(patient, in_edge_counts, q):
            """
            Can patient be added, given the list of edges already traversed
            and the queue of patients to traverse from (past, present and future) 
            """
            if patient.p_id < start_p_id:
                # Patient's ID must be at least as great as the the ID of the
                # first patient visited
                return False
            if patient in q:
                return in_edge_counts[patient.p_id] <= 1
            else:
                return len(q) < max_size

        def find_unfull_patients(q, in_edge_counts):
            """
            Returns a list of two (possibly equal) patients who still need an in-edge.
            This is only intended to work for the last two edges to make up an exchange
            of the maximum permiteed size.
            """
            retval = []
            for patient in q:
                in_edge_count = in_edge_counts[patient.p_id]
                if in_edge_count == 1:
                    retval.append(patient)
                elif in_edge_count == 0:
                    return (patient, patient)
            return tuple(retval)

        def find(edges, q, q_current, in_edge_counts):
            # edges: list of edges already traversed
            # q: a list of patients that have been explored or are to be explored
            # q_current: the index of the current patient in q,
            #            from which edges are to be explored
            patient = q[q_current]

            if q_current == max_size-1:
                # An optimised way to check if the last two edges can be added
                p0, p1 = find_unfull_patients(q, in_edge_counts)
                for d0, d1 in combinations(patient.paired_donors, 2):
                    if adj_mat[d0.d_id][p0.p_id] and adj_mat[d1.d_id][p1.p_id]:
                        exchanges.append(Exchange(
                                edges + [adj_mat[d0.d_id][p0.p_id], adj_mat[d1.d_id][p1.p_id]]))
                    if p0!=p1 and adj_mat[d0.d_id][p1.p_id] and adj_mat[d1.d_id][p0.p_id]:
                        exchanges.append(Exchange(
                                edges + [adj_mat[d0.d_id][p1.p_id], adj_mat[d1.d_id][p0.p_id]]))

            else:
                for d0, d1 in combinations(patient.paired_donors, 2):
                    for p0 in d0.compat_patients:
                        if can_be_added(p0, in_edge_counts, q):
                            edges.append(DPMatch(d0, p0))
                            in_edge_counts[p0.p_id] += 1
                            p0_appended = p0 not in q
                            if p0_appended: q.append(p0)
                            for p1 in d1.compat_patients:
                                if can_be_added(p1, in_edge_counts, q):
                                    edges.append(DPMatch(d1, p1))
                                    in_edge_counts[p1.p_id] += 1
                                    p1_appended = p1 not in q
                                    if p1_appended: q.append(p1)
                                    if len(edges)==len(q)*2:
                                        exchanges.append(Exchange(edges[:]))
                                    elif len(q) > q_current+1:
                                        find(edges, q, q_current+1, in_edge_counts)
                                    del edges[-1]
                                    in_edge_counts[p1.p_id] -= 1
                                    if p1_appended: del q[-1]
                            del edges[-1]
                            in_edge_counts[p0.p_id] -= 1
                            if p0_appended: del q[-1]

        for start_p in self.patients:
            start_p_id = start_p.p_id  # This is frequently used, so keep it in a variable
            find([], [start_p], 0, [0]*len(self.patients))

        return exchanges

    def solve_cycle_formulation(self, max_size, edge_score_fun):
        s = create_solver()

        exchanges = self.find_exchanges(max_size)

        # For each set of patients, keep track of the best available exchange
        best_exchs = {}
        for exch in exchanges:
            patients = frozenset(e.patient for e in exch.edges)
            exch.score = sum(edge_score_fun(e) for e in exch.edges)
            if patients not in best_exchs or best_exchs[patients].score < exch.score:
                best_exchs[patients] = exch

        exchanges = best_exchs.values()

        for p in self.patients:
            p.mip_vars = []

        for i, exch in enumerate(exchanges):
            exch.mip_var = s.IntVar(0, 1, 'e' + str(i))
            for p in set(edge.patient for edge in exch.edges):
                p.mip_vars.append(exch.mip_var)
        
        for p in self.patients:
            s.Add(s.Sum(p.mip_vars) <= 1)

        z = s.Sum(sum(edge_score_fun(e) for e in exch.edges)*exch.mip_var
                        for exch in exchanges)
        s.Maximize(z)
        solve_status = s.Solve()
        if solve_status != s.OPTIMAL:
            raise OptimisationException("Solver status was " + str(solve_status))

        optimal_exchanges = [exch for exch in exchanges
                if exch.mip_var.SolutionValue() > 0.5]
        optimal_exchanges.sort()
        return s.Objective().Value(), optimal_exchanges

    def _create_edges(self):
        edges = []
        for d in self.donors:
            for p in d.compat_patients:
                edges.append(DPMatch(d, p))

        return edges

    def _create_exchanges_from_edge_list(self, edges):
        # Each element of edge_lists will be a pair,
        # with the first element being the set of patients and the
        # second element being the list of edges
        edge_lists = []
        def try_to_add_to_existing_list(edge):
            for l in edge_lists:
                if edge.patient in l[0] or edge.donor.patient in l[0]:
                    l[1].append(edge)
                    l[0].add(edge.patient)
                    l[0].add(edge.donor.patient)
                    return True
            return False
        for edge in edges:
            if not try_to_add_to_existing_list(edge):
                edge_lists.append(({edge.patient, edge.donor.patient}, [edge]))

        # Combine edge lists with overlapping patents
        combined_edge_lists = []
        for l in edge_lists:
            found_list = False
            for l1 in combined_edge_lists:
                if l[0] & l1[0]: # If patient-sets overlap
                    found_list = True
                    l1[0].update(l[0])
                    l1[1].extend(l[1])
            if not found_list:
                combined_edge_lists.append(l)

        return [Exchange(edge_list[1]) for edge_list in combined_edge_lists]

    def solve_eef(self, max_size, edge_score_fun):
        for patient in self.patients:
            if len(patient.paired_donors) > 2:
                raise InvalidInputException("Patient " + str(patient.p_id) +
                        " has more than two donors.")

        n_patients = len(self.patients)
        n_donors = len(self.donors)

        s = create_solver()

        # in_vars[i][j] is a list of variables representing edges to
        # patient i graph copy j
        in_vars  = [[[] for low_id in range(n_patients)] for p_id in range(n_patients)]
        out_vars = [[[] for low_id in range(n_patients)] for p_id in range(n_patients)]

        # out_from_d_in_copy[i][j] is a list of variables representing
        # edges from donor i in graph copy j
        out_from_d_in_copy = [[[] for low_id in range(n_patients)] for d_id in range(n_donors)]
        
        in_vars_to_p = [[] for p_id in range(n_patients)]
        out_vars_from_d = [[] for d_id in range(n_donors)]
        
        edges = self._create_edges()
        for edge in edges:
            # Get the lower patient ID in this edge
            low_id = min(edge.patient.p_id, edge.donor.patient.p_id)
            # Add a MIP var corresponding to this edge in each copy
            # of the graph up to low_id
            edge.mip_vars = {}
            for i in range(0, low_id+1):
                mip_var = s.IntVar(0, 1,
                    'edge' + str(i) + "-" + str(edge.donor.d_id) +
                    "-" + str(edge.patient.p_id))
                edge.mip_vars[i] = mip_var
                in_vars[edge.patient.p_id][i].append(mip_var)
                out_vars[edge.donor.patient.p_id][i].append(mip_var)
                in_vars_to_p[edge.patient.p_id].append(mip_var)
                out_vars_from_d[edge.donor.d_id].append(mip_var)
                out_from_d_in_copy[edge.donor.d_id][i].append(mip_var)

        # Constraint 1:
        # Within each graph copy, for each patient we have
        # count of in-edges == count of out-edges
        for low_p_id in range(n_patients):
            for p_id in range(low_p_id, n_patients):
                s.Add(s.Sum(in_vars[p_id][low_p_id]) == s.Sum(out_vars[p_id][low_p_id]))

        # Constraint 2:
        # At most 2*max_size edges are chosen in each graph copy
        for low_p_id in range(n_patients):
            mip_vars = []
            for edge in edges:
                if low_p_id in edge.mip_vars:
                    mip_vars.append(edge.mip_vars[low_p_id])
            s.Add(s.Sum(mip_vars) <= 2*max_size)

        # Constraint 3:
        # Each donor has <=1 outgoing arc in total
        for d_id in range(0, n_donors):
            s.Add(s.Sum(out_vars_from_d[d_id]) <= 1)

        # Constraint 4:
        # In each graph copy, for each patient p, the number of outgoing arcs from p's
        # first donor equals the number of outgoing arcs from p's
        # second donor
        for low_p_id in range(n_patients):
            for patient in self.patients:
                d_id0 = patient.paired_donors[0].d_id
                d_id1 = patient.paired_donors[1].d_id
                s.Add(s.Sum(out_from_d_in_copy[d_id0][low_p_id]) ==
                      s.Sum(out_from_d_in_copy[d_id1][low_p_id]))
        
        # Objective
        objective_terms = []
        for edge in edges:
            for mip_var in edge.mip_vars.values():
                objective_terms.append(edge_score_fun(edge) * mip_var)
        z = s.Sum(objective_terms)
        s.Maximize(z)
        solve_status = s.Solve()
        if solve_status != s.OPTIMAL:
            raise OptimisationException("Solver status was " + str(solve_status))

        optimal_edges = []
        for edge in edges:
            for low_id in edge.mip_vars:
                if edge.mip_vars[low_id].SolutionValue() > 0.5:
                    optimal_edges.append(edge)
        
        optimal_exchanges = self._create_exchanges_from_edge_list(optimal_edges)
        optimal_exchanges.sort()
        return s.Objective().Value(), optimal_exchanges


    def solve_uef(self, edge_score_fun):
        n_patients = len(self.patients)
        n_donors = len(self.donors)

        s = create_solver()

        out_vars_from_p = [[] for p_id in range(n_patients)]
        in_vars_to_p = [[] for p_id in range(n_patients)]
        out_vars_from_d = [[] for d_id in range(n_donors)]
        
        edges = self._create_edges()
        for edge in edges:
            mip_var = s.IntVar(0, 1,
                'edge' + str(edge.donor.d_id) +
                "-" + str(edge.patient.p_id))
            edge.mip_var = mip_var
            out_vars_from_p[edge.donor.patient.p_id].append(mip_var)
            in_vars_to_p[edge.patient.p_id].append(mip_var)
            out_vars_from_d[edge.donor.d_id].append(mip_var)

        # Constraint 1:
        # For each patient we have
        # count of in-edges == count of out-edges
        for p_id in range(0, n_patients):
            s.Add(s.Sum(in_vars_to_p[p_id]) == s.Sum(out_vars_from_p[p_id]))

        # Constraint 2:
        # Each patient has 2 or 0 incoming arcs
        for p_id in range(0, n_patients):
            b = s.IntVar(0, 1, 'b' + str(p_id))
            s.Add(2*b + s.Sum(in_vars_to_p[p_id]) == 2)

        # Constraint 3:
        # Each donor has <=1 outgoing arc
        for d_id in range(0, n_donors):
            s.Add(s.Sum(out_vars_from_d[d_id]) <= 1)

        # Objective
        objective_terms = []
        for edge in edges:
            objective_terms.append(edge_score_fun(edge) * edge.mip_var)
        z = s.Sum(objective_terms)
        s.Maximize(z)
        solve_status = s.Solve()
        if solve_status != s.OPTIMAL:
            raise OptimisationException("Solver status was " + str(solve_status))

        optimal_edges = []
        for edge in edges:
            if edge.mip_var.SolutionValue() > 0.5:
                optimal_edges.append(edge)
        
        optimal_exchanges = self._create_exchanges_from_edge_list(optimal_edges)
        optimal_exchanges.sort()
        return s.Objective().Value(), optimal_exchanges

    def show(self):
        for patient in self.patients:
            print patient.p_id, "---", str(patient)
        for donor in self.donors:
            print str(donor), ":", [p.p_id for p in donor.compat_patients]


