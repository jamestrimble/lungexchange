import sys
from pool import Pool

def start(filename, max_size):
    pool = Pool.from_file(filename)

    pool.create_donor_patient_arcs(True)

    pool.show()

    exchanges = pool.find_exchanges(max_size)
    for exch in exchanges:
        exch.show()

    print "Number of exchanges: ", len(exchanges)

    # Find out how many exchanges there are if two exchanges with the same
    # set of patients are considered equal
    unique_exch = set()
    for exch in exchanges:
        unique_exch.add(frozenset(e.patient for e in exch.edges))
    print "....", len(unique_exch)

    def edge_score_fun(e):
        return 1 + (e.patient.weight - e.donor.weight)/999.0

    print "Optimal exchanges with restricted cycle size (cycle formulation):"
    obj_val, optimal_exchanges = pool.solve_cycle_formulation(max_size, edge_score_fun)
    print obj_val
    for exch in optimal_exchanges:
        exch.show()

    print "Optimal exchanges with restricted cycle size:"
    obj_val, optimal_exchanges = pool.solve_eef(max_size, edge_score_fun)
    print obj_val
    for exch in optimal_exchanges:
        exch.show()

    print "Optimal exchanges with unrestricted cycle size:"
    obj_val, optimal_exchanges = pool.solve_uef(edge_score_fun)
    print obj_val
    for exch in optimal_exchanges:
        exch.show()

    print "Number of patients with 2 compatible paired donors:"
    print pool.n_patients_with_2_compat_donors()
    # TODO: put this in a unit test
    print pool.solve_eef(1, lambda e: 1)[0]/2

if __name__=="__main__":
    start(sys.argv[1], int(sys.argv[2]))

