import sys
from pool import Pool

def go(filename, max_size, use_donor_weight):
    def edge_score_fun(e):
        return 1# + (e.patient.weight - e.donor.weight)/999.0

    sys.stdout.write(filename + "\t")
    sys.stdout.write("1\t" if use_donor_weight else "0\t")

    pool = Pool.from_file(filename)
    pool.create_donor_patient_arcs(use_donor_weight)
    sys.stdout.write(str(len(pool.patients)) + "\t")

    unrestricted_obj_val, optimal_exchanges = pool.solve_uef(edge_score_fun)
    sys.stdout.write(str(unrestricted_obj_val) + "\t")

    obj_val = None
    for i in range(2, max_size+1):
        if obj_val==unrestricted_obj_val:
            pass
        elif i < 6:
            obj_val, optimal_exchanges = pool.solve_cycle_formulation(i, edge_score_fun)
        else:
            obj_val, optimal_exchanges = pool.solve_eef(i, edge_score_fun)

        sys.stdout.write(str(obj_val) + "\t")
        sys.stdout.flush()

    sys.stdout.write(str(pool.n_patients_with_2_compat_donors()) + "\n")
    

if __name__=="__main__":
    filename = sys.argv[1]
    max_size = int(sys.argv[2])
    use_donor_weight = sys.argv[3].lower() == "true"

    go(filename, max_size, use_donor_weight)
