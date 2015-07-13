from nose.tools import *
from os import listdir
import os.path
from lungexchange.pool import *

def edge_score_fun(edge):
    return (1 + edge.donor.weight/1000.0 + edge.patient.weight/1000.0 +
            0.1*edge.donor.weight/edge.patient.weight)

def test_restricted_optimisation():
    test_dir = "test_instances"
    files = listdir(test_dir)
    for file in files:
        pool = Pool.from_file(os.path.join(test_dir, file))
        # Only test small instances to save time
        if len(pool.patients) < 20:
            pool.create_donor_patient_arcs()
            for max_size in range(2,6):
                opt, exchs = pool.solve_cycle_formulation(max_size, edge_score_fun)
                opt2, exchs2 = pool.solve_eef(max_size, edge_score_fun)
                assert_almost_equal(opt, opt2)
                assert_equal(exchs, exchs2)

def test_unrestricted_optimisation():
    test_dir = "test_instances"
    files = listdir(test_dir)
    for file in files:
        pool = Pool.from_file(os.path.join(test_dir, file))
        if len(pool.patients) < 20:
            pool.create_donor_patient_arcs()
            opt, exchs = pool.solve_uef(edge_score_fun)
            opt2, exchs2 = pool.solve_eef(len(pool.patients), edge_score_fun)
            assert_almost_equal(opt, opt2)
            assert_equal(exchs, exchs2)

def test_exch_finding_1_patient():
    # For simplicity, some patient and donor attributes aren't
    # specified here
    pool = Pool()
    p = Patient(None, None, None, None, None, 0)
    d0 = Donor(None, None, None, None, 0, p)
    d1 = Donor(None, None, None, None, 1, p)
    pool.patients.append(p)
    pool.donors.extend([d0, d1])
    p.paired_donors.extend([d0, d1])
    d0.compat_patients.append(p)
    d1.compat_patients.append(p)

    exchanges = pool.find_exchanges(10)
    assert_equal(len(exchanges), 1)
    assert_equal(exchanges, [Exchange([DPMatch(d0, p), DPMatch(d1, p)])])

    exchanges = pool.find_exchanges(1)
    assert_equal(len(exchanges), 1)
    assert_equal(exchanges, [Exchange([DPMatch(d0, p), DPMatch(d1, p)])])


def test_exch_finding_2_patients_a():
    # For simplicity, some patient and donor attributes aren't
    # specified here
    pool = Pool()
    p = Patient(None, None, None, None, None, 0)
    d0 = Donor(None, None, None, None, 0, p)
    d1 = Donor(None, None, None, None, 1, p)
    pool.patients.append(p)
    pool.donors.extend([d0, d1])
    p.paired_donors.extend([d0, d1])
    d0.compat_patients.append(p)
    d1.compat_patients.append(p)

    p1 = Patient(None, None, None, None, None, 1)
    d2 = Donor(None, None, None, None, 2, p1)
    d3 = Donor(None, None, None, None, 3, p1)
    pool.patients.append(p1)
    pool.donors.extend([d2, d3])
    p1.paired_donors.extend([d2, d3])
    d2.compat_patients.append(p)
    d3.compat_patients.append(p)
    d0.compat_patients.append(p1)
    d1.compat_patients.append(p1)

    exchanges = pool.find_exchanges(10)
    assert_equal(len(exchanges), 2)
    assert_equal(
        sorted(exchanges),
        sorted([
            Exchange([DPMatch(d0, p), DPMatch(d1, p)]),
            Exchange([DPMatch(d0, p1), DPMatch(d1, p1), DPMatch(d2, p), DPMatch(d3, p)])
        ]))

    exchanges = pool.find_exchanges(1)
    assert_equal(len(exchanges), 1)
    assert_equal(exchanges, [Exchange([DPMatch(d0, p), DPMatch(d1, p)])])

def test_exch_finding_2_patients_b():
    pool = Pool()
    p = Patient(None, None, None, None, None, 0)
    d0 = Donor(None, None, None, None, 0, p)
    d1 = Donor(None, None, None, None, 1, p)
    pool.patients.append(p)
    pool.donors.extend([d0, d1])
    p.paired_donors.extend([d0, d1])

    p1 = Patient(None, None, None, None, None, 1)
    d2 = Donor(None, None, None, None, 2, p1)
    d3 = Donor(None, None, None, None, 3, p1)
    pool.patients.append(p1)
    pool.donors.extend([d2, d3])
    p1.paired_donors.extend([d2, d3])

    d0.compat_patients.append(p)
    d2.compat_patients.append(p)
    d1.compat_patients.append(p1)
    d3.compat_patients.append(p1)

    exchanges = pool.find_exchanges(10)
    assert_equal(len(exchanges), 1)
    assert_equal(
        exchanges,
        [Exchange([DPMatch(d0, p), DPMatch(d1, p1), DPMatch(d2, p), DPMatch(d3, p1)])]
        )

    exchanges = pool.find_exchanges(1)
    assert_equal(len(exchanges), 0)

def test_exch_finding_3_patients_a():
    pool = Pool()
    p = Patient(None, None, None, None, None, 0)
    d0 = Donor(None, None, None, None, 0, p)
    d1 = Donor(None, None, None, None, 1, p)
    pool.patients.append(p)
    pool.donors.extend([d0, d1])
    p.paired_donors.extend([d0, d1])

    p1 = Patient(None, None, None, None, None, 1)
    d2 = Donor(None, None, None, None, 2, p1)
    d3 = Donor(None, None, None, None, 3, p1)
    pool.patients.append(p1)
    pool.donors.extend([d2, d3])
    p1.paired_donors.extend([d2, d3])

    p2 = Patient(None, None, None, None, None, 2)
    d4 = Donor(None, None, None, None, 4, p2)
    d5 = Donor(None, None, None, None, 5, p2)
    pool.patients.append(p2)
    pool.donors.extend([d4, d5])
    p2.paired_donors.extend([d4, d5])

    d0.compat_patients.append(p)
    d1.compat_patients.append(p1)
    d2.compat_patients.append(p)
    d3.compat_patients.append(p2)
    d4.compat_patients.append(p2)
    d5.compat_patients.append(p1)

    exchanges = pool.find_exchanges(10)
    assert_equal(len(exchanges), 1)
    assert_equal(
        exchanges,
        [Exchange([DPMatch(d0, p), DPMatch(d1, p1),
                   DPMatch(d2, p), DPMatch(d3, p2),
                   DPMatch(d4, p2), DPMatch(d5, p1)])]
        )

    exchanges = pool.find_exchanges(1)
    assert_equal(len(exchanges), 0)

