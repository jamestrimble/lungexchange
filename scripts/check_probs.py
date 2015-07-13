import sys
import pandas as pd

def start(patients_file, donors_file):
    patients = pd.read_table(patients_file, sep=' ',
            header=None, names=['condition',
            'sex', 'age', 'blood_type', 'weight'])

    donors = pd.read_table(donors_file, sep=' ',
            header=None, names=['patient_id',
            'sex', 'age', 'blood_type', 'weight'])

    print " *** Ratio of cystic fibrosis to pulm. hypertension"
    print " *** should be approximately 2.66 to 1."
    print (1. * patients.condition.value_counts() /
            patients.condition.value_counts()['PH'])

    weight_breaks = [
            99, 109, 119, 129, 139, 149, 159, 169, 179, 189, 199,
                209, 219, 229, 239, 249, 259, 269, 279, 289, 299,
            319, 339, 359, 379, 399, 419, 439, 500]

    for condition in ['CF', 'PH']:
        print " *** Check gender distribution for patients with {0}".format(condition)
        p_with_cond = patients[patients['condition']==condition].copy()
        print 1. * p_with_cond.sex.value_counts() / len(p_with_cond.sex)

        print " *** Check age distribution for patients with {0}".format(condition)
        age_breaks = [0, 17, 34, 49, 64, 100]
        p_with_cond['age_band'] = pd.cut(p_with_cond['age'], age_breaks)
        pt = pd.pivot_table(p_with_cond, index=['age_band'],
                columns=['sex'], values=['age'], aggfunc=[len])
        print 1. * pt / pt.sum()

        print " *** Check blood type distribution for patients with {0}".format(condition)
        print 1. * p_with_cond.blood_type.value_counts() / len(p_with_cond.blood_type)
        print

    def print_wt_dist(df):
        age_breaks = [0, 19, 29, 39, 49, 59, 69, 79, 100]
        df['age_band'] = pd.cut(df['age'], age_breaks)
        df['weight_band'] = pd.cut(df['weight'], weight_breaks)
        pt = pd.pivot_table(df, index=['weight_band'],
                columns=['sex', 'age_band'], values=['age'], aggfunc=[len])
        print 1. * pt.cumsum(0) / pt.sum()

    print " *** Check weight distribution for patients"
    print_wt_dist(patients)

    print " *** Check weight distribution for donors"
    print_wt_dist(donors)

    print " *** approximately 50 per cent of donors should be female"
    print (1. * donors.sex.value_counts() /
            len(donors.sex))

    print " *** Check donor blood type distribution"
    print (1. * donors.blood_type.value_counts() /
            len(donors.blood_type))

    print " *** Check donor age distribution"
    age_breaks = [0, 19, 24, 29, 34, 39, 44, 49, 54, 59, 64, 69, 74, 79, 100]
    donors['age_band'] = pd.cut(donors['age'], age_breaks)
    pt = pd.pivot_table(donors, index=['age_band'],
            columns=['sex'], values=['age'], aggfunc=[len])
    print 1. * pt / pt.sum()

if __name__=="__main__":
    start(sys.argv[1], sys.argv[2])
