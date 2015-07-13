import random
import json
import sys
from discrete_prob import DiscreteDist
from pool import Patient, Donor

class WeightDistributions(object):
    def __init__(self, weight_dist_data):
        self.weight_dist_data = weight_dist_data
        for w_dist in weight_dist_data:
            w_dist["wt_dist"] = DiscreteDist(w_dist["weight_bands"],
                                   w_dist["cum_probs"],
                                   are_cumulative=True)

    def draw_weight(self, gender, age):
        for w_dist in self.weight_dist_data:
            if w_dist["sex"]==gender and age <= w_dist["max_age"]:
                wt_band = w_dist["wt_dist"].sample1()
                return randint(wt_band["min"], wt_band["max"])

def generate_patients(d, n_patients, weight_dists):
    patients = []

    for p_id in range(n_patients):
        if rand() < d["prob_cystic_fibrosis"]:
            disease = "CF"
            dd = d["cystic_f_dists"]
        else:
            disease = "PH"
            dd = d["pulm_h_dists"]

        if rand() < dd["prob_female"]:
            gender = "F"
            age_d = dd["age_dist_f"]
        else:
            gender = "M"
            age_d = dd["age_dist_m"]

        age_dist = DiscreteDist(age_d["age_bands"], age_d["probs"])
        age_band = age_dist.sample1()
        age = randint(age_band["min"], age_band["max"])

        bt_d = dd["bt_dist"]
        bt_dist = DiscreteDist(bt_d["bts"], bt_d["probs"])
        bt = bt_dist.sample1()

        weight = weight_dists.draw_weight(gender, age)

        patient = Patient(disease, gender, age, bt, weight, p_id)
        patients.append(patient)
    return patients

def generate_donors(d, patients, weight_dists):
    donors = []
    d_id = 0
    for patient in patients:
        n_donors = patient.draw_number_of_donors()
        for i in range(n_donors):
            if random.random() < 0.5:
                sex = "F"
                age_d = d["donor_age_dist_f"]
            else:
                sex = "M"
                age_d = d["donor_age_dist_m"]
            age_dist = DiscreteDist(age_d["ages"], age_d["probs"])
            age = age_dist.sample1()

            bt_d = d["donor_bt_dist"]
            bt_dist = DiscreteDist(bt_d["bts"], bt_d["probs"])
            bt = bt_dist.sample1()               
            weight = weight_dists.draw_weight(sex, age)

            donor = Donor(sex, age, bt, weight, d_id, patient)
            donors.append(donor)
            d_id += 1

    return donors

def generate(d, n_patients):
    weight_dists = WeightDistributions(d["weight_dists"])

    patients = generate_patients(d, n_patients, weight_dists)
    donors = generate_donors(d, patients, weight_dists)

    return patients, donors

if __name__=="__main__":
    n_patients = int(sys.argv[1])

    rand = random.random
    randint = random.randint

    with open('data/gen_data.json') as f:
        d = json.load(f)
        patients, donors = generate(d, n_patients)
        print len(patients), len(donors)
        for patient in patients:
            print patient
        for donor in donors:
            print donor
