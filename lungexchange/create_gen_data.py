import pandas as pd
import json

d = {
  "prob_cystic_fibrosis": 14.57/(14.57+5.47),

  "cystic_f_dists": {
    "prob_female": 0.4874,
    "age_dist_f": {
      "probs": [0.7624, 0.2078, 0.0298],
      "age_bands": [
        {"min": 18, "max": 34},
        {"min": 35, "max": 49},
        {"min": 50, "max": 64}
      ]
    },
    "age_dist_m": {
      "probs": [0.7030, 0.2597, 0.0373],
      "age_bands": [
        {"min": 18, "max": 34},
        {"min": 35, "max": 49},
        {"min": 50, "max": 64}
      ]
    },
    "bt_dist": {
      "probs": [0.4535, 0.4178, 0.0977, 0.0314],
      "bts": ["O", "A", "B", "AB"]
    }
  },

  "pulm_h_dists": {
    "prob_female": 0.7618,
    "age_dist_f": {
      "probs": [0.3030, 0.4376, 0.2595],
      "age_bands": [
        {"min": 18, "max": 34},
        {"min": 35, "max": 49},
        {"min": 50, "max": 64}
      ]
    },
    "age_dist_m": {
      "probs": [0.2938, 0.4613, 0.2448],
      "age_bands": [
        {"min": 18, "max": 34},
        {"min": 35, "max": 49},
        {"min": 50, "max": 64}
      ]
    },
    "bt_dist": {
      "probs": [0.4696, 0.3669, 0.1249, 0.0397],
      "bts": ["O", "A", "B", "AB"]
    }
  },

  "donor_bt_dist": {
    "probs": [0.44, 0.42, 0.1, 0.04],
    "bts": ["O", "A", "B", "AB"]
  },
}

age_df = pd.read_csv("data/ages.csv")
weight_df = pd.read_csv("data/weights.csv")

ages = list(age_df.Age)
male_probs = list(age_df.M)
female_probs = list(age_df.F)

d["donor_age_dist_f"] = {
  "probs": female_probs,
  "ages": ages
}

d["donor_age_dist_m"] = {
  "probs": male_probs,
  "ages": ages
}

weight_dists = []
weight_bands = [{"min": min_wt, "max": max_wt} for min_wt, max_wt
                  in zip(weight_df.min_wt, weight_df.max_wt)]
for col_name in list(weight_df):
    if len(col_name)==3:
        sex = col_name[0]
        max_age = int(col_name[1:])
        cum_probs = list(weight_df[col_name])
        weight_dists.append({
            "sex": sex,
            "max_age": max_age,
            "cum_probs": cum_probs,
            "weight_bands": weight_bands
        })

d["weight_dists"] = weight_dists

print json.dumps(d, sort_keys=True,
                  indent=4, separators=(',', ': '))
