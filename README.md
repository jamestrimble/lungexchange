# lungexchange

The concept of a _lung exchange_ scheme is similar to that of a kidney
exchange scheme, with a key difference being that each patient in a
lung exchange scheme requires two donors rather than one. To the
best of my knowledge, lung exchange has not yet been implemented in
practice, but the following recent papers discuss theoretical foundations:

- _Lung Exchange_. Haluk Ergin, Tayfun Sönmez, M. Utku Ünver. _Working
  paper_. [Link](http://www.tayfunsonmez.net/wp-content/uploads/2014/09/lung-exchange-6.pdf
)

- _Mechanism design and implementation for lung exchange_. Suiqian Luo and
  Pingzhong Tang. *IJCAI-2015*, Buenos Aires, Argentina. [Link](http://iiis.tsinghua.edu.cn/~kenshin/lung.pdf)

The theory of lung exchange was also discussed in a recent
[_Economist_ article](http://www.economist.com/blogs/freeexchange/2014/09/lung-exchanges).

This repository contains my work so far on implementing algorithms for
lung exchange. In particular, it includes integer programming formulations
for:

- Optimisation with unbounded exchange size, using an edge formulation
- Optimisation with bounded exchange size, using a variant of Constantino et al.'s
  [extended edge formulation](http://www.sciencedirect.com/science/article/pii/S0377221713004244) for kidney exchange
- Optimisation with bounded exchange size, with a variable for each feasible
  exchange.

## Requirements

[Google or-tools](https://developers.google.com/optimization/installing) is required for solving integer programs.

## Notes on usage

### Generator

To generate a lung-exchange instance with 2 patients and 4 donors, run the following command.

    python lungexchange/lung_generator.py 2

I have tried to closely imitate the generator described by Ergin et al. The generator outputs an instance like this:

    2 4
    CF M 40 O 236
    CF F 18 A 177
    0 F 62 O 185
    0 M 42 O 187
    1 F 27 A 139
    1 F 27 A 145

The first row gives (number of patients), (number of donors). The patients follow, with each row showing condition (cystic fibrosis or pulmonary hypertension), gender, age, blood type, and weight. The remaining rows are donors in the format paired-patient-ID, gender, age, blood type, weight.

### Optimisation

To find the number of transplants in an optimal solution, run:

    python lungexchange/optimise.py instance.txt 5 True

The command-line arguments are instance file name, maximum permitted cycle size, and whether weights should be taken into account whether a donor is compatible with a patient. The results look like:

    instance.txt    1   50  48.0    18.0    28.0    34.0    40.0    5

The first four columns are (1) file name, (2) an indictor for whether donor and patient weights were used to determine compatibility, (3) number of patients in the instance, and (4) number of donors used in an optimal solution with unrestricted exchange size. The final column is the number of patients in the instance who have two compatible paired donors. The remaining columns are the number of donors used in an optimal solution with an exchange size limit of 2, 3, ..., n. Note that if a patient's two donors are able to give directly to the patient, these both count towards the optimal value.

More verbose output can be viewed using:
    
    python lungexchange/optimise_verbose.py instance.txt 4

Note that if a patient P is compatible with two of P's own paired donors, then the program does not identify compatibilities between P and any other patient's paired donors, or between P's paired donors and any other patient.
