import pandas as pd
import sys

def start():
    if len(sys.argv) < 2:
        print "Usage: python " + sys.argv[0] + " results_file.txt"
        return

    d = pd.read_table(sys.argv[1], header=None, names=['file',
        'use_wt', 'n_patients', 'optu', 'opt2', 'opt3', 'opt4',
        'opt5', 'direct'])

    d = d[['file', 'use_wt', 'n_patients', 'direct', 'opt2',
        'opt3', 'opt4', 'opt5', 'optu']]

    d['optu'] = d['optu'] / 2 - d['direct']
    d['opt2'] = d['opt2'] / 2 - d['direct']
    d['opt3'] = d['opt3'] / 2 - d['direct']
    d['opt4'] = d['opt4'] / 2 - d['direct']
    d['opt5'] = d['opt5'] / 2 - d['direct']

    print "MEANS"
    print "-----"
    print d.groupby(['n_patients', 'use_wt']).mean()

    print
    print "STANDARD DEVIATIONS"
    print "-------------------"
    print d.groupby(['n_patients', 'use_wt']).std()

if __name__=='__main__':
    start()
