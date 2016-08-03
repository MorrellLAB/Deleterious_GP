#!/usr/bin/env python
"""Associate the PROVEAN score with the effects table."""

import sys
import os
import math

effects_table = sys.argv[1]
pred_dir = sys.argv[2]

#   Get the prediction files out of the predictions directory
pred_dir_contents = os.listdir(pred_dir)
pred_files = [
    fname
    for fname
    in pred_dir_contents
    if fname.endswith('_prov_out.txt')]
    #   Note that _prov_out.txt is my file extension. When others run PROVEAN
    #   they will have different extensions

provean_predictions = {}
#   Then, iterate through the files and open them up, and save the prediction
#   information in them. We use the transcript ID and cds position as keys
for f in pred_files:
    txid = '.'.join(f.rstrip('_prov_out.txt').rsplit('_', 1))
    if txid not in provean_predictions:
        provean_predictions[txid] = {}
    with open(os.path.join(pred_dir, f), 'r') as g:
        for line in g:
            if not line.startswith('#') or not line.startswith('['):
                tmp = line.strip().split()
                cdspos = tmp[0][1:-1]
                score = tmp[1]
                if cdspos not in provean_predictions[txid]:
                    provean_predictions[txid][cdspos] = score


#   Now associate the scores with the effects table
with open(effects_table, 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            print line.strip() + '\t' + 'PROVEAN'
        else:
            tmp = line.strip().split()
            silent = tmp[3]
            if silent == 'Yes':
                print line.strip() + '\t' + 'NA'
            else:
                txid = tmp[4]
                cdspos = tmp[10]
                if txid not in provean_predictions:
                    print line.strip() + '\t' + 'NA'
                elif cdspos not in provean_predictions[txid]:
                    print line.strip() + '\t' + 'NA'
                else:
                    print line.strip() + '\t' + provean_predictions[txid][cdspos]
