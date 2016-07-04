#!/usr/bin/env python
"""Associate the BAD_Mutations prediction information with the effects table.
Relies on the transcript ID and the codon position to be correct."""

import sys

effect_table = sys.argv[1]
bad_mutations_table = sys.argv[2]

bad_mutations_data = {}
with open(bad_mutations_table, 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            header = line.strip().split()
        else:
            tmp = line.strip().split()
            #   Key on transcript ID and position. Since HyPhy chokes on . in
            #   gene names, we had to replace . with _. We fix that now
            txid = '.'.join(tmp[0].rsplit('_', 1))
            #   And get the CDS position
            cdspos = tmp[1]
            if txid not in bad_mutations_data:
                bad_mutations_data[txid] = {}
            if cdspos not in bad_mutations_data[txid]:
                bad_mutations_data[txid][cdspos] = tmp


#   Then, iterate through the effect table, printing out the proper info
with open(effect_table, 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            print line.strip() + '\t' + '\t'.join(header[2:])
        else:
            tmp = line.strip().split()
            silent = tmp[3]
            if silent == 'Yes':
                print line.strip() + '\t' + '\t'.join(['NA']*11)
            else:
                txid = tmp[4]
                cdspos = tmp[10]
                if txid not in bad_mutations_data:
                    print line.strip() + '\t' + '\t'.join(['NA']*11)
                elif cdspos not in bad_mutations_data[txid]:
                    print line.strip() + '\t' + '\t'.join(['NA']*11)
                else:
                    print line.strip() + '\t' + '\t'.join(bad_mutations_data[txid][cdspos][2:])
