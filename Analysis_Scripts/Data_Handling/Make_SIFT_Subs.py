#!/usr/bin/env python
"""Use the output from SNP_Effect_Predictor.py to generate an input
substitutions file for SIFT."""

import sys

effects = {}
with open(sys.argv[1], 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split()
            ns = tmp[3]
            if ns == 'No':
                eff = tmp[8] + tmp[10] + tmp[9]
                if '*' in eff:
                    continue
                if tmp[4] in effects:
                    effects[tmp[4]].append(eff)
                else:
                    effects[tmp[4]] = [eff]

for tid in effects.keys():
    handle = open(tid + '.subs', 'w')
    handle.write('\n'.join(effects[tid]))
    handle.close()
