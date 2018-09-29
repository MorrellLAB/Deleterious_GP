#!/usr/bin/env/python
"""Make a data frame for plotting that has only the random and selected lines,
which cycle they are from, which panel they are in, and the derived homozygous
counts of various functional classes of SNPs. Takes two arguments:
    1) Homozygous derived file (gzipped)
    2) Adjusted phenotypic data file with ran/sel/cycle information
"""

import sys
import gzip

# Make the sets used for membership tests
c0 = set()
c1_r = set()
c1_s = set()
c2_r = set()
c2_s = set()
c3_r = set()
c3_s = set()
with open(sys.argv[2], 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split(',')
            lid = tmp[0]
            cyc = tmp[1]
            panel = tmp[2]
            if cyc == 'C0':
                c0.add(lid)
            elif cyc == 'C1':
                if panel == 'ran':
                    c1_r.add(lid)
                elif panel == 'sel':
                    c1_s.add(lid)
                else:
                    continue
            elif cyc == 'C2':
                if panel == 'ran':
                    c2_r.add(lid)
                elif panel == 'sel':
                    c2_s.add(lid)
                else:
                    continue
            elif cyc == 'C3':
                if panel == 'ran':
                    c3_r.add(lid)
                elif panel == 'sel':
                    c3_s.add(lid)
                else:
                    continue
            else:
                continue

# Then, iterate through the burden file and slice out the lines that we are
# interested in keeping
with gzip.open(sys.argv[1], 'rt') as f:
    for index, line in enumerate(f):
        if index == 0:
            # Print a new header
            print(','.join([
                'Line', 'Cycle', 'Type', 'NC', 'SY', 'NS', 'DE'
                ]))
        else:
            tmp = line.strip().split()
            if tmp[0] in c0:
                cyc = 'C0'
                panel = 'NA'
            elif tmp[0] in c1_r:
                cyc = 'C1'
                panel = 'Ran'
            elif tmp[0] in c1_s:
                cyc = 'C1'
                panel = 'Sel'
            elif tmp[0] in c2_r:
                cyc = 'C2'
                panel = 'Ran'
            elif tmp[0] in c2_s:
                cyc = 'C2'
                panel = 'Sel'
            elif tmp[0] in c3_r:
                cyc = 'C3'
                panel = 'Ran'
            elif tmp[0] in c3_s:
                cyc = 'C3'
                panel = 'Sel'
            else:
                continue
            # Print the line
            toprint = [tmp[0], cyc, panel] + tmp[1:]
            print(','.join(toprint))
