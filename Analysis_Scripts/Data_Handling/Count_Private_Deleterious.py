#!/usr/bin/env python
"""Counts the number of private deleterious SNPs in each of the parental
lines. Takes two arguments:
    1) The derived allele counts file
    2) The deleterious IDs file.
"""

import sys

daf = sys.argv[1]
del_ids = sys.argv[2]

deleterious = []
with open(del_ids, 'r') as f:
    for line in f:
        deleterious.append(line.strip())

private = {}
with open(daf, 'r') as f:
    for line in f:
        tmp = line.strip().split()
        snpid = tmp[0]
        if snpid not in deleterious or len(tmp) != 4:
            continue
        else:
            if len(tmp[3].split(',')) != 1:
                continue
            else:
                if tmp[3] == '-':
                    continue
                else:
                    if tmp[3] not in private:
                        private[tmp[3]] = 1
                    else:
                        private[tmp[3]] += 1

for lineid in sorted(private):
    print lineid, private[lineid]
