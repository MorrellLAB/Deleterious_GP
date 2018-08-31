#!/usr/bin/env python
"""Add the intervening inbreeding steps to the pedigree for AlphaPeel. This is
because we need to model the selfing of the progeny to the F3 generation. Takes
one argument:
    1) AlphaPeel pedigree without inbreeding
"""

import sys
import pprint


def parse_ped(p):
    """Parse the pedigree and return the family-level pedigree."""
    parents = []
    fam_ped = {}
    order = []
    progeny = {}
    with open(p, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            # For the fam-level pedigree, we just want the progeny family ID
            if tmp[1] != '0' and tmp[2] != '0':
                famid = tmp[0].split('-')[0]
                if famid in progeny:
                    progeny[famid].append(tmp[0])
                else:
                    progeny[famid] = [tmp[0]]
                if famid in fam_ped:
                    pass
                else:
                    order.append(famid)
                    fam_ped[famid] = (tmp[1], tmp[2])
            else:
                parents.append(tmp[0])
    return (parents, progeny, fam_ped, order)


def make_selfs(famped, progeny, order, famsize=24):
    """Make the bulk F1, and SSD F2 individuals that give rise to the genotyped
    individuals in the experiment. We pass the order of families just so that
    the numbers look nice."""
    full_ped = []
    for index, fam in enumerate(order):
        famno = 'Fam_' + str(index + 1).zfill(3)
        # First get the parents of the F1
        f1_p1 = famped[fam][0]
        f1_p2 = famped[fam][1]
        f1 = []
        for i in range(0, famsize):
            idno = str(i + 1).zfill(2)
            f1_id = famno + '_F1_I' + idno
            full_ped.append([f1_id, f1_p1, f1_p2])
            f1.append(f1_id)
        f2 = []
        for i in range(0, famsize):
            idno = str(i + 1).zfill(2)
            f2_id = famno + '_F2_I' + idno
            f1_par = famno + '_F1_I' + idno
            full_ped.append([f2_id, f1_par, f1_par])
            f2.append(f2_id)
        for f2par, f3prog in zip(f2, progeny[fam]):
            full_ped.append([f3prog, f2par, f2par])
    return full_ped


def main(noinb_ped):
    """Main function."""
    par, pro, fam, order = parse_ped(noinb_ped)
    full_ped = make_selfs(fam, pro, order)
    # Then print out the full pedigree
    for p in par:
        print(p, '0', '0')
    for l in full_ped:
        print(' '.join(l))
    return


if len(sys.argv) != 2:
    print("""Add the intervening inbreeding steps to the pedigree for AlphaPeel. This is
because we need to model the selfing of the progeny to the F3 generation. Takes
one argument:
    1) AlphaPeel pedigree without inbreeding""")
    exit(1)
else:
    main(sys.argv[1])
