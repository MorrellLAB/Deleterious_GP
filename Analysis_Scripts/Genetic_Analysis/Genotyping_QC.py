#!/usr/bin/env python
"""Set the Mendel errors to missing in the dataset. We do this on an individual
genotype basis. Takes three arguments:
    1) Mendel Error summary
    2) PED file
    3) MAP file"""

import sys
import pprint


def plink(p, m):
    """Parse the PED/MAP files and return a complicated dictionary structure
    that has named SNP genotypes, and the marker order."""
    ped_data = {}
    snp_order = []
    ped_meta = []
    # Iterate through the MAP file and save the SNP order
    with open(m, 'r') as f:
        for line in f:
            snp_order.append(line.strip().split()[1])
    # Then, iterate through the PED and keep track of the genotypes
    with open(p, 'r') as f:
        for line in f:
            # Save the PED metadata. We need this to produce a valid PED file
            # later
            tmp = line.strip().split()
            ped_meta.append(tmp[0:6])
            # Then, append the genotypes into the ped_data dict
            ped_data[tmp[1]] = {}
            for index, snp in enumerate(snp_order):
                genotype = (tmp[6 + (2*index)], tmp[6 + (2*index + 1)])
                ped_data[tmp[1]][snp] = genotype
    return (ped_data, snp_order, ped_meta)


def set_me_missing(mendel, ped):
    """Iterate through the mendel errors file and set the proper genotypes to
    missing."""
    ped_me = ped
    with open(mendel, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            # Skip the header lines. Some are interspersed throughout the file,
            # so skip them all.
            if tmp[0] == 'FID':
                continue
            else:
                indid = tmp[1]
                snpid = tmp[3]
                # If the individual ID is not in the PED data, skip it; it has
                # already been filtered
                if indid not in ped_me:
                    continue
                # If snpid is not in the PED data, skip it, because it may have
                # been dropped already
                if snpid not in ped_me[indid]:
                    continue
                else:
                    ped_me[indid][snpid] = ('0', '0')
    return ped_me


def print_ped(ped_dat, ped_meta, order):
    """Print a PED file."""
    # First, iterate through the ped meta data
    for l in ped_meta:
        toprint = l
        # Then, get the individual ID. It's the second field
        iid = l[1]
        # Iterate through the genotype data, in order, and append the calls
        for snp in order:
            toprint += list(ped_dat[iid][snp])
        print ' '.join(toprint)
    return


def main(mendel, pedfile, mapfile):
    """Main function."""
    pdat, order, pmeta = plink(pedfile, mapfile)
    flt_ped = set_me_missing(mendel, pdat)
    print_ped(flt_ped, pmeta, order)
    return


if len(sys.argv) != 4:
    print """Set the Mendel errors to missing in the dataset. We do this on an individual
genotype basis. Takes three arguments:
    1) Mendel Error summary
    2) PED file
    3) MAP file"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3])
