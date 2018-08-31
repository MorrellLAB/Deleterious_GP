#!/usr/bin/env python
"""Reorder a PLINK PED file so that it has markers listed in genetic map order,
rather than in physical map order. This is necessary because the barley genome
currently does not have physical position for SNPs.

Usage:
Fix_Map_Order.py Genetic_Order.map Physical_Order.map Physical_PED

will reorder a PED according to the order listed in Genetic_Order.map. The SNPs
listed in the Genetic_Order.map and Physical_Order.map must be the same, just
with a different sorting. Assumes the PED file is sorted by the order in
Physical_Order.map.

Also assumes that the PED is "full" = it has columns for Family ID, Individual
ID, Maternal ID, Paternal ID, Sex, and Phenotype."""

import sys


def reorder_map(orig_map, new_map):
    """Returns a list of integers that is the the index of each SNP in new_map
    in orig_map."""
    #   First, do a check to make sure that the maps have the same SNPs.
    same_snp = [n in orig_map for n in new_map]
    if len(orig_map) != len(new_map) and not all(same_snp):
        print("Error! The maps do not have the same SNPs!")
        exit(1)
    #   Then determine the order
    reorder = [orig_map.index(s) for s in new_map]
    return reorder


def reorder_ped(fixed_order, ped):
    """Prints a PED markers in genetic map order."""
    fixed_ped = []
    with open(ped, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            #   Strip out the genetic data
            gen = tmp[6:]
            #   Then, reorder the markers
            new_order = [gen[2*i] + '\t' + gen[2*i+1] for i in fixed_order]
            #   And append it all to the fixed PED output
            fixed_ped.append(tmp[0:6] + new_order)
    return fixed_ped


def usage():
    "Print a usage message."
    usage = """
Usage:
Fix_Map_Order.py Genetic_Order.map Physical_Order.map Physical_PED

will reorder a PED according to the order listed in Genetic_Order.map. The SNPs
listed in the Genetic_Order.map and Physical_Order.map must be the same, just
with a different sorting. Assumes the PED file is sorted by the order in
Physical_Order.map.

Also assumes that the PED is "full" = it has columns for Family ID, Individual
ID, Maternal ID, Paternal ID, Sex, and Phenotype."""
    print(usage)
    return


def main(orig_map, new_map, ped):
    """Main function"""
    #   Read the maps and get the SNP order from them
    orig_order = [line.strip().split('\t')[1]
                  for line
                  in open(orig_map, 'r')]
    new_order = [line.strip().split('\t')[1]
                 for line
                 in open(new_map, 'r')]
    fixed_order = reorder_map(orig_order, new_order)
    fixed_ped = reorder_ped(fixed_order, ped)
    for l in fixed_ped:
        print('\t'.join(l))
    return

if len(sys.argv) != 4:
    usage()
    exit(2)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3])
