#!/usr/bin/env python
"""Make the BOPA SNPs consistent between a reference VCF and a PLINK
pedigree file. Will check for strand errors, but cannot handle actual
genotyping errors. You're on your own for that... Takes three arguments:
    1) Reference VCF with true allele states
    2) PED to fix
    3) MAP with marker order. Names must match the VCF.
"""

import sys
import pprint


# define a dictionary for reverse complement lookups
REVCOMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def store_alleles(vcf):
    """Reads a VCF and stores a dictionary of the reference and alternate
    alleles for a SNP. SNP IDs are keys, and tuples of (ref, alt) are values."""
    states = {}
    with open(vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split()
                states[tmp[2]] = (tmp[3], tmp[4])
    return states


def parse_map(map_file):
    """Parse the map file and store the order in which the SNPs occur. We just
    use this to know which order to do the lookup for later, when checking the
    alleles in the PED file. """
    map_order = []
    with open(map_file, 'r') as f:
        for line in f:
            tmp = line.strip().split('\t')
            map_order.append(tmp[1])
    return map_order


def order_snps(map_data, alleles):
    """Takes the SNP order from the PLINK MAP and the VCF data to make an
    ordered list of alleles for lookup."""
    ordered = []
    for snp in map_data:
        if snp in alleles:
            ordered.append(alleles[snp])
        else:
            ordered.append(('N', 'N'))
    return ordered


def t_ped(ped_file):
    """Transpose a PED file for easy iteration."""
    ped_to_trans = []
    metadata = []
    with open(ped_file, 'r') as f:
        for line in f:
            ped_to_trans.append(line.strip().split()[6:])
            metadata.append(line.strip().split()[:6])
    return metadata, zip(*ped_to_trans)


def fix_alleles(ped_column, states):
    """Check the states of the SNPs passed, and fix them accordingly. If they
    cannot be resolved, set them to 0 (missing in PED)"""
    fixed_states = []
    for call in zip(*ped_column):
        a1 = call[0]
        a2 = call[1]
        if a1 in states and a2 in states:
            fixed_states.append((a1, a2))
        elif a1 == '0' and a2 == '0':
            fixed_states.append(('0', '0'))
        elif REVCOMP[a1] in states and REVCOMP[a2] in states:
            fixed_states.append((REVCOMP[a1], REVCOMP[a2]))
        else:
            fixed_states.append(('0', '0'))
    return fixed_states


def main(vcf, ped, mapf):
    """Main function. Do work here."""
    v_alleles = store_alleles(vcf)
    map_snps = parse_map(mapf)
    snps_in_order = order_snps(map_snps, v_alleles)
    metadata, trans_ped = t_ped(ped)
    fixed_ped = []
    counter = 0
    for snp in snps_in_order:
        fixed_ped.append(fix_alleles(trans_ped[(2*counter):(2*counter)+2], snp))
        counter += 1
    # Transpose it again
    fixed_ped = zip(*fixed_ped)
    for index, pedrow in enumerate(metadata):
        print ' '.join(pedrow) + ' ' + ' '.join([' '.join(a) for a in fixed_ped[index]])


if len(sys.argv) != 4:
    print """
Make the BOPA SNPs consistent between a reference VCF and a PLINK
pedigree file. Will check for strand errors, but cannot handle actual
genotyping errors. You're on your own for that... Takes three arguments:
    1) Reference VCF with true allele states
    2) PED to fix
    3) MAP with marker order. Names must match the VCF."""
    exit(1)

main(sys.argv[1], sys.argv[2], sys.argv[3])
