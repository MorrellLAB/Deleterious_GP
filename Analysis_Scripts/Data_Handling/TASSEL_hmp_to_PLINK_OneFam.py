#!/usr/bin/env python
"""Convert from the HMP format from FSFHap (TASSEL) to PLINK PED/MAP. This
script is different than the other HapMap to PLINK because FSFHap outputs one
letter genotypes. This script is written for single families that will be
checked for Mendel errors. Takes two arguments:
    1) FSFHap (TASSEL) hmp file for one family
    2) PLINK MAP
 """

import sys

# Define a dictionary that links one-letter genotype calls to their
# two-character meanings.
IUPAC = {
    'R': ('A', 'G'),
    'Y': ('C', 'T'),
    'S': ('G', 'C'),
    'W': ('A', 'T'),
    'K': ('G', 'T'),
    'M': ('A', 'C'),
    'A': ('A', 'A'),
    'T': ('T', 'T'),
    'C': ('C', 'C'),
    'G': ('G', 'G'),
    'N': ('0', '0')
    }

def parse_hapmap(hmp):
    """Parse the HapMap file and store it as a nested dictionary. The first key
    is the sample name and the second key is the SNP name."""
    hap_dat = {}
    with open(hmp, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                samples = line.strip().split('\t')[11:]
                hap_dat['Samples'] = samples
                for s in samples:
                    hap_dat[s] = {}
            else:
                tmp = line.strip().split('\t')
                snpid = tmp[0]
                for samp, geno in zip(samples, tmp[11:]):
                    call = IUPAC[geno]
                    hap_dat[samp].update({snpid: (call[0], call[1])})
    return hap_dat


def generate_ped(ped_data, maporder, fam):
    """Generate the PED in the proper order."""
    ped = []
    # How many parents are in the dataset?
    npar = len(ped_data['Samples']) - 24
    # This is a dumb workaround - PLINK needs sexes for the parents
    sexes = ['1', '2']
    # Iterate over all samples. If it's a parent, then treat it differently
    for index, sample in enumerate(ped_data['Samples']):
        if index+1 <= npar:
            pid = '0'
            mid = '0'
            sex = sexes[index]
        else:
            pid = ped_data['Samples'][0]
            mid = ped_data['Samples'][1]
            sex = '0'
        fid = fam
        pheno = '-9'
        geno = []
        for snp in maporder:
            geno += list(ped_data[sample][snp])
        ped.append([fid, sample, pid, mid, sex, pheno] + geno)
    return ped


def main(hmp, plinkmap):
    """Main function."""
    hap = parse_hapmap(hmp)
    # Get the family name from the filename
    fam = hmp.replace('imputed_genotypes_Chr1H_', '').replace('.hmp.txt', '')
    # Get the map order
    order = []
    with open(plinkmap, 'r') as f:
        for line in f:
            order.append(line.strip().split()[1])
    pedfile = generate_ped(hap, order, fam)
    for row in pedfile:
        sys.stdout.write('\t'.join(row) + '\n')
    return


main(sys.argv[1], sys.argv[2])
