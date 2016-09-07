#!/usr/bin/env python
"""Convert from the HMP format from T3 into PLINK PED files."""

import sys

#   The missing value to use in the output file. PLINK uses 0 as missing
#   genotypes, and -9 as missing phenotypes
missing = '0'
pheno_missing = '-9'
samplenames = []
genotype_matrix = {}

with open(sys.argv[1], 'r') as f:
    for index, line in enumerate(f):
        tmp = line.strip().split('\t')
        if index == 0:
            #   In the first row, get the sample names.
            samplenames = tmp[11:]
        else:
            gt = []
            #   Convert to the tab-separated PLINK format
            for g in tmp[11:]:
                if g == 'NN':
                    gt.append(missing + '\t' + missing)
                else:
                    gt.append(g[0] + '\t' + g[1])
            #   Tack the genotypes onto our matrix
            genotype_matrix[tmp[0]] = gt

#   Which order should the SNPs be printed in? This should be from the PLINK
#   .map file
snporder = []
with open(sys.argv[2], 'r') as f:
    for line in f:
        snporder.append(line.strip())

ordered_geno = []
for s in snporder:
    if s in genotype_matrix:
        ordered_geno.append(genotype_matrix[s])
    else:
        continue

#   Print the .PED file, with missing values for family ID, maternal ID,
#   paternal ID, and sex. zip() will transpose a matrix, now indiviuals are
#   rows
for l in zip(samplenames, *ordered_geno):
    towrite = [
        missing,
        l[0],
        missing,
        missing,
        missing,
        pheno_missing,
        '\t'.join(l[1:])
        ]
    #   Print it out
    print '\t'.join(towrite)
