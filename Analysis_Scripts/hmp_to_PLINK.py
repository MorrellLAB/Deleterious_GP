#!/usr/bin/env python
"""Convert from the hapmap(?) format into PLINK PED format. In the original
file, SNPs are rows and individuals are columns, and genotype calls are coded
as -1, 0, 1."""

import sys

#   The missing value to use in the output file. PLINK uses 0 as missing
#   genotypes, and -9 as missing phenotypes
missing = '0'
pheno_missing = '-9'
samplenames = []
genotype_matrix = []

with open(sys.argv[1], 'r') as f:
    for index, line in enumerate(f):
        tmp = line.strip().split('\t')
        if index == 0:
            #   In the first row, get the sample names.
            samplenames = tmp[4:]
        else:
            #   Alleles are the second column, separated by a slash
            alleles = tmp[1].split('/')
            gt = []
            #   Convert from the -1/0/1 genotype calls to ATCG calls for PLINK
            for g in tmp[4:]:
                if g == '-1':
                    #   Homozygous for the first allele
                    gt.append(alleles[0] + '\t' + alleles[0])
                elif g == '0':
                    #   Heterozygous
                    gt.append(alleles[0] + '\t' + alleles[1])
                elif g == '1':
                    #   Homozygous for the second allele
                    gt.append(alleles[1] + '\t' + alleles[1])
                elif g == 'NA':
                    #   Missing
                    gt.append(missing + '\t' + missing)
            #   Tack the genotypes onto our matrix
            genotype_matrix.append(gt)

#   Print the .PED file, with missing values for family ID, maternal ID,
#   paternal ID, and sex. zip() will transpose a matrix, now indiviuals are
#   rows
for l in zip(samplenames, *genotype_matrix):
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
