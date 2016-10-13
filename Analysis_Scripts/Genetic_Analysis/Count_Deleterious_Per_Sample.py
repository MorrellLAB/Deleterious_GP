#!/usr/bin/env python
"""Counts the number of deleterious SNPs per progeny line, outputs a two column
table. Takes three arguments:
    1) A file with deleterious SNP IDs, one per line, no header.
    2) Ancestral state assignments table - 4 columns
    3) VCF of population"""

import sys

#   Store the deleterious SNPs in a list
deleterious = []
with open(sys.argv[1], 'r') as f:
    for line in f:
        deleterious.append(line.strip())

#   Then, get the ancestral state for the deleterious SNPs.
anc = {}
with open(sys.argv[2], 'r') as f:
    for line in f:
        tmp = line.strip().split('\t')
        #   If the SNP ID (first column) is in the list of deleterious SNP IDs,
        #   save the ancestral state (second column). Note that if the
        #   ancestral state is N, then we skip it.
        if tmp[0] in deleterious and tmp[1] != 'N':
            anc[tmp[0]] = tmp[1]

#   Then, read through the VCF, and save the genotype calls for each line.
genotypes = {}
with open(sys.argv[3], 'r') as f:
    for line in f:
        #   Skip the header lines, which start with ##
        if line.startswith('##'):
            continue
        elif line.startswith('#CHROM'):
            #   VCF sample info line starts with #CHROM, and sample names start
            #   on the 10th column.
            samples = line.strip().split()[9:]
            for s in samples:
                genotypes[s] = 0
        else:
            tmp = line.strip().split()
            snpid = tmp[2]
            #   If we do not have good ancestral state information for this
            #   variant, then we exclude it
            if snpid not in anc:
                continue
            #   Otherwise, we append the calls to the matrix
            else:
                ref = tmp[3]
                alt = tmp[4]
                #   Set which is the derived allele. Reference allele is 0 in
                #   a VCF, and alternate is 1.
                if ref == anc[snpid]:
                    der = '1'
                else:
                    der = '0'
                sample_info = tmp[9:]
                gt = [g.split(':')[0] for g in sample_info]
                #   Iterate through the sample calls and count up how many are
                #   deleterious.
                for index, s in enumerate(samples):
                    genotypes[s] += gt[index].count(der)

#   Then, for each sample, print out the number of deleterious SNPs it carries.
#   Print out the cycle the line is part of, too.
for s in samples:
    if s.startswith('MS10S3'):
        cycle = 'C1'
    elif s.startswith('MS11S2'):
        cycle = 'C2'
    elif s.startswith('MS12_'):
        cycle = 'C3'
    else:
        cycle = 'C0'
    print '\t'.join([cycle, s, str(genotypes[s])])
