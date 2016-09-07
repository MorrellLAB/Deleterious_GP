#!/usr/bin/env python
"""Apply stringent filtering criteria to a VCF to make a set of high confidence
variants for recalibration. Note that this is not used to filter to a set of
analysis-ready variants. This is only for the GATK pipeline."""

import sys

#   Define the filters here
#       Variant quality >= 60
#       Depth is [105, 3150] (5 to 150 reads per sample, on average)
#       At most 2 missing calls
#       At most two heterozygous calls
#       Phred-scaled p-value for excess heterozygosity > 13 (p < 0.05)
#       No length polymorphisms
minqual = 60
mindepth = 105
maxdepth = 3150
maxmissing = 2
maxhet = 2
excesshet = 13

with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.startswith('#'):
            print line.strip()
        else:
            tmp = line.strip().split()
            #   Get the genotype calls
            genotypes = [g.split(':')[0] for g in tmp[9:]]
            #   Count the number of missing calls and heterozygous calls
            nmiss = genotypes.count('./.')
            nhet = genotypes.count('0/1') + genotypes.count('1/0')
            #   And get the depth and Fisher's Exact Test for excess het.
            for field in tmp[7].split(';'):
                if field.startswith('DP'):
                    dp = int(field[3:])
                if field.startswith('ExcessHet'):
                    ehet = float(field[10:])
            #   If REF or ALT are not 1 character long, it is a length
            #   polymorphism. Skip it.
            if len(tmp[3]) != 1 or len(tmp[4]) != 1:
                continue
            #   Skip low quality variants
            elif float(tmp[5]) < minqual:
                continue
            #   Skip low-depth variants
            elif dp < mindepth or dp > maxdepth:
                continue
            #   Skip variants with a lot of missing data
            elif nmiss > maxmissing:
                continue
            #   Skip variants with a lot of heterozygosity
            elif nhet > maxhet:
                continue
            elif ehet > excesshet:
                continue
            else:
                print line.strip()
