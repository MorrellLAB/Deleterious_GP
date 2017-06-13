#!/usr/bin/env python
"""Convert from the HMP format from T3 into PLINK PED files. Will order markers
according to a supplied PLINK Map file, and check alleles for RC errors with
a VCF. Takes three agruments:
    1) T3 hmp file
    2) PLINK MAP
    3) BOPA VCF"""

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


def parse_hapmap(hmp):
    """Parse the HapMap file and store it as a nested dictionary. The first key
    is the sample name and the second key is the SNP name."""
    hap_dat = {}
    with open(hmp, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                samples = line.strip().split('\t')[11:]
                for s in samples:
                    hap_dat[s] = {}
            else:
                tmp = line.strip().split('\t')
                snpid = tmp[0]
                for samp, geno in zip(samples, tmp[11:]):
                    if geno == 'NN':
                        call = '00'
                    else:
                        call = geno
                    hap_dat[samp].update({snpid: (call[0], call[1])})
    return hap_dat


def fix_alleles(ped_data, bopa_data):
    """Check the alleles in the PED data against the BOPA alleles."""
    fixed_ped = {}
    # First, iterate through the samples
    for samp in ped_data:
        fixed_ped[samp] = {}
        # Then through the SNPs
        for snp in ped_data[samp]:
            alc_1 = ped_data[samp][snp][0]
            alc_2 = ped_data[samp][snp][1]
            # If the BOPA SNP does not have a physical position (and thus, no
            # reference and alternate allele calls), then we can't do anything
            # and just append it...
            if snp not in bopa_data:
                fixed_ped[samp].update({snp: (alc_1, alc_2)})
            elif alc_1 == '0' and alc_2 == '0':
                # If tehy are both missing, then we just append missing calls
                fixed_ped[samp].update({snp: ('0', '0')})
            elif alc_1 in bopa_data[snp] and alc_2 in bopa_data[snp]:
                # If the alleles are in the same orientation as the VCF data,
                # then we're all good.
                fixed_ped[samp].update({snp: (alc_1, alc_2)})
            elif REVCOMP[alc_1] == alc_2:
                # If the two alleles are reverse complements, i.e., A/T or C/G,
                # then we just keep them as-is, since those are hard to detect
                # as strand-flip errors.
                fixed_ped[samp].update({snp: (alc_1, alc_2)})
            elif REVCOMP[alc_1] in bopa_data[snp] and REVCOMP[alc_2] in bopa_data[snp]:
                # Else, if the ALCHEMY calls are RCed with respect to the
                # reference genome, we flip them
                fixed_ped[samp].update({snp: (REVCOMP[alc_1], REVCOMP[alc_2])})
            else:
                # If they cannot be resolved by any of these, we set them to
                # missing.
                fixed_ped[samp].update({snp: ('0', '0')})
    return fixed_ped


def generate_ped(ped_data, maporder):
    """Generate the PED in the proper order."""
    ped = []
    for sample in sorted(ped_data):
        fid = '0'
        pid = '0'
        mid = '0'
        sex = '0'
        pheno = '-9'
        geno = []
        for snp in maporder:
            geno += list(ped_data[sample][snp])
        ped.append([fid, sample, pid, mid, sex, pheno] + geno)
    return ped


def main(hmp, plinkmap, bopa_vcf):
    """Main function."""
    hap = parse_hapmap(hmp)
    bopa = store_alleles(bopa_vcf)
    fixed_hap = fix_alleles(hap, bopa)
    # Get the map order
    order = []
    with open(plinkmap, 'r') as f:
        for line in f:
            order.append(line.strip().split()[1])
    pedfile = generate_ped(fixed_hap, order)
    for row in pedfile:
        sys.stdout.write('\t'.join(row) + '\n')
    return


main(sys.argv[1], sys.argv[2], sys.argv[3])

# #   The missing value to use in the output file. PLINK uses 0 as missing
# #   genotypes, and -9 as missing phenotypes
# missing = '0'
# pheno_missing = '-9'
# samplenames = []
# genotype_matrix = {}

# with open(sys.argv[1], 'r') as f:
#     for index, line in enumerate(f):
#         tmp = line.strip().split('\t')
#         if index == 0:
#             #   In the first row, get the sample names.
#             samplenames = tmp[11:]
#         else:
#             gt = []
#             #   Convert to the tab-separated PLINK format
#             for g in tmp[11:]:
#                 if g == 'NN':
#                     gt.append(missing + '\t' + missing)
#                 else:
#                     gt.append(g[0] + '\t' + g[1])
#             #   Tack the genotypes onto our matrix
#             genotype_matrix[tmp[0]] = gt

# #   Which order should the SNPs be printed in? This should be from the PLINK
# #   .map file
# snporder = []
# with open(sys.argv[2], 'r') as f:
#     for line in f:
#         snporder.append(line.strip())

# ordered_geno = []
# for s in snporder:
#     if s in genotype_matrix:
#         ordered_geno.append(genotype_matrix[s])
#     else:
#         continue

# #   Print the .PED file, with missing values for family ID, maternal ID,
# #   paternal ID, and sex. zip() will transpose a matrix, now indiviuals are
# #   rows
# for l in zip(samplenames, *ordered_geno):
#     towrite = [
#         missing,
#         l[0],
#         missing,
#         missing,
#         missing,
#         pheno_missing,
#         '\t'.join(l[1:])
#         ]
#     #   Print it out
#     print '\t'.join(towrite)
