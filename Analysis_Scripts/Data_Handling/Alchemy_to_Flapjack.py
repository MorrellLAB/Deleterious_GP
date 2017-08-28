#!/usr/bin/env python
"""Convert from ALCHEMY report to PLINK PED files. Translate the names from the
old names to the BOPAC names, and put them in physical map order. Check the
alleles and fix as necessary. Takes four arguments:
    1) Alchemy report
    2) BOPA Physical positions VCF
    3) Genetic map CSV
    4) SNP Name Translation CSV"""

import sys


def parse_vcf(vcf):
    """Read the VCF, store the SNPs in physical order."""
    vcf_data = {}
    bopa_alleles = {}
    with open(vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split('\t')
                snpid = tmp[2]
                ref = tmp[3]
                alt = tmp[4]
                pos = int(tmp[1])
                chrom = tmp[0]
                # Store just the alleles for checking flip errors
                bopa_alleles[snpid] = (ref, alt)
                if chrom not in vcf_data:
                    vcf_data[chrom] = [(snpid, pos, ref, alt)]
                else:
                    vcf_data[chrom].append((snpid, pos, ref, alt))
    # Sort it by position
    for c in vcf_data:
        vcf_data[c] = sorted(vcf_data[c], key=lambda x: x[1])
    return (bopa_alleles, vcf_data)


def parse_translation(name_trans):
    """Store the translation table as a dictionary."""
    trans = {}
    with open(name_trans, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split(',')
                trans[tmp[3]] = tmp[2]
    return trans


def parse_alchemy(trans, alchemy):
    """Parse the ALCHEMY report, and store the genotyping data in a dictionary
    of dictionaries. The first key is the individual ID, and the second key is
    the SNP ID, translated to BOPAC."""
    missing = '-'
    thresh = 0.7
    alc_data = {}
    with open(alchemy, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split('\t')
                if tmp[0] in trans:
                    snpid = trans[tmp[0]]
                else:
                    snpid = tmp[0]
                sample = tmp[1]
                ab_call = tmp[2]
                nuc_call = tmp[3]
                prob = float(tmp[4])
                # Check the sample ID and the probability. If the sample ID is
                # 'Blank', we skip it
                if sample == 'Blank':
                    continue
                else:
                    # Fix the sample names so that they match up with other
                    # data sources. Cycle 1 names are MS10S30XXX, Cycle 2 names
                    # are MS11S2XXXX, and Cycle 3 names are MS12_XXXX.
                    if sample.startswith('G'):
                        sample = sample.replace('G10W', 'MS10S3')
                        sample = sample.replace('-', '-0')
                    elif sample.startswith('MS11'):
                        sample = sample
                    elif sample.startswith('MS12'):
                        sample = sample
                    # Then check the probability. If it is less than our
                    # threshold (0.7) we append missing. Else, we append the
                    # call
                    if prob < thresh:
                        if sample in alc_data:
                            alc_data[sample].update({snpid: missing*2})
                        else:
                            alc_data[sample] = {snpid: missing*2}
                    else:
                        #calls = (nuc_call[0], nuc_call[1])
                        calls = ab_call[0] + ab_call[1]
                        if sample in alc_data:
                            alc_data[sample].update({snpid: calls})
                        else:
                            alc_data[sample] = {snpid: calls}
    return alc_data


def order_snps(ped_data, p_map):
    """Generate the order that the SNPs should go in by iterating through the
    chromosomes of the physical map."""
    snp_order = []
    for chrom in sorted(p_map):
        for snp in p_map[chrom]:
            if snp[0] in ped_data[ped_data.keys()[0]]:
                snp_order.append(snp[0])
    return snp_order


def generate_ped(g_mat, map_order):
    """Generate the PED data with the genotype matrix and the map order of
    the SNPs."""
    ped_file = []
    # Append the header
    ped_file.append([' # fjFile = GENOTYPE'])
    # Append the names of the SNPs
    ped_file.append([''] + map_order)
    for sample in sorted(g_mat):
        # Use missing codes for the family ID and the maternal and paternal IDs
        # these will be filled in later.
        geno = []
        for marker in map_order:
            geno.append(g_mat[sample][marker])
        ped_file.append([sample] + geno)
    return ped_file


def main(alchemy, bopa, genetic_map, name_trans):
    """Main function."""
    bopa_alleles, phys_map = parse_vcf(bopa)
    name_key = parse_translation(name_trans)
    alchemy_calls = parse_alchemy(name_key, alchemy)
    ordered = order_snps(alchemy_calls, phys_map)
    flapjack_file = generate_ped(alchemy_calls, ordered)
    # Then write out the data
    for row in flapjack_file:
        print '\t'.join(row)
    return


main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
