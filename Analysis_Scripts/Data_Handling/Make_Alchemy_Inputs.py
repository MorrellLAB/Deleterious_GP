#!/usr/bin/env python
"""Generate the required SNP map and sample map files for ALCHEMY input. Be sure
to set the flags/constants at the head of the script to properly generate the
files. Takes three arguments:
    1) Alchemy input intensities file
    2) SNP A-B mapping CSV file
    3) Prefix for output file name
Note that files will be written into the currect directory."""

import sys
import os


# Define constants here
#   Alternate sample ID for sample map
ALT_ID = 'FAKE'
#   Inbreeding coefficient. We are working with F3 individuals.
INB_COEFF = '0.75'
#   Dummy values for chromosome and position in SNP map
DUMMY = '1'


def get_samples(intensities):
    """Parse the ALCHEMY input intensities file to extract the list of samples
    that will be used."""
    samples = []
    with open(intensities, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split()
                if tmp[0] in samples:
                    continue
                else:
                    samples.append(tmp[0])
    return samples


def translate_snps(tl_table):
    """Generate a dictionary of name translations from SNP name to another. We
    will translate all SNPs to their 'SNP_Name', because that is what is present
    in the Genome Studio export"""
    t = {}
    with open(tl_table, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split(',')
                t[tmp[2]] = tmp[3]
    return t


def store_ab(ab_map):
    """Store the A/B allele calls in a dictionary. We no longer need to use
    the translation table, because we will convert all GenomeStudio reports into
    BOPAC names."""
    ab = {}
    with open(ab_map, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split(',')
                ab[tmp[0]] = (tmp[1], tmp[2])
    return ab


def generate_inputs(samples, ab_alleles, prefix):
    """Generate the input files."""
    # Replace slashes with dots to keep users from trying to write files all
    # over the place
    pref = prefix.replace('/', '.')
    snpmap_out = open(pref + '_snp_map.txt', 'w')
    sampmap_out = open(pref + '_samp_map.txt', 'w')
    # Generate the SNP map
    snpmap_out.write('#snp_id\tchromosome\tposition\tA_allele\tB_allele\n')
    for snp in sorted(ab_alleles):
        a = ab_alleles[snp][0]
        b = ab_alleles[snp][1]
        snpmap_out.write('\t'.join([snp, DUMMY, DUMMY, a, b]) + '\n')
    # Generate the sample map
    sampmap_out.write('#assay_id\tsample_id\talt_id\tF\n')
    for samp in samples:
        sampmap_out.write('\t'.join([samp, samp, ALT_ID, INB_COEFF]) + '\n')
    # Close the handles, print a little message
    snpmap_out.close()
    sys.stderr.write('Wrote ' + pref + '_snp_map.txt\n')
    sampmap_out.close()
    sys.stderr.write('Wrote ' + pref + '_samp_map.txt\n')
    return


def main(intensities, ab_map, prefix):
    """Main function."""
    # Get the samples
    samp = get_samples(intensities)
    # Store the A/B alleles
    ab_alleles = store_ab(ab_map)
    # Then, generate the inputs
    generate_inputs(samp, ab_alleles, prefix)
    return


if len(sys.argv) != 4:
    print """Generate the required SNP map and sample map files for ALCHEMY input. Be sure
to set the flags/constants at the head of the script to properly generate the
files. Takes three arguments:
    1) Alchemy input intensities file
    2) SNP A-B mapping CSV file
    3) Prefix for output file name"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3])
