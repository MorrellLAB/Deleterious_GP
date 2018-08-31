#!/usr/bin/env python
"""Calculate the derived site frequency spectra for the deleterious variants in
the three cycles of the genomic prediction and selection experiment. Takes four
arguments:
    1) Ancestral state assignment
    2) C1 Frequencies
    3) C2 Frequencies
    4) C3 Frequencies
"""

import sys
import gzip


def store_ancestral(a):
    """Parse the gzipped ancestal alleles file and store it in a dictionary."""
    anc_alleles = {}
    snp_order = []
    with gzip.open(a, 'rt') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split()
                snpid = tmp[2]
                ancestral = tmp[3]
                derived = tmp[4]
                anc_alleles[snpid] = ancestral
                snp_order.append(snpid)
    return (anc_alleles, snp_order)


def calc_daf(a, f):
    """Use the ancestral state to calculate the derived allele frequency for
    each variant. If the ancestral state is N, then we will return NA for the
    derived allele frequency."""
    dafs = {}
    with gzip.open(f, 'rt') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split()
                # First check if the SNP ID is not in the ancestral states. If
                # not, then we toss it
                if tmp[1] not in a:
                    continue
                else:
                    snpid = tmp[1]
                    a1 = tmp[2]
                    a2 = tmp[3]
                    na1 = 2*int(tmp[4]) + int(tmp[5])
                    na2 = 2*int(tmp[6]) + int(tmp[5])
                    na = 2*(int(tmp[4]) + int(tmp[5]) + int(tmp[6]))
                    # get the ancestral state
                    astate = a[snpid]
                    if astate == 'N':
                        d = 'NA'
                    elif astate == a1:
                        d = str(float(na2)/float(na))
                    elif astate == a2:
                        d = str(float(na1)/float(na))
                    else:
                        d = 'NA'
                    dafs[snpid] = d
    return dafs


def main(anc, c1, c2, c3):
    """Main function."""
    # store the ancestral alleles
    a_a, snp_order = store_ancestral(anc)
    # Then, parse the frequency files and get the derived allele frequencies
    c1_daf = calc_daf(a_a, c1)
    c2_daf = calc_daf(a_a, c2)
    c3_daf = calc_daf(a_a, c3)
    # Print a header
    print('SNP_ID\tAnc\tC1_DAF\tC2_DAF\tC3_DAF')
    for snp in snp_order:
        print('\t'.join([
            snp,
            a_a.get(snp, 'N'),
            c1_daf.get(snp, 'NA'),
            c2_daf.get(snp, 'NA'),
            c3_daf.get(snp, 'NA')]))
    return


if len(sys.argv) != 5:
    print("""Calculate the derived site frequency spectra for the deleterious variants in
the three cycles of the genomic prediction and selection experiment. Takes four
arguments:
    1) Ancestral state assignment
    2) C1 Frequencies
    3) C2 Frequencies
    4) C3 Frequencies""")
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
