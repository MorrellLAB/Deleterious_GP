#!/usr/bin/env python
"""Calculate the dosage of deleterious variants across all lines in the
experiment. Filters to those that only have ancestral state. Takes four
arguments:
    1) SNP IDs in AlphaPeel dosages file
    2) Inferred ancestal states for each SNP (gzipped)
    3) Deleterious SNP IDs
    4) Combined AlphaPeel dosages output (gzipped)
"""


import sys
import gzip
import pprint


def anc_del_ids(i, a, d):
    """Dumb function to intersect three lists. We want to return a list of SNPs
    that are deleterious and have inferred ancestral states. Further, we want
    to return the IDs and the column offsets in the dosages file."""
    all_ids = []
    anc_ids = []
    polarity = []
    del_ids = []
    # Store all IDs in a list
    with open(i, 'r') as f:
        for line in f:
            snpid = line.strip().split()[2]
            if ';' in snpid:
                snpid = snpid.split(';')[1]
            all_ids.append(snpid)
    # Then those with ancestral state
    with gzip.open(a, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split()
                snpid = tmp[2]
                anc = tmp[3]
                morex = tmp[5]
                if anc == 'N':
                    continue
                else:
                    anc_ids.append(snpid)
                    # If reference state is ancestral, then we append 0 - the
                    # deleterious dosage for these is d becaues the dosage
                    # reported is the dosage of alt allele.
                    if morex == 'A':
                        polarity.append(0)
                    elif morex == 'D':
                        polarity.append(1)
    # Then the deleterious IDs
    with open(d, 'r') as f:
        for line in f:
            del_ids.append(line.strip())
    # Next, we intersect the lists
    k_off = []
    k_pol = []
    d_flt = []
    for snpid in del_ids:
        # if the SNP does not have ancestral state, we skip it
        if snpid not in anc_ids or snpid not in all_ids:
            continue
        else:
            d_flt.append(snpid)
            # Get the offset of it from the main alphapeel file
            ap_off = all_ids.index(snpid)
            # And get its polarity
            pol = polarity[anc_ids.index(snpid)]
            # Append these to the list
            k_off.append(ap_off)
            k_pol.append(pol)
    # Return the package
    return (d_flt, k_off, k_pol)


def del_dosage(cols, pols, dos):
    """Iterate through the dosage matrix, keep the columns that are interesting
    and adjust the polarity accordingly."""
    mat = []
    samp_ids = []
    with gzip.open(dos, 'rb') as f:
        for line in f:
            tmp = line.strip().split()
            # This is a huge list. 400k ish entries.
            flt_row = []
            for index, field in enumerate(tmp):
                if index == 0:
                    # Keep the sample IDs.
                    samp_ids.append(field)
                elif index-1 in cols:
                    # If the field is in the list of columns that we wish to
                    # save, process it.
                    offs = cols.index(index-1)
                    p = pols[offs]
                    if p == 0:
                        # If the polarity is 0, then reference is ancestral,
                        # and we just append the dosage value
                        flt_row.append(field)
                    else:
                        # Else, morex is derived, and we have to use 2-d
                        flt_row.append(str(2.0-float(field)))
                else:
                    continue
            # Then append the filtered row to the matrix
            mat.append(flt_row)
    # Return the pair
    return (samp_ids, mat)


def main(ids, anc, dsnps, dosages):
    """Main function."""
    # Get the list of deleterious SNP IDs, offsets we want to keep, and the
    # polarities of the SNPs
    d_ids, offsets, polarities = anc_del_ids(ids, anc, dsnps)
    # Then, trim through the matrix
    s_ids, f_dos = del_dosage(offsets, polarities, dosages)
    # And print out the matrix
    print(' '.join(['IID'] + d_ids + ['TotDos']))
    for index, row in enumerate(f_dos):
        print(s_ids[index], ' '.join(row), sum([float(x) for x in row]))
    return


if len(sys.argv) != 5:
    print("""Calculate the dosage of deleterious variants across all lines in the
experiment. Filters to those that only have ancestral state. Takes four
arguments:
    1) SNP IDs in AlphaPeel dosages file
    2) Inferred ancestal states for each SNP (gzipped)
    3) Deleterious SNP IDs
    4) Combined AlphaPeel dosages output (gzipped)""")
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
