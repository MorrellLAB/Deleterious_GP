#!/usr/bin/env python
"""Takes a VCF, effects table, and a directory of multiple sequence alignments,
and generates a crude estimate of ancestral state. It will read the multuple
sequence alignment, and determine majority state among outgroup sequences. Takes
three arguments:
    1) VCF
    2) Effects table
    3) MSA directory
Requires Biopython.
"""

import sys
import os
from collections import Counter
try:
    from Bio import AlignIO
except ImportError:
    print 'This script requires Biopython.'
    exit(1)


def read_effects(eff):
    """Read the effects table and store the data in a dictionary. We will only
    store coding sites."""
    eff_dat = {}
    with open(eff, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split('\t')
                # Check if the SNP is coding or not. If it is, skip it
                if tmp[4] == '-':
                    continue
                else:
                    snpid = tmp[0]
                    txid = tmp[4]
                    unaligned_pos = tmp[10]
                    eff_dat[snpid] = (txid, unaligned_pos)
    return eff_dat


def get_outgroup_majority(msa_dir, txid, u_pos):
    """Get the outgroup bases at the supplied *unaligned* position."""
    # First, parse the alignment
    aln = os.path.join(
        os.path.expanduser(os.path.abspath(msa_dir)),
        txid + '_MSA.fasta')
    # If the file doesn't exist, we can't infer ancestral state
    if not os.path.isfile(aln):
        return ('N', 0)
    # Same for if the file is empty
    if os.stat(aln).st_size == 0:
        return ('N', 0)
    aln = AlignIO.read(aln, 'fasta')
    # Define a "safe name" that would be present in the alignment
    safe_name = txid.replace('.', '_')
    # Get the sequence ID of the query gene - we will not consider this position
    # for ancestral state inference.
    for index, rec in enumerate(aln):
        if rec.id == safe_name:
            ref_seq = index
            break
    # Then, identify the proper column to read.
    key = 1
    correct_col = -1
    for index, base in enumerate(str(aln[ref_seq].seq)):
        if key == u_pos:
            correct_col = index
            break
        elif base == '-':
            continue
        else:
            key += 1
    if correct_col == -1:
        return ('N', 0)
    # Grab the correct column
    col = aln[:, correct_col]
    # then, remove the query sequnce
    col = [c for i, c in enumerate(col) if i != ref_seq]
    # Turn it into a Counter object
    states = Counter(col)
    # Remove gaps
    del states['-']
    # If there are no bases left after removing gaps, then we return N
    if len(states) < 1:
        return ('N', 0)
    maj = states.most_common(1)[0][0]
    num_seqs = sum(states.values())
    return (maj, num_seqs)


def main(vcf, eff, msa_dir):
    """Main function."""
    print 'SNP_ID\tAnc_Base'
    eff_table = read_effects(eff)
    for snp in sorted(eff_table.keys()):
        anc, count = get_outgroup_majority(msa_dir, eff_table[snp][0], int(eff_table[snp][1]))
        print snp + '\t' + anc
    return


main(sys.argv[1], sys.argv[2], sys.argv[3])
