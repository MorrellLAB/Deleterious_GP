#!/usr/bin/env python
"""Takes a VCF, effects table, and a directory of multiple sequence alignments
to output an input file for use with Peter Keightley's program. The multiple
sequence alignments must have two outgroups, and one ingroup. Takes three
arguments:
    1) VCF
    2) Effects table
    3) MSA directory
Requires Biopython.
"""

import sys
import os
import pprint
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


def get_outgroup_bases(msa_dir, txid, u_pos):
    """Get the outgroup bases at the supplied *unaligned* position."""
    # First, parse the alignment
    aln = os.path.join(
        os.path.expanduser(os.path.abspath(msa_dir)),
        txid + '_MSA.fasta')
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
    pass


def main(vcf, eff, msa_dir):
    """Main function."""
    eff_table = read_effects(eff)
    pprint.pprint(eff_table)
    return


main(sys.argv[1], sys.argv[2], sys.argv[3])
