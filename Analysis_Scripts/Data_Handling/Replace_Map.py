#!/usr/bin/env python
"""Replace the GoldenPath physical positions with the RefSeqv1 positions."""

import sys


def main(goldenpath, refseq):
    """Main function. Replaces the physical positions column of the GoldenPath
    map with the physical positions from the RefSeq map. Keeps the marker order
    of the GoldenPath map."""

    #   Read the RefSeq positions, store in dict
    refseq_positions = {}
    with open(refseq, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split()
                #   Key is the SNPID, the value is (Chr, Physical Pos)
                refseq_positions[tmp[2]] = (tmp[0], tmp[1])

    #   Iterate through the GoldenPath map, print the new positions
    with open(goldenpath, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            if tmp[1] not in refseq_positions:
                #   Assign basepair position of 0, but keep the chromosome
                chrom = tmp[0]
                physical_pos = '0'
            else:
                chrom = refseq_positions[tmp[1]][0]
                physical_pos = refseq_positions[tmp[1]][1]
            print '\t'.join([chrom, tmp[1], tmp[2], physical_pos])
    return


if len(sys.argv) != 3:
    print """
Usage:
Replace_Map.py [GoldenPath Map] [RefSeq VCF]

will replace the physical positions in [GoldenPath Map] with those in
[RefSeq Map], while keeping the same order as in [GoldenPath Map]. If a marker
is not found in [RefSeq VCF], it will be assigned a basepair position of 0."""
    exit(1)

main(sys.argv[1], sys.argv[2])
