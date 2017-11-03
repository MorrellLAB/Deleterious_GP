#!/usr/bin/env python
"""Add the genetic positions from the Munoz et al 2011 consensus BOPA SNP map
to a PLINK MAP. Assumes that the genetic positions are stored in a CSV file
with the following columns: SNP_Name, Chromosome, Position. Takes two
arguments:
    1) MAP file to edit
    2) Geneic map positions
"""

import sys

def parse_map(gmap):
    """Parse the consensus genetic map and store it in a dictionary."""
    gpos = {}
    with open(gmap, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split(',')
                gpos[tmp[0]] = (tmp[1], tmp[2])
    return gpos


def main(plinkmap, gmap):
    """Main function."""
    # get the genetic map positions
    genetic_pos = parse_map(gmap)
    # Iterate through the PLINK MAP and priunt out new chromosomes and
    # genetic positions
    with open(plinkmap, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            chrom = tmp[0]
            snpid = tmp[1]
            if chrom == 'chrUn':
                phys = '-9'
            else:
                phys = tmp[3]
            if snpid in genetic_pos:
                toprint = [genetic_pos[snpid][0], snpid, genetic_pos[snpid][1], phys]
            else:
                toprint = [chrom.replace('chr', ''), snpid, '-9', phys]
            print '\t'.join(toprint)
    return


main(sys.argv[1], sys.argv[2])
