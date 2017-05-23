#!/usr/bin/env python
"""Script to generate parent-of-origin BED files for each progeny line in the
genomic prediction experiment. Parses the output from GERMLINE to find segments
that are identical by descent. Takes two arguments:
    1) Pedigree CSV from Tyler
    2) GERMLINE match file
"""

import sys
import gzip
import pprint


# Define some constants for the lengths of chromosomes, and the first and
# last markers on each. The first and last positions are dependent on the
# genotyping platform. The ones listed here are for the 294 QC-passing BOPA
# SNPs for the genomic prediction experiment.
CHROMOSOMES = {
    '1H': {
        'Len': 558535432,
        'First': 2514495,
        'Last': 555426554
    },
    '2H': {
        'Len': 768075024,
        'First': 6420971,
        'Last': 764420842
    },
    '3H': {
        'Len': 699711114,
        'First': 3981172,
        'Last': 693425600
    },
    '4H': {
        'Len': 647060158,
        'First': 2789908,
        'Last': 646103540
    },
    '5H': {
        'Len': 670030160,
        'First': 1154365,
        'Last': 668031986
    },
    '6H': {
        'Len': 583380513,
        'First': 343390,
        'Last': 583345651
    },
    '7H': {
        'Len': 657224000,
        'First': 6314541,
        'Last': 654900595
    }
    }

def parse_pedigree(ped):
    """Parse the pedigree CSV, and save the founder IDs that contributed to each
    of the progeny lines."""
    pedigree = {}
    with open(ped, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split(',')
                famid = tmp[1].split('-')[0]
                if famid in pedigree:
                    continue
                else:
                    # We check the family ID. If it starts with MS10, then we
                    # are in cycle 1, and we save the parents. If it starts with
                    # MS11, it is cycle 2, and we need the grandparents. If it
                    # starts with MS12, then it is cycle 3 and we need the great
                    # grandparents.
                    if famid.startswith('MS10'):
                        pedigree[famid] = tmp[7:9]
                    elif famid.startswith('MS11'):
                        pedigree[famid] = tmp[9:13]
                    elif famid.startswith('MS12'):
                        pedigree[famid] = tmp[13:]
    return pedigree


def get_parent_of_origin(fam_id, parents, germline):
    """Parse the germline file and idetnify only matches between the parents
    and the lines in the specified family. This is not very efficient, since it
    requires us to parse the germline file (huge) for each family, but it'll
    have to do..."""
    # We will store the info in a dictionary with the following format:
    #   {
    #       line_id: {founder_1: [(chr, start, stop), (chr, start, stop), ... },
    #                {founder_2: [(chr, start, stop), (chr, start, stop), ... },
    #                ...
    #       line_id: ...
    #   }
    fam_ibd = {}
    with gzip.open(germline, 'rb') as f:
        for line in f:
            tmp = line.strip().split()
            # The fields are as follows:
            #   Family ID 1
            #   Individual ID 1
            #   Family ID 2
            #   Individual ID 2
            #   Chromosome
            #   Segment start (bp)
            #   Segment end (bp)
            #   Segment start (SNP)
            #   Segment end (SNP)
            #   Total SNPs in segment
            #   Genetic length of segment
            #   Units for genetic length (cM or MB)
            #   Mismatching SNPs in segment
            #   1 if Individual 1 is homozygous in match; 0 otherwise
            #   1 if Individual 2 is homozygous in match; 0 otherwise
            # We are only really interested in columns 0, 1, 2, 3, 4, 5, and 6
            fid1 = tmp[0]
            iid1 = tmp[1]
            fid2 = tmp[2]
            iid2 = tmp[3]
            chrom = tmp[4]
            start = tmp[5]
            end = tmp[6]
            # If start and end are identical, then we set them to be the
            # entire chromosome.
            if start == end:
                adj_start = 0
                adj_end = CHROMOSOMES[chrom]['Len']
            # Check the boundaries of the region. Extend them to the ends
            # of the chromosome if necessary
            if int(start) == CHROMOSOMES[chrom]['First']:
                adj_start = 0
            else:
                adj_start = int(start)
            if int(end) == CHROMOSOMES[chrom]['Last']:
                adj_end = CHROMOSOMES[chrom]['Len']
            else:
                adj_end = int(end)
            # Check the family IDs. One of them must be 0, and the other must
            # be the supplied family ID.
            fams = [fid1, fid2]
            if '0' in fams and fam_id in fams:
                # Next, check the IDs to make sure that founder is in the list
                # of parents of that family. The progeny all have individual
                # IDs that start with MS
                f_id = [i for i in [iid1, iid2] if not i.startswith('MS')][0]
                # If this founder is not in the list of parents, then we
                # continue
                if f_id not in parents:
                    continue
                else:
                    # Then, we want to get the progeny ID
                    p_id = [i for i in [iid1, iid2] if i != f_id][0]
                    # And we want to keep track of the IBD segments. First check
                    # if we are encountering a new progeny ID
                    if p_id not in fam_ibd:
                        fam_ibd[p_id] = {f_id: [(chrom, adj_start, adj_end)]}
                    else:
                        # Then, check if we are encountering a new founder ID
                        if f_id not in fam_ibd[p_id]:
                            fam_ibd[p_id][f_id] = [(chrom, adj_start, adj_end)]
                        else:
                            fam_ibd[p_id][f_id].append((chrom, adj_start, adj_end))
    return fam_ibd


def process_ibd(ibd):
    """Collapse any identical intervals and remove segments that overlap between
    parents."""
    
    pass


def main(ped, germline):
    """Main function."""
    pedigree = parse_pedigree(ped)
    t = get_parent_of_origin('MS10S3012', ['FEG183-52', 'FEG141-20'], germline)
    pprint.pprint(t)


main(sys.argv[1], sys.argv[2])
