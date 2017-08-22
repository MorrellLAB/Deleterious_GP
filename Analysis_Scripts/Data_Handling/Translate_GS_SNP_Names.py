#!/usr/bin/env python
"""Translate the BOPA SNP names in a GenomeStudio export file into BOPAC names.
Takes two arguments:
    1) Translation table
    2) GenomeStudio report
"""

import sys
import pprint


def parse_translation(trans):
    """Parse the translation table and store it as a dictionary."""
    t = {}
    with open(trans, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split(',')
                # Store all names as keying to the BOPAC name. The BOPAC name
                # is in column 3
                for col in [0, 1, 3]:
                    if tmp[col]:
                        t[tmp[col]] = tmp[2]
    return t


def main(translation, gs):
    """Main function."""
    tab = parse_translation(translation)
    # Iterate through the GenomeStudio report and replace the SNP names
    with open(gs, 'r') as f:
        for line in f:
            if line.startswith('#'):
                print line.strip()
            else:
                tmp = line.strip().split('\t')
                snpid = tmp[1]
                # Get the SNP ID from the translation table
                new_id = tab.get(snpid, snpid)
                # Print the new rows
                print '\t'.join([tmp[0], new_id, tmp[2], tmp[3]])
    return


if len(sys.argv) != 3:
    print """Translate the BOPA SNP names in a GenomeStudio export file into BOPAC names.
Takes two arguments:
    1) Translation table
    2) GenomeStudio report"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
