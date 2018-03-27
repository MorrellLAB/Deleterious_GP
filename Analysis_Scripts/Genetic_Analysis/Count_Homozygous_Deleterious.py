#!/usr/bin/env python
"""A really simple script to count the number of homozygous deleterious SNPs
per line in the genomic prediction experiment. A dosage of 2 is treated as
homozygous deleterious. Takes one argument:
    1) Deleterious dosages file
"""

import sys


def main(dosages):
    """Main function."""
    with open(dosages, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split()
                lid = tmp[0]
                nhom = 0
                for site in tmp[1:]:
                    if round(float(site), 1) == 2:
                        nhom += 1
                print lid, nhom
    return


if len(sys.argv) != 2:
    print """A really simple script to count the number of homozygous deleterious SNPs
per line in the genomic prediction experiment. A dosage of 2 is treated as
homozygous deleterious. Takes one argument:
    1) Deleterious dosages file"""
    exit(1)
else:
    main(sys.argv[1])
