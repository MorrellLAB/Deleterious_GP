#!/usr/bin/env python
"""Trim down the ancestral state file based on SNP ID. Takes two arguments:
    1) SNP ID list
    2) Ancestral table (gzipped)
"""

import sys
import gzip


def main(snpids, anc):
    """Main function."""
    s = []
    with open(snpids, 'r') as f:
        for line in f:
            s.append(line.strip())
    with gzip.open(anc, 'rb') as f:
        for line in f:
            if line.strip().split()[2] in s:
                print line.strip()
    return


if len(sys.argv) != 3:
    print """Trim down the ancestral state file based on SNP ID. Takes two arguments:
    1) SNP ID list
    2) Ancestral table (gzipped)"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
