#!/usr/bin/env python
"""Order a PLINK map file by physical position, and interpolate the genetic
positions that are out of order. BEAGLE requires that the genetic positions
be montonically increasing with physical position. Takes one argument:
    1) MAP file
"""

import sys
import math


def parse_map(plinkmap):
    """Parse the PLINK map and store the markers in a dictioanry of lists. The
    keys are the chromosomes, and the lists hold the marker orders."""
    plink = {}
    with open(plinkmap, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            chrom = tmp[0]
            if chrom not in plink:
                plink[chrom] = []
            plink[chrom].append(tmp)
    # Then sort on physical position
    for c in plink:
        plink[c] = sorted(plink[c], key=lambda x: int(x[3]))
    return plink


def find_bad_order(plinkmap):
    """Iterate through the sorted map and identify markers that need to be
    interpolated (inverted) or jittered (equal pos)."""
    modif = {}
    for chrom in plinkmap:
        # For each pair of adjacent markers, identify the type of modification
        # that needs to happen
        if chrom not in modif:
            modif[chrom] = []
        # Iterate from 1 to the number of markers on the chromosome.
        for i in range(1, len(plinkmap[chrom])):
            g2 = float(plinkmap[chrom][i][2])
            g1 = float(plinkmap[chrom][i-1][2])
            if g2 - g1 < 0:
                modif[chrom].append('I')
            elif g2 - g1 == 0:
                modif[chrom].append('J')
            else:
                modif[chrom].append(None)
    return modif


def mod_map(mods, plinkmap):
    """Apply the modificiations that are necessary to the map."""
    modmap = {}
    for chrom in plinkmap:
        if chrom not in modmap:
            modmap[chrom] = []
        markers = plinkmap[chrom]
        modif = mods[chrom]
        for i, m in enumerate(modif):
            if m == 'I':
                p2 = float(markers[i+1][3])
                p1 = float(markers[i-1][3])
                pk = float(markers[i][3])
                g2 = float(markers[i+1][2])
                g1 = float(markers[i-1][2])
                d = (p2 - pk) / (p2 - p1)
                gu = g2 - d*(g2 - g1)
                if g2 == gu:
                    gi = str(round((g2 + g1)/2, ndigits=2))
                else:
                    gi = str(round(gu, ndigits=2))
                modmar = [markers[i][0], markers[i][1], gi, markers[i][3]]
            elif m == 'J':
                jgpos = marker[i][2] + '1'
                modmar = [markers[i][0], markers[i][1], jgpos, markers[i][3]]
            else:
                modmar = markers[i]
            modmap[chrom].append(modmar)
    return modmap

def main(plinkmap):
    """Main function."""
    pmap = parse_map(plinkmap)
    # Then, idetnify markers that are out of order
    modifications = find_bad_order(pmap)
    fixedmap = mod_map(modifications, pmap)
    for c in sorted(fixedmap):
        for m in fixedmap[c]:
            print '\t'.join(m)
    return


if len(sys.argv) != 2:
    print """Order a PLINK map file by physical position, and interpolate the genetic
positions that are out of order. BEAGLE requires that the genetic positions
be montonically increasing with physical position. Takes one argument:
    1) MAP file"""
    exit(1)
else:
    main(sys.argv[1])
