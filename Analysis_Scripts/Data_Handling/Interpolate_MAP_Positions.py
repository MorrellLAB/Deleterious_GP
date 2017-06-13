#!/usr/bin/env python
"""Interpolate the genetic and physical map positions of the SNPs with unknown
positions in the PLINK MAP file.
    1) MAP file
"""

import sys
import pprint


def parse_map(mapfile):
    """Return the MAP data as a list of lists. We do this so that we don't have
    to hit the disk over and over again."""
    map_data = []
    with open(mapfile, 'r') as f:
        for line in f:
            map_data.append(line.strip().split())
    return map_data


def missing_gen(mapdata):
    """Return a list of row numbers that have missing genetic positions in the
    MAP file."""
    miss = [
        idx
        for idx, row
        in enumerate(mapdata)
        if row[2] == '-9'
        ]
    return miss


def missing_phys(mapdata):
    """Return a list of row numbers that have missing physical positions in the
    MAP file."""
    miss = [
        idx
        for idx, row
        in enumerate(mapdata)
        if row[3] == '-9'
        ]
    return miss


def interpolate_pos(rows, mapdata, pos_type):
    """Interpolate the positions. Uses the information from the flanking
    markers to fill in values. If the genetic positions are identical for
    interpolating the physical position, we will use the midpoint of the
    flanking markers as the physical location."""
    interpolated = []
    for r in rows:
        upstream = mapdata[r-1]
        downstream = mapdata[r+1]
        if pos_type == 'G':
            g1 = float(upstream[2])
            g2 = float(downstream[2])
            p1 = float(upstream[3])
            pk = float(mapdata[r][3])
            p2 = float(downstream[3])
            d = (p2 - pk) / (p2 - p1)
            gu = g2 - d*(g2 - g1)
            if g2 == gu:
                interpolated.append(str(round((g2 + g1)/2, ndigits=2)))
            else:
                interpolated.append(str(round(gu, ndigits=2)))
        elif pos_type == 'P':
            # Extract upstream and downstream marker info
            g1 = float(upstream[2])
            gk = float(mapdata[r][2])
            g2 = float(downstream[2])
            p1 = float(upstream[3])
            p2 = float(downstream[3])
            # Calculate the interpolated position, as a proportion of the
            # physical distance between the flanking markers
            d = (g2 - gk) / (g2 - g1)
            pu = p2 - d*(p2 - p1)
            # If pu == p2, then return the midpoint of the physical pos
            if pu == p2:
                interpolated.append(str(int(round((p2 + p1)/2))))
            else:
                interpolated.append(str(int(round(pu))))
    return interpolated


def main(mapfile):
    """Main function."""
    map_data = parse_map(mapfile)
    no_phys = missing_phys(map_data)
    i_phys = interpolate_pos(no_phys, map_data, 'P')
    # Put the interpolated values in, and re-sort
    for r, i in zip(no_phys, i_phys):
        map_data[r][3] = i
    resorted = sorted(map_data, key=lambda x: (x[0], int(x[3])))
    # Then, find the rows with missing genetic positions
    no_gen = missing_gen(resorted)
    # And interpolate them
    i_gen = interpolate_pos(no_gen, resorted, 'G')
    for r, i in zip(no_gen, i_gen):
        resorted[r][2] = i
    # Then re-sort on genetic positions
    reresorted = sorted(resorted, key=lambda x: (x[0], float(x[2])))
    # And print it out
    for row in reresorted:
        print '\t'.join(row)


main(sys.argv[1])
