#!/usr/bin/env python
"""Make the AlphaPeel map file. It is a somewhat strange format. See the manual
and the accompanying MarkDown document to learn more. Takes two arguments:
    1) Combined BOPA+ExCap VCF with BOPA names
    2) Chromosome
"""

import sys
import gzip


def read_positions(v, c):
    """Read the positions out of the VCF. Keep only markers on the chromosome
    of interest. We will keep a tuple of marker information:
    (phys pos, interpolate)"""
    positions = []
    with gzip.open(v, 'rb') as f:
        for line in f:
            if line.startswith(c):
                tmp = line.strip().split()
                pos = int(tmp[1])
                snp_id = tmp[2]
                # This is somewhat kludge-y, but if a SNP has 11_ or 12_ in the
                # name, then we know its position, and it will be an "anchor"
                # in the map
                if '11_' in snp_id or '12_' in snp_id:
                    interpolate = False
                else:
                    interpolate = True
                positions.append((pos, interpolate))
            else:
                continue
    return positions


def find_flanks(pos):
    """Identify the flanking markers with known position from the tuple of VCF
    position information."""
    left = []
    right = []
    known = 0
    max_known = sum([not p[1] for p in pos])
    for p in pos:
        if p[1]:
            l = known
            r = known+1
        else:
            l = known+1
            r = known+1
            known += 1
        if l < 1:
            l = 1
        if r < 1:
            r = 1
        if l > max_known:
            l = max_known
        if r > max_known:
            r = max_known
        left.append(l)
        right.append(r)
    return (zip(left, right))


def interpolate(pos, flank):
    """Give a relative position between markers of known position. If the
    two flanking marker indices are equal, then return 0."""
    # First, extract the merkers of known position to look up later
    known_pos = []
    for p in pos:
        if not p[1]:
            known_pos.append(p[0])
    # Then, iterate through the positions and interpolate the ones that need it
    interpolated = []
    for p, f in zip(pos, flank):
        if f[0] == f[1]:
            interpolated.append('0')
        else:
            # Get the known positions out of the known list. Because the marker
            # indices are 1-based, we have to subtract 1 to get the correct
            # position.
            left = known_pos[f[0] - 1]
            right = known_pos[f[1] - 1]
            i = (p[0] - left) / float((right - left))
            interpolated.append('{:.10f}'.format(i))
    return interpolated


def main(vcf, chromosome):
    """Main function."""
    # Get a list of the variant positions
    p = read_positions(vcf, chromosome)
    flanks = find_flanks(p)
    relpos = interpolate(p, flanks)
    # Then print it out
    for f, r in zip(flanks, relpos):
        print f[0], f[1], r


if len(sys.argv) != 3:
    print """Make the AlphaPeel map file. It is a somewhat strange format. See the manual
and the accompanying MarkDown document to learn more. Takes two arguments:
    1) Combined BOPA+ExCap VCF with BOPA names
    2) Chromosome"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
