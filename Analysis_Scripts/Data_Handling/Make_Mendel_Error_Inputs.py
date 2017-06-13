#!/usr/bin/env python
"""Take a large PLINK PED file and generate a bunch of family-specific files
that are properly formatted for PLINK Mendelian error checking. The parents
will have a family code that matches that of their progeny. Valid sex codes
for the parents will also be listed. Takes one argument:
    1) PED file with family info
"""

import sys


def parse_ped(pedfile):
    """Parse the PED file and store it in some complicated data structures. The
    first is a dictionary of {lineid: [PED_data]}. The other will be
    {family ID: [p1, p2, ind1, ind2, ...]}."""
    ped_data = {}
    families = {}
    with open(pedfile, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            fid = tmp[0]
            lid = tmp[1]
            # add the line to the ped_data
            ped_data[lid] = tmp
            # Then build up the pedigree
            if fid != '0':
                if fid not in families:
                    # Start off the list with the parents and the current
                    # progeny
                    families[fid] = [tmp[2], tmp[3], tmp[1]]
                else:
                    families[fid].append(lid)
    return (ped_data, families)


def generate_ped(fid, fam_members, ped_data):
    """Generate the family PED data to write for Mendel error checking. Replace
    the sex codes and the family IDs of the parents."""
    family_ped = []
    for m in fam_members:
        if m not in ped_data:
            continue
        p_dat = ped_data[m]
        if p_dat[0] != fid:
            p_dat[0] = fid
            p_dat[2] = '0'
            p_dat[3] = '0'
            # Append the maternal parent first
            if len(family_ped) == 0:
                p_dat[4] = '2'
            else:
                p_dat[4] = '1'
        family_ped.append(p_dat)
    return family_ped


def main(pedfile):
    """Main function."""
    ped_data, family_data = parse_ped(pedfile)
    for fam in sorted(family_data):
        famped = generate_ped(fam, family_data[fam], ped_data)
        # Send it to a file
        handle = open(fam + '_Mendel.ped', 'w')
        for ind in famped:
            handle.write(' '.join(ind) + '\n')
        handle.flush()
        handle.close()
    return


main(sys.argv[1])
