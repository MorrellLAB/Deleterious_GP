#!/usr/bin/env python3
"""Takes PED files for the parental lines, and a CSV describing the pedigrees
of the individuals in the genomic prediction experiment, and writes r/QTL input
files (csvs format), with the genotypes polarized by parent. That is, the first
parent gets the A allele, and the second parent gets the B allele. Will also
remove markers that are monomorphic within familes. This script treats
heterozyogosity in the parents as missing data. This will have to be examined
more closely."""

import sys


def usage():
    """Usage message. Prints if no arguments."""
    msg = """
Usage:
    Assign_Parental_Calls.py [pedigrees] [map] [parents.ped] [progeny.ped]

Will write csv format data file for r/QTL, separated by family. Each data file
will have A, B, and H calls polarized by parent, with the maternal parent
represented by A, paternal by B, and H for heterozygotes."""
    print(msg)
    exit(1)


def parse_pedigrees(pedigree_file):
    """Reads the pedigrees file, and returns a dictionary of the parental IDs,
    keyed on line ID."""
    pedigree = {}
    with open(pedigree_file, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                #   Should be a comma separated values file
                tmp = line.strip().split(',')
                #   The second field is the line ID, the eighth is parent 1,
                #   and the ninth is parent 2
                lineid = tmp[1]
                par1 = tmp[7]
                par2 = tmp[8]
                pedigree[lineid] = (par1, par2)
    return pedigree


def parse_peds(ped, cycle):
    """Reads the PLINK ped files, and returns a dictionary of genotyping
    matrices, keyed on line ID."""
    ped_data = {}
    with open(ped, 'r') as f:
        for line in f:
            tmp = line.strip().split('\t')
            if cycle == 0:
                #   founder generation, just take the individual ID field out
                lineid = tmp[1]
            elif cycle == 1:
                #   Cycle 1 names are MS10S3 + family ID + individual number
                lineid = 'MS10S3' + tmp[0][1:] + '-0' + tmp[1]
            elif cycle == 2:
                lineid = 'MS11S2' + tmp[0] + '-' + tmp[1]
            elif cycle == 3:
                lineid = 'MS12_2' + tmp[0] + '-' + tmp[1]
            ped_data[lineid] = tmp[6:]
    return ped_data


def parse_map(p_map):
    """Parse a PLINK .map file. Returns the names and cM positions of the SNPs
    contained within it."""
    marker_names = []
    marker_chr = []
    marker_cm = []
    with open(p_map, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            #   Build up lists for outputting the map data as the header columns
            #   in an r/QTL csvs file. The columns in a MAP file are, from
            #   left to right:
            #       Chromosome
            #       Marker name
            #       Genetic map position
            #       Physical position
            marker_names.append(tmp[1])
            marker_chr.append(tmp[0])
            marker_cm.append(tmp[2])
    return(marker_names, marker_chr, marker_cm)


def assign_genotypes(p1_geno, p2_geno, progeny):
    """Assigns A for p1 allele, B for p2 for allele, and H for heterozygous
    sites. Returns a copy of `progeny` with the calls replaced."""
    #   Iterate through the three genotype rows and assign the appropriate
    #   state to the progeny
    new_genotypes = []
    for i in xrange(0, len(p1_geno), 2):
        p1_bases = (p1_geno[i], p1_geno[i+1])
        p2_bases = (p2_geno[i], p2_geno[i+1])
        prog_bases = (progeny[i], progeny[i+1])
        #   If the progeny have missing genotypes, then they should stay missing
        if '0' in prog_bases:
            polarized = 'NA'
        #   Then, if either of the parents are heterozygous, then call it
        #   missing
        elif len(set(p1_bases)) != 1 or len(set(p2_bases)) != 1:
            polarized = 'NA'
        #   If the parents are monomorphic, skip the marker
        elif p1_bases == p2_bases:
            #   Give it a dummy value
            polarized = 'monomorphic'
        elif prog_bases == p1_bases:
            polarized = 'A'
        elif prog_bases == p2_bases:
            polarized = 'B'
        #   We will assume that if the progeny has two different alleles, then
        #   one came from each parent, and it is a true heterozygote
        elif len(set(prog_bases)) != 1:
            polarized = 'H'
        new_genotypes.append(polarized)
    return new_genotypes


def print_csv(s_names, chrom, cm, family, progeny):
    """Prints the genotyping matrix in csv format for r/QTL."""
    to_drop = []
    #   We will iterate through the progeny genotypes, and build up a list of
    #   positions to drop from the genotyping matrix. Markers with 'NA' or
    #   'monomorphic' calls will be removed.
    for p in progeny:
        for index, marker in enumerate(p):
            if marker == 'monomorphic':
                to_drop.append(index)
    to_drop = sorted(list(set(to_drop)))
    #   Then, write the data out to a csv file
    famname = family[0] + '_x_' + family[1]
    handle = open(famname + '.csv', 'w')
    #   Get the marker names to keep
    final_snps = [s for i, s in enumerate(s_names) if i not in to_drop]
    #   And the chromosomes
    final_chr = [c for i, c in enumerate(chrom) if i not in to_drop]
    final_cm = [p for i, p in enumerate(cm) if i not in to_drop]
    genotypes = []
    #   Then iterate through the progeny matrix, and build the genotypes
    for i, p in enumerate(progeny):
        #   Start it off with NA, no phenotypic data.
        prog_gen = ['NA']
        for marker, c in enumerate(p):
            if marker in to_drop:
                continue
            elif c == 'NA':
                prog_gen.append('-')
            else:
                prog_gen.append(c)
        genotypes.append(prog_gen)
    #   Then, print everything out
    handle.write(
        ','.join(['Phe'] + final_snps) + '\n'
        )
    handle.write(
        ',' + ','.join(final_chr) + '\n'
        )
    handle.write(
        ',' + ','.join(final_cm) + '\n'
        )
    for g in genotypes:
        handle.write(','.join(g) + '\n')

    #   Close that handle
    handle.close()
    return


def main(pedigree, plink_map, parent, progeny):
    """Main function."""
    crosses = parse_pedigrees(pedigree)
    parental_genotypes = parse_peds(parent, cycle=2)
    progeny_genotypes = parse_peds(progeny, cycle=3)
    snpnames, snpchrom, snppos = parse_map(plink_map)
    families = {}
    for p in progeny_genotypes.items():
        prog_id = p[0]
        if prog_id not in crosses:
            #   We need to include this check because the genotyping data
            #   contains individuals that are not part of the genomic selection
            #   experiment. They will not show up in the pedigree, even though
            #   they are in the genotyping data.
            continue
        else:
            parents = crosses[prog_id]
        #   If, for some reason the parents aren't in the genotype matrix (it
        #   happens for M138, since it is missing every genotype call)
        if parents[0] not in parental_genotypes:
            continue
        if parents[1] not in parental_genotypes:
            continue
        polarized_progeny = assign_genotypes(
            parental_genotypes[parents[0]],
            parental_genotypes[parents[1]],
            p[1])
        if (parents[0], parents[1]) not in families:
            families[(parents[0], parents[1])] = [polarized_progeny]
        else:
            families[(parents[0], parents[1])].append(polarized_progeny)
    #   Then, iterate through each family and remove monomorphic markers, 
    for f, p in families.iteritems():
        print_csv(snpnames, snpchrom, snppos, f, p)

if len(sys.argv) != 5:
    usage()
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
