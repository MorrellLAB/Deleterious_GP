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
    Assign_Parental_Calls.py [parents.ped] [parents.map] [prog.ped] [prog.map]

Will write csv format data file for r/QTL, separated by family. Each data file
will have A, B, and H calls polarized by parent, with the maternal parent
represented by A, paternal by B, and H for heterozygotes."""
    print(msg)
    exit(1)


def parse_pedigrees(pedigree_file):
    """Reads a PED file, and returns a dictionary of the parental IDs,
    keyed on progeny ID."""
    pedigree = {}
    with open(pedigree_file, 'r') as f:
        for line in f:
            tmp = line.strip().split('\t')
            #   The second field is the line ID, the eighth is parent 1,
            #   and the ninth is parent 2
            famid = tmp[0]
            lineid = tmp[1]
            par1 = tmp[2]
            par2 = tmp[3]
            pedigree[lineid] = (par1, par2, famid)
    return pedigree


def parse_peds(ped, cycle):
    """Reads the PLINK ped files, and returns a dictionary of genotyping
    matrices, keyed on line ID."""
    ped_data = {}
    with open(ped, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            lineid = tmp[1]
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


def assign_genotypes(p1_geno, p2_geno, progeny, par_snpnames, prog_snpnames):
    """Assigns A for p1 allele, B for p2 for allele, and H for heterozygous
    sites. Returns a copy of `progeny` with the calls replaced."""
    #   Build a list of the indexes of the parental SNPs that are in the
    #   progeny.
    seg = [i
           for i in xrange(0, len(par_snpnames))
           if par_snpnames[i] in prog_snpnames]
    #   Iterate through the three genotype rows and assign the appropriate
    #   state to the progeny
    new_genotypes = []
    for i, j in enumerate(seg):
        p1_bases = (p1_geno[2*j], p1_geno[2*j + 1])
        p2_bases = (p2_geno[2*j], p2_geno[2*j + 1])
        prog_bases = (progeny[2*i], progeny[2*i + 1])
        print p1_bases, p2_bases, prog_bases
        if '0' in prog_bases:
            polarized = 'NA'
        #   Then, if either of the parents are heterozygous, then call it
        #   missing
        elif len(set(p1_bases)) != 1 or len(set(p2_bases)) != 1:
            polarized = 'NA'
        elif prog_bases == p1_bases:
            polarized = 'A'
        elif prog_bases == p2_bases:
            polarized = 'B'
        #   We will assume that if the progeny has two different alleles, then
        #   one came from each parent, and it is a true heterozygote
        elif len(set(prog_bases)) != 1:
            polarized = 'H'
        else:
            polarized = '-'
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
    handle = open(family + '.csv', 'w')
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


def main(parent_ped, parent_map, progeny_ped, progeny_map):
    """Main function."""
    #   Get the pedigree from the progeny PED file. These will just be the
    #   two parent IDs from the PED file.
    crosses = parse_pedigrees(progeny_ped)
    parental_genotypes = parse_peds(parent_ped, cycle=0)
    progeny_genotypes = parse_peds(progeny_ped, cycle=1)
    snpnames_par, snpchrom_par, snppos_par = parse_map(parent_map)
    snpnames_prog, snpchrom_prog, snppos_prog = parse_map(progeny_map)
    #   Assign A/B genotypes to each progeny
    family = {}
    for p in progeny_genotypes.iteritems():
        polarized_progeny = assign_genotypes(
            parental_genotypes[crosses[p[0]][0]],
            parental_genotypes[crosses[p[0]][1]],
            p[1],
            snpnames_par,
            snpnames_prog)
        if crosses[p[0]][2] not in family:
            family[crosses[p[0]][2]] = [polarized_progeny]
        else:
            family[crosses[p[0]][2]].append(polarized_progeny)
#   Then, iterate through each family and remove monomorphic markers
    for f, p in family.iteritems():
        print_csv(snpnames_prog, snpchrom_prog, snppos_prog, f, p)

if len(sys.argv) != 5:
    usage()
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
