#!/usr/bin/env python3
"""Reads a single family PED file (with parental lines) and writes a csvs input
file for r/QTL. The genotypes are polarized by parent, with parent 1 getting the
A allele, and parent 2 getting the B allele. The parents are not included in the
output. The parents should be the first two lines in the PED file."""

import sys


def usage():
    """Usage message. Prints if no arguments."""
    msg = """
Usage:
    Assign_Parental_Calls.py [ped] [map]

Will write csv format data file for r/QTL, separated by family. Each data file
will have A, B, and H calls polarized by parent, with the maternal parent
represented by A, paternal by B, and H for heterozygotes."""
    print(msg)
    exit(1)


def parse_peds(ped):
    """Reads a PED file and returns a genotyping matrix with diploid calls."""
    ped_data = []
    with open(ped, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            calls = tmp[6:]
            nsnps = len(calls)/2
            genotypes = []
            for i in range(0, nsnps):
                dip = calls[2*i] + calls[1 + 2*i]
                genotypes.append(dip)
            ped_data.append(genotypes)
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


def assign_calls(ped_column):
    """Assigns A/B/H to each progeny given its genotype and its parents'
    genotypes. Heterozygous sites in the parents are treated as missing data."""
    parent_1 = list(set(ped_column[0]))
    parent_2 = list(set(ped_column[1]))
    #   If either of the parents are heterozygous, then the whole column gets
    #   missing data.
    if len(parent_1) != 1 or len(parent_2) != 1:
        return ['-'] * len(ped_column[2:])
    #   Then, get the parental alleles
    parent_1 = parent_1[0]
    parent_2 = parent_2[0]
    new_calls = []
    for prog in ped_column[2:]:
        #   If the parents are identical, skip the marker. But we have to use
        #   a placeholder value so that our genotypes do not get out of sync
        #   with the map
        if parent_1 == parent_2:
            new_calls.append('DROP')
        #   If the progeny is the same as p1, it gets 'A'
        elif prog == 2*parent_1:
            new_calls.append('A')
        #   p2 gets B
        elif prog == 2*parent_2:
            new_calls.append('B')
        #   hets get H
        elif prog == parent_1 + parent_2 or prog == parent_2 + parent_1:
            new_calls.append('H')
        #   If the progeny is NOT p1, but has p2, then we know it's not p1/p1,
        #   and gets C
        elif parent_1 not in prog and parent_2 in prog:
            new_calls.append('C')
        #   if not p2, then D
        elif parent_2 not in prog and parent_1 in prog:
            new_calls.append('D')
        #   Everything else is missing
        else:
            new_calls.append('-')
    return new_calls


def clean_matrix(genotypes, snpnames, snpchrom, snppos):
    """Removes markers that are monomorphic."""
    new_geno = []
    new_names = []
    new_chrom = []
    new_pos = []
    for g, n, c, p in zip(genotypes, snpnames, snpchrom, snppos):
        if len(set(g)) == 1:
            continue
        else:
            new_geno.append(g)
            new_names.append(n)
            new_chrom.append(c)
            new_pos.append(p)
    return (new_geno, new_names, new_chrom, new_pos)


def print_csv(genotypes, snpnames, snpchrom, snppos):
    """Prints the genotyping matrix in csv format for r/QTL."""
    #   First, print the marker names, with a 'Phe' in front for Phenotype. It
    #   will be missing, but it is part of the r/QTL format
    print ','.join(['Phe'] + snpnames)
    #   Then print the marker chromosomes. We need an empty field in front for
    #   the phenotype.
    print ','.join([''] + snpchrom)
    #   And then the marker positions
    print ','.join([''] + snppos)
    #   And then the genotyping data. We have to transpose it, since it is still
    #   column-oriented.
    for g in zip(*genotypes):
        print ','.join(['NA'] + list(g))
    return


def main(ped, snp_map):
    """Main function."""
    genotypes = parse_peds(ped)
    snpnames, snpchrom, snppos = parse_map(snp_map)
    #   Transpose the genotyping matrix
    genotypes_t = zip(*genotypes)
    #   Then assign A/B/H/- calls
    polarized = [assign_calls(col) for col in genotypes_t]
    cleaned_matrix, cleaned_names, cleaned_chrom, cleaned_pos = clean_matrix(
        polarized, snpnames, snpchrom, snppos)
    print_csv(cleaned_matrix, cleaned_names, cleaned_chrom, cleaned_pos)
    return

if len(sys.argv) != 3:
    usage()
else:
    main(sys.argv[1], sys.argv[2])
