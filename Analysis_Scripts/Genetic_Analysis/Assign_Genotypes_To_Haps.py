#!/usr/bin/env python
"""Print a PLINK PED file with the imputed genotypes from AlphaPeel. We will
use the most likely haplotype configurating for the PED file. Takes four
arguments:
    1) VCF with ref/alt for the exome capture SNPs (Gzipped)
    2) Pedigree that does not include inbreeding (Gzipped)
    3) AlphaPeel haps output file (Gzipped)
    4) Chromosome
"""

import sys
import gzip


def parse_vcf(v, c):
    """Parse the VCF and store SNP ID, physical position, ref, and alt alleles.
    Restrict it to only the specified chromosome."""
    vcf_dat = []
    with gzip.open(v, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            elif not line.startswith(c):
                continue
            else:
                tmp = line.strip().split()
                pos = tmp[1]
                snpid = tmp[2]
                ref = tmp[3]
                alt = tmp[4]
                vcf_dat.append((pos, snpid, ref, alt))
    return vcf_dat


def parse_ped(p):
    """Parse the pedigree and return a dictionary that gives the family ID,
    maternal ID, and paternal ID keyed on individual ID."""
    ped = {}
    with gzip.open(p, 'rt') as f:
        for line in f:
            tmp = line.strip().split()
            if tmp[0].startswith('MS'):
                famid = tmp[0].split('-')[0]
            else:
                famid = '0'
            ped[tmp[0]] = (famid, tmp[1], tmp[2])
    return ped


def make_ped_string(h, v, p):
    """Iterate through the haps file, four at a time and taket he most likely
    haplotype to print. Use the ref and alt alleles from the VCF data, and the
    family info from the pedigree to make a complete PED file."""
    with gzip.open(h, 'rt') as f:
        # Iterate through the file four lines at a time. The first element is
        # the probability of homozygous reference, then the two hetereozygous
        # phases, then homozygous alt.
        for rr, ra, ar, aa in zip(f, f, f, f):
            # Skip internediate "nuisance" generations
            if rr.startswith('Fam'):
                continue
            ped_string = []
            for index, hap_probs in enumerate(zip(rr.split(), ra.split(), ar.split(), aa.split())):
                if index == 0:
                    iid = hap_probs[0]
                    if iid not in p:
                        famid = '0'
                        mid = '0'
                        pid = '0'
                    else:
                        famid, mid, pid = p[iid]
                    sex = '0'
                    phen = '0'
                    ped_string += [famid, iid, pid, mid, sex, phen]
                else:
                    # Make the genotypes from the VCF info
                    ref = v[index-1][2]
                    alt = v[index-1][3]
                    genotypes = [
                        ref + ' '+ ref,
                        ref + ' ' + alt,
                        alt + ' ' + alt]
                    probs = [
                        float(hap_probs[0]),
                        float(hap_probs[1]) + float(hap_probs[2]),
                        float(hap_probs[3])]
                    # Append the correct haplotype to the ped
                    highest = max(probs)
                    ped_string.append(genotypes[probs.index(highest)])
            print(' '.join(ped_string))
    return


def make_map(chrom, m):
    """Make a PLINK map file for Mendel error checking. Write to stderr."""
    for line in m:
        snpid = line[1]
        gpos = '-9'
        ppos = line[0]
        sys.stderr.write(' '.join([chrom, snpid, gpos, ppos]) + '\n')
    return


def main(vcf, ped, haps, chrom):
    """Main function."""
    # First, parse the VCF
    v = parse_vcf(vcf, chrom)
    # Then parse the pedigree
    p = parse_ped(ped)
    # Then iterate through the haplotypes and print the PED string. We nedd to
    # know the ref and alt for each
    make_ped_string(haps, v, p)
    # Next, make the map file. Write this to stderr.
    make_map(chrom, v)

if len(sys.argv) != 5:
    print("""Print a PLINK PED file with the imputed genotypes from AlphaPeel. We will
use the most likely haplotype configurating for the PED file. Takes four
arguments:
    1) VCF with ref/alt for the exome capture SNPs (Gzipped)
    2) Pedigree that does not include inbreeding (Gzipped)
    3) AlphaPeel haps output file (Gzipped)
    4) Chromosome""")
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
