#!/usr/bin/env python
"""Convert from Alchemy output format to PLINK format for easy data handling.
Takes two arguments:
    1) The Alchemy report
    2) the physical/genetic map information from HarvEST (SNP_BAC.txt).
Note: the SNPs listed in the Alchemy report and those from HarvEST should have
the same naming scheme (probably BOPAC). This script will not check that."""

import sys


class PedFile(object):
    """Data for the PED file. Will store the six required columns and will
    build genotype columns from the input matrix."""
    #   These are the six fixed columns, which must be populated for each
    #   individual
    ped_data = {}
    #   Some constants use to define missing values
    pheno_miss = '-9'
    missing = '0'

    def __init__(self):
        pass

    def initialize_ped(self, ind_id):
        """Create the six required fields from the sample name, and start a
        dictionary to hold genetic data. We use a dictionary so that we can
        ensure agreement between the PED and the MAP files."""
        #   Watch out for duplicate entries
        if ind_id in self.ped_data or ind_id == 'Blank':
            return
        else:
            #   Note: the way we slice up the sample name to get family ID and
            #   individual ID will vary from cycle to cycle. Naming conventions
            #   are slightly different.
            #       Cycle 1: G10WXXX-YY
            #           G: ?
            #           10: 2010
            #           W: ?
            #           XXX: Family ID
            #           YY: Line ID
            #       Cycle 2: MS11S2XXX-YYY
            #           M: Minnesota
            #           S: Spring
            #           11: 2011
            #           S2: 2nd gen. selfing (F3)
            #           XXX: Family ID
            #           YYY: Line ID
            #       Cycle 3: MS12_2XXX-YYY
            #           M: Minnesota
            #           S: Spring
            #           12: 2012
            #           2: S2 (F3)
            #           XXX: Family ID
            #           YYY: Line ID
            self.ped_data[ind_id] = {
                'familyid': ind_id[6:9],
                'individualid': ind_id[-3:],
                'paternalid': self.missing,
                'maternalid': self.missing,
                'sex': self.missing,
                'phenotype': self.pheno_miss,
                'genotypes': {}}
            return

    def assign_genotypes(self, snpid, sample, call, prob, cutoff=0.8):
        """Assign calls to the SNPs for each sample."""
        #   First, check that the sample has a row. If not, continue
        if sample not in self.ped_data:
            return
        else:
            #   Then, start building up the genotypes dictionary
            #   Check that the probability of a call is above the cutoff. If it
            #   is less than the cutoff, we assign a missing genotype. We use
            #   the two calls separated in a list since PED files store diploid
            #   data.
            if float(prob) < cutoff:
                gen = [self.missing, self.missing]
            else:
                gen = [call[0], call[1]]
            self.ped_data[sample]['genotypes'][snpid] = gen
            return


class MapFile(object):
    """A class for the MAP file used by PLINK. Stores SNP ID, the genetic map
    position, the physical map position, and the chromosome."""

    def __init__(self, harvest):
        """Uses the data from HarvEST to store the necessary information for a
        PLINK MAP file:
            Chromosome
            SNPID
            Genetic Map Pos
            Physical Map Pos
        """
        self.genetic_data = {}
        with open(harvest, 'r') as f:
            for index, line in enumerate(f):
                if index == 0:
                    continue
                else:
                    tmp = line.strip().split('\t')
                    #   The columns we want are
                    #   1: SNP ID
                    #   14: chromosome (2011 map)
                    #   15: cM (2011 map)
                    #   11: bp pos (golden path)
                    if tmp[0] in self.genetic_data:
                        continue
                    else:
                        snpid = tmp[0]
                        chrom = tmp[13]
                        cm = tmp[14]
                        gp = tmp[10]
                        #   Then go through and replace empty values with
                        #   missing values.
                        if chrom == '':
                            chrom = '0'
                        if cm == '':
                            cm = '0'
                        if gp == '':
                            gp = '0'
                        self.genetic_data[snpid] = (
                            chrom,
                            cm,
                            gp)
        return

    def sort_map(self):
        """Sorts the genetic map by chromosome, then by genetic map position.
        Returns a list of SNPs in that order."""
        sorted_snps = sorted(
            self.genetic_data.items(),
            key=lambda x: (x[1][0], float(x[1][1]))
            )
        return(sorted_snps)


def main(alchemy, pg_map):
    """The main function. Builds output names for the .ped and the .map files,
    and fills them with the appropriate data."""
    #   Start new classes for the Ped and Map files
    p = PedFile()
    m = MapFile(pg_map)
    with open(alchemy, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split()
                p.initialize_ped(tmp[1])
                p.assign_genotypes(
                    tmp[0],
                    tmp[1],
                    tmp[3],
                    tmp[4])
    s = m.sort_map()
    #   Then, write the PED and MAP files
    ped_out = open(alchemy.replace('.txt', '.ped'), 'w')
    map_out = open(alchemy.replace('.txt', '.map'), 'w')
    for sample in sorted(p.ped_data):
        fid = p.ped_data[sample]['familyid']
        iid = p.ped_data[sample]['individualid']
        pid = p.ped_data[sample]['paternalid']
        mid = p.ped_data[sample]['maternalid']
        sex = p.ped_data[sample]['sex']
        phe = p.ped_data[sample]['phenotype']
        gen = []
        for snp, gmap in s:
            gen += p.ped_data[sample]['genotypes'][snp]
        #   Make the genotypes a string so that we can use it with join()
        gen = '\t'.join(gen)
        #   Write the PED line
        ped_out.write(
            '\t'.join(
                [
                    fid,
                    iid,
                    pid,
                    mid,
                    sex,
                    phe,
                    gen]
                ) + '\n'
            )
    #   And then write the MAP file
    for snp, gmap in s:
        map_out.write(
            '\t'.join(
                [
                    gmap[0],
                    snp,
                    gmap[1],
                    gmap[2]]
                ) + '\n'
            )
    #   Flush and close the handles
    ped_out.flush()
    ped_out.close()
    map_out.flush()
    map_out.close()
    return


main(sys.argv[1], sys.argv[2])
