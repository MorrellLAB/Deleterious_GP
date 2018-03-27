#!/usr/bin/env python
"""Read through the cycle-by-cycle genotypic counts and produce a file that
gives the cycle, functional class, and proportion of heterozygous genotypes
observed for each variant. Does not take arguments, but instead reads the files
in from hard-coded paths. We will print NAs for those that are monomorphic."""

import gzip

# Define paths here
C1_FRQ = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Genotype_Freqs/C1_Freq.frqx.gz'
C2_FRQ = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Genotype_Freqs/C2_Freq.frqx.gz'
C3_FRQ = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Genotype_Freqs/C3_Freq.frqx.gz'
NONC = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Noncoding.txt'
SYN = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Synonymous.txt'
NONSYN = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Nonsynonymous.txt'
DEL = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Deleterious.txt'

# Read the functional classes
# Parse the lists and store SNP IDs
nonc = []
with open(NONC, 'r') as f:
    for line in f:
        nonc.append(line.strip())

syn = []
with open(SYN, 'r') as f:
    for line in f:
        syn.append(line.strip())

deleterious = []
with open(DEL, 'r') as f:
    for line in f:
        deleterious.append(line.strip())


nonsyn = []
with open(NONSYN, 'r') as f:
    for line in f:
        if line.strip() in deleterious:
            continue
        else:
            nonsyn.append(line.strip())

c1_dat = {}
c2_dat = {}
c3_dat = {}
order = []
fcs = []
# Then iterate through the frequency files and calculate the proportion of
# heterozygous genotypes
with gzip.open(C1_FRQ, 'rb') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split()
            snpid = tmp[1]
            order.append(snpid)
            # Get the functional class
            if snpid in nonc:
                fc = 'Noncoding'
            elif snpid in syn:
                fc = 'Synonymous'
            elif snpid in nonsyn:
                fc = 'Nonsynonymous'
            elif snpid in deleterious:
                fc = 'Deleterious'
            else:
                fc = 'NA'
            fcs.append(fc)
            hom1 = int(tmp[4])
            het = int(tmp[5])
            hom2 = int(tmp[6])
            if (hom1 == 0 and het == 0) or (hom2 == 0 and het == 0):
                phet = 'NA'
            else:
                phet = str(float(het)/float(hom1 + hom2 + het))
            c1_dat[snpid] = phet

with gzip.open(C2_FRQ, 'rb') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split()
            snpid = tmp[1]
            # Get the functional class
            if snpid in nonc:
                fc = 'Noncoding'
            elif snpid in syn:
                fc = 'Synonymous'
            elif snpid in nonsyn:
                fc = 'Nonsynonymous'
            elif snpid in deleterious:
                fc = 'Deleterious'
            else:
                fc = 'NA'
            hom1 = int(tmp[4])
            het = int(tmp[5])
            hom2 = int(tmp[6])
            if (hom1 == 0 and het == 0) or (hom2 == 0 and het == 0):
                phet = 'NA'
            else:
                phet = str(float(het)/float(hom1 + hom2 + het))
            c2_dat[snpid] = phet

with gzip.open(C3_FRQ, 'rb') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split()
            snpid = tmp[1]
            # Get the functional class
            if snpid in nonc:
                fc = 'Noncoding'
            elif snpid in syn:
                fc = 'Synonymous'
            elif snpid in nonsyn:
                fc = 'Nonsynonymous'
            elif snpid in deleterious:
                fc = 'Deleterious'
            else:
                fc = 'NA'
            hom1 = int(tmp[4])
            het = int(tmp[5])
            hom2 = int(tmp[6])
            if (hom1 == 0 and het == 0) or (hom2 == 0 and het == 0):
                phet = 'NA'
            else:
                phet = str(float(het)/float(hom1 + hom2 + het))
            c3_dat[snpid] = phet

# Print a header
print 'SNPID\tFC\tCycle\tPHet'
for s, f in zip(order, fcs):
    print '\t'.join([s, f, 'C1', c1_dat.get(s, 'NA')])
    print '\t'.join([s, f, 'C2', c2_dat.get(s, 'NA')])
    print '\t'.join([s, f, 'C3', c3_dat.get(s, 'NA')])
