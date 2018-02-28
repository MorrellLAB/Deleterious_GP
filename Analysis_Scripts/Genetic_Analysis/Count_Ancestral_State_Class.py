#!/usr/bin/env python
"""Dumb script to count up the number of SNPs with inferred ancestral state
across functional categories. Takes no arguments: paths are hardcoded."""

import gzip
import sys

ANCESTRAL = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Ancestral_State/GP_Ancestral.txt.gz'
NONCODING = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Noncoding.txt'
SYNONYMOUS = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Synonymous.txt'
NONSYNONYMOUS = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Nonsynonymous.txt'
DELETERIOUS = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Deleterious.txt'

# Parse the lists and store SNP IDs
nonc = []
with open(NONCODING, 'r') as f:
    for line in f:
        nonc.append(line.strip())

syn = []
with open(SYNONYMOUS, 'r') as f:
    for line in f:
        syn.append(line.strip())

deleterious = []
with open(DELETERIOUS, 'r') as f:
    for line in f:
        deleterious.append(line.strip())


nonsyn = []
with open(NONSYNONYMOUS, 'r') as f:
    for line in f:
        if line.strip() in deleterious:
            continue
        else:
            nonsyn.append(line.strip())

# Keep track of the counts
nc_tot = 0
nc_anc = 0
s_tot = 0
s_anc = 0
ns_tot = 0
ns_anc = 0
del_tot = 0
del_anc = 0
with gzip.open(ANCESTRAL, 'rb') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split()
            snpid = tmp[2]
            sys.stderr.write(snpid + '\n')
            anc = tmp[3]
            if snpid in nonc:
                nc_tot += 1
            elif snpid in syn:
                s_tot += 1
            elif snpid in nonsyn:
                ns_tot += 1
            elif snpid in deleterious:
                del_tot += 1
            if anc == 'N':
                continue
            else:
                if snpid in nonc:
                    nc_anc += 1
                elif snpid in syn:
                    s_anc += 1
                elif snpid in nonsyn:
                    ns_anc += 1
                elif snpid in deleterious:
                    del_anc += 1

print 'Class Tot N_Anc Prop_Anc'
print 'All', nc_tot+s_tot+ns_tot+del_tot, nc_anc+s_anc+ns_anc+del_anc, float(nc_anc+s_anc+ns_anc+del_anc)/(nc_tot+s_tot+ns_tot+del_tot)
print 'Noncoding', nc_tot, nc_anc, float(nc_anc)/nc_tot
print 'Synonymous', s_tot, s_anc, float(s_anc)/s_tot
print 'Nonsynonymous', ns_tot, ns_anc, float(ns_anc)/ns_tot
print 'Deleterious', del_tot, del_anc, float(del_anc)/del_tot
