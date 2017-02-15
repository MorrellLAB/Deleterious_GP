#!/usr/bin/env python
"""Translates a directory of coverage files from 'bedtools coverage' to a single
table for nice plotting with R. Mostly used for plotting coverage over well
characterized loci in 21 samples. Takes a diretory of coverage files as input.
"""

import sys
import os
import re
import pprint

#   Get all the files called '*Coverage.txt'
cov_files = os.listdir(sys.argv[1])
cov_files = [f for f in cov_files if f.endswith('Coverage.txt')]
samplenames = [s.split('_')[0] for s in cov_files]

#   Start a big data structure to hold the coverage data. We key on gene, since
#   that will our major key in the output. This is trickly because the original
#   data are keyed on sample.
#   {
#       Gene: {
#           Sample: [cov, cov, cov, ...],
#           Sample: [cov, cov, cov, ...],
#           ...
#           },
#       Gene: ...
#
coverage_data = {}
for f in cov_files:
    #   Get the sample name from the file
    samplename = f.split('_')[0]
    with open(os.path.join(sys.argv[1], f), 'r') as g:
        for line in g:
            #   Search for text in parenteses, remove the parens, make it upper
            genename = re.search(r'\(.+\)', line).group(0)[1:-1].upper()
            #   Start building the dictionary
            if genename not in coverage_data:
                coverage_data[genename] = {}
            #   Check the sample name, too
            if samplename not in coverage_data[genename]:
                coverage_data[genename][samplename] = []
            #   Append the stuff to the list
            coverage_data[genename][samplename].append(line.strip().split('\t')[-1])

print 'Gene\tPos\t' + '\t'.join(samplenames)
for gene in coverage_data:
    #   Then for each position in this gene...
    #       Since each gene length should be same across all samples (they were
    #       generated with the same BED file), this should be safe...
    genelength = len(coverage_data[gene]['6B01-2218'])
    for pos in range(1, genelength + 1):
        sitecov = []
        for sample in samplenames:
            sitecov.append(coverage_data[gene][sample][pos-1])
        print gene + '\t' + str(pos) + '\t' + '\t'.join(sitecov)
