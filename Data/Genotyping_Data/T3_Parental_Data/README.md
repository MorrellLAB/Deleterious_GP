# Parental Genotype Data From T3
The parental BOPA data was downloaded as a series of genotyping experiments
because it was the only way to ensure that the genotypes would be in the same
orientiation (strand) between the parents and progeny. Genotyping experiments
were downloaded with no filtering applied: the maximum per-marker missing data
should be set to 100%, and the minimum minor allele frequency should be set
to 0%.

To download genotyping experiments, go to the T3 home page and click on
"Genotype Experiments" on the left. Then choose the batch of genotypes to
download. Set the maximum missing data to 100%, and the minimum MAF to 0%, and
click the "refresh" button to apply the filters. Then click "download allele
data" then "download Zip file of results" to get the data. It will serve a
zip file with a hapmap file (genotype.hmp.txt) inside with the genotype data.

Note that MN in 2009 has only the 'BOPA1' chip associated with it - this is an
error in T3, it looks like. The 'BOPA1' data has both BOPA1 and BOPA2 SNPs.

Also note: the 'ND' data was the 6-row data (N6) in T3.

## Chips With the Parents
The genotpying experiments that contain the parents of the genomic prediction
population are as follows:

- 2009 MN
- 2009 BA
- 2009 ND
- 2008 MN
    - FEG183-52
- 2008 BA
    - 6B06-1132
- 2008 ND
    - ND20448
    - ND24906
    - ND25160
    - ND25986
    - ND26036
    - ND26104
- 2007 MN
    - FEG153-58
    - FEG154-47
    - FEG175-57
- 2007 BA
    - 6B04-0290
    - 6B05-0922
- 2007 ND
    - ND25652
    - ND25728
- 2006 MN
    - FEG141-20
- 2006 BA
    - 6B01-2218
    - 6B03-4304
    - 6B03-4478
- 2006 ND

M122 and M138 are not present in this genotyping data. M138 does not have any
genotyping data in T3 at all. M122 has only been genotyped on the 9K and BOPA1
chips.
