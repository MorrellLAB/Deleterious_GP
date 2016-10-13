#   Prints out a list of SNP IDs that are predicted to be deleterious by all
#   three approaches.

#   Take arguments
args <- commandArgs(TRUE)
snp_table <- args[1]
snp_table <- read.table(snp_table, header=TRUE)

#   For the LRT, we apply a Bonferroni correction, with the number of codons
#   tested as the denominator. We use the "classic" 0.05 threshold.
n_codons <- sum(snp_table$Silent == "No")
sig <- 0.05 / n_codons
#   Additionally, we apply some heuristics to the masked constraint and nseqs
#   The constraint value is 95% specificity, as per LL's work.
max_const <- 0.91347390
minseq <- 10

#   For PROVEAN, the score should be less than or equal to -4.1528, for 95%
#   specificity, from LL's work.
provean_cutoff <- -4.1528

#   We don't have a similar set of filters for PolyPhen2, as it just prints a
#   Y/N prediction.

#   Get the deleterious SNPs by each approach. We also include premature stop
lrt <- snp_table[(snp_table$MaskedP.value <= sig & snp_table$SeqCount >= minseq & snp_table$MaskedConstraint <= max_const), "SNP_ID"]
provean <- snp_table[snp_table$PROVEAN <= provean_cutoff, "SNP_ID"]
pph2 <- snp_table[snp_table$PPH2 == "deleterious", "SNP_ID"]
nonsense <- snp_table[(snp_table$Silent == "No" & snp_table$AA1 != "*" & snp_table$AA2 == "*"), "SNP_ID"]

#   Remove NA, cast to a character
lrt <- as.character(lrt[!is.na(lrt)])
provean <- as.character(provean[!is.na(provean)])
pph2 <- as.character(pph2[!is.na(pph2)])
nonsense <- as.character(nonsense[!is.na(nonsense)])

#   Then, intersect them all!
a <- intersect(lrt, provean)
b <- intersect(a, pph2)

print(c(LRT=length(lrt), PROVEAN=length(provean), PPH2=length(pph2), Intersect=length(b), Nonsense=length(nonsense)))

#   append the nonsense SNPs to the vector of deleterious IDs
del <- c(b, nonsense)
write(del, file="Genomic_Prediction_Deleterious_IDs.txt", ncolumns=1)
