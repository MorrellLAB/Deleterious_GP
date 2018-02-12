# Use a chi-squared test to calculate the probability of segregation distortion
# in the genomic prediction panel.

# Read in the data files
c1_obs <- read.table(gzfile("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Genotype_Freqs/C1_Freq.frqx.gz"), skip=1, header=FALSE)
c1_exp <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Genotype_Freqs/Cycle_1_ExpGeno.txt", header=TRUE)
c2_obs <- read.table(gzfile("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Genotype_Freqs/C2_Freq.frqx.gz"), skip=1, header=FALSE)
c2_exp <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Genotype_Freqs/Cycle_2_ExpGeno.txt", header=TRUE)
c3_obs <- read.table(gzfile("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Genotype_Freqs/C3_Freq.frqx.gz"), skip=1, header=FALSE)
c3_exp <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Genotype_Freqs/Cycle_3_ExpGeno.txt", header=TRUE)
bim <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Imputation/AlphaPeel/AlphaPeel_Imputed_Merged.bim", header=FALSE)

# Chromosomes
chroms <- c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H")

# Set column names for the observed counts
names(c1_obs) <- c("CHR", "SNP", "A1", "A2", "AA", "AB", "BB", "HAPAA", "HAPBB", "MISS")
names(c2_obs) <- c("CHR", "SNP", "A1", "A2", "AA", "AB", "BB", "HAPAA", "HAPBB", "MISS")
names(c3_obs) <- c("CHR", "SNP", "A1", "A2", "AA", "AB", "BB", "HAPAA", "HAPBB", "MISS")

# Trim down the observed data
c1_obs <- c1_obs[,c("CHR", "AA", "AB", "BB")]
c2_obs <- c2_obs[,c("CHR", "AA", "AB", "BB")]
c3_obs <- c3_obs[,c("CHR", "AA", "AB", "BB")]

# Define a function to run a chi-sq test and return a P-value
seg_dist <- function(x) {
    obs <- as.numeric(x[c("AA", "AB", "BB")])
    expect <- as.numeric(x[c("ExpAA", "ExpAB", "ExpBB")])
    # pval <- chisq.test(x=obs, p=expect, rescale.p=TRUE)$p.value
    pval <- fisher.test(obs, expect)$p.value
    return(pval)
}

# Put the observed and expected together
c1 <- cbind(c1_obs, c1_exp)
c1_pval <- apply(c1, 1, seg_dist)
# Multiply the p-values by the number of tests
c1_pval <- c1_pval * length(c1_pval)
# Take the -log of the p-value
logc1_pval <- -log(c1_pval)
# Clip the values
logc1_pval[logc1_pval > 10] <- 10
logc1_pval[logc1_pval < 0] <- 0
# Save the values for lookup
tosave <- cbind(c1_obs, c1_exp, c1_pval/length(c1_pval))
names(tosave) <- c("CHR", "AA", "AB", "BB", "SNPName", "ExpAA", "ExpAB", "ExpBB", "PVal")
write.csv(tosave, file="Cycle_1_SegDist.csv", quote=F, row.names=F)

c2 <- cbind(c2_obs, c2_exp)
c2_pval <- apply(c2, 1, seg_dist)
# Multiply the p-values by the number of tests
c2_pval <- c2_pval * length(c2_pval)
# Take the -log of the p-value
logc2_pval <- -log(c2_pval)
# Clip the values
logc2_pval[logc2_pval > 10] <- 10
logc2_pval[logc2_pval < 0] <- 0
tosave <- cbind(c2_obs, c2_exp, c2_pval/length(c2_pval))
names(tosave) <- c("CHR", "AA", "AB", "BB", "SNPName", "ExpAA", "ExpAB", "ExpBB", "PVal")
write.csv(tosave, file="Cycle_2_SegDist.csv", quote=F, row.names=F)

c3 <- cbind(c3_obs, c3_exp)
c3_pval <- apply(c3, 1, seg_dist)
# Multiply the p-values by the number of tests
c3_pval <- c3_pval * length(c3_pval)
# Take the -log of the p-value
logc3_pval <- -log(c3_pval)
# Clip the values
logc3_pval[logc3_pval > 10] <- 10
logc3_pval[logc3_pval < 0] <- 0
tosave <- cbind(c3_obs, c3_exp, c3_pval/length(c3_pval))
names(tosave) <- c("CHR", "AA", "AB", "BB", "SNPName", "ExpAA", "ExpAB", "ExpBB", "PVal")
write.csv(tosave, file="Cycle_3_SegDist.csv", quote=F, row.names=F)

# Plot it by chromosome
pdf(file="C1_Seg_Dist.pdf", height=10.5, width=8)
par(mfrow=c(7, 1), mar=c(3, 4, 1, 2))
for(i in chroms) {
    sd_chr <- logc1_pval[bim$V1 == i]
    bppos <- bim$V4[bim$V1 == i]/1000000
    plot(sd_chr ~ bppos, main=paste("Cycle 1 ", i, sep=""), type="n", ylab="Adjusted -log(P)", xlab="Physical Position (Mb)", xlim=c(0, 750), xaxt="n")
    lines(sd_chr~bppos, lwd=0.25)
    lines(lowess(sd_chr ~ bppos, f=0.01), col="red", lwd=2)
    abline(h=-log(0.01), col="blue", lwd=2)
    axis(side=1, at=seq(0, 750, by=50), labels=seq(0, 750, by=50))
}
dev.off()

pdf(file="C2_Seg_Dist.pdf", height=10.5, width=8)
par(mfrow=c(7, 1), mar=c(3, 4, 1, 2))
for(i in chroms) {
    sd_chr <- logc2_pval[bim$V1 == i]
    bppos <- bim$V4[bim$V1 == i]/1000000
    plot(sd_chr ~ bppos, main=paste("Cycle 2 ", i, sep=""), type="n", ylab="Adjusted -log(P)", xlab="Physical Position (Mb)", xlim=c(0, 750), xaxt="n")
    lines(sd_chr~bppos, lwd=0.25)
    lines(lowess(sd_chr ~ bppos, f=0.01), col="red", lwd=2)
    abline(h=-log(0.01), col="blue", lwd=2)
    axis(side=1, at=seq(0, 750, by=50), labels=seq(0, 750, by=50))
}
dev.off()

pdf(file="C3_Seg_Dist.pdf", height=10.5, width=8)
par(mfrow=c(7, 1), mar=c(3, 4, 1, 2))
for(i in chroms) {
    sd_chr <- logc3_pval[bim$V1 == i]
    bppos <- bim$V4[bim$V1 == i]/1000000
    plot(sd_chr ~ bppos, main=paste("Cycle 3 ", i, sep=""), type="n", ylab="Adjusted -log(P)", xlab="Physical Position (Mb)", xlim=c(0, 750), xaxt="n")
    lines(sd_chr~bppos, lwd=0.25)
    lines(lowess(sd_chr ~ bppos, f=0.01), col="red", lwd=2)
    abline(h=-log(0.01), col="blue", lwd=2)
    axis(side=1, at=seq(0, 750, by=50), labels=seq(0, 750, by=50))
}
dev.off()
