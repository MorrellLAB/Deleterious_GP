# Plot empirical cumulative distribution functions for the P-values from
# PLINK quantitative association analyses. The SNPs have been separated by
# functional class.

yld.nc <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/PLINK_qAssoc/Yield/Noncoding/Noncoding_Assoc.qassoc", header=TRUE)
yld.syn <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/PLINK_qAssoc/Yield/Synonymous/Synonymous_Assoc.qassoc", header=TRUE)
yld.nonsyn <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/PLINK_qAssoc/Yield/Nonsynonymous/Nonsynonymous_Assoc.qassoc", header=TRUE)
yld.del <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/PLINK_qAssoc/Yield/Deleterious/Deleterious_Assoc.qassoc", header=TRUE)

don.nc <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/PLINK_qAssoc/DON/Noncoding/Noncoding_Assoc.qassoc", header=TRUE)
don.syn <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/PLINK_qAssoc/DON/Synonymous/Synonymous_Assoc.qassoc", header=TRUE)
don.nonsyn <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/PLINK_qAssoc/DON/Nonsynonymous/Nonsynonymous_Assoc.qassoc", header=TRUE)
don.del <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/PLINK_qAssoc/DON/Deleterious/Deleterious_Assoc.qassoc", header=TRUE)

pht.nc <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/PLINK_qAssoc/Height/Noncoding/Noncoding_Assoc.qassoc", header=TRUE)
pht.syn <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/PLINK_qAssoc/Height/Synonymous/Synonymous_Assoc.qassoc", header=TRUE)
pht.nonsyn <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/PLINK_qAssoc/Height/Nonsynonymous/Nonsynonymous_Assoc.qassoc", header=TRUE)
pht.del <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/PLINK_qAssoc/Height/Deleterious/Deleterious_Assoc.qassoc", header=TRUE)


# Make plots for the ECDF
pdf(file="qAssoc_ECDF.pdf", height=4, width=12)
par(mfrow=c(1, 3), mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
plot(
    ecdf(yld.nc$P),
    xlab="P-value",
    ylab="Proportion of SNPs",
    main="\nYield",
    col="black",
    lwd=2)
lines(ecdf(yld.syn$P), lwd=2, col="blue")
lines(ecdf(yld.nonsyn$P), lwd=2, col="green")
lines(ecdf(yld.del$P), lwd=2, col="red")

plot(
    ecdf(don.nc$P),
    xlab="P-value",
    ylab="Proportion of SNPs",
    main="\nDON",
    col="black",
    lwd=2)
lines(ecdf(don.syn$P), lwd=2, col="blue")
lines(ecdf(don.nonsyn$P), lwd=2, col="green")
lines(ecdf(don.del$P), lwd=2, col="red")

plot(
    ecdf(pht.nc$P),
    xlab="P-value",
    ylab="Proportion of SNPs",
    main="\nHeight",
    col="black",
    lwd=2)
lines(ecdf(pht.syn$P), lwd=2, col="blue")
lines(ecdf(pht.nonsyn$P), lwd=2, col="green")
lines(ecdf(pht.del$P), lwd=2, col="red")

legend("bottomright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"),
    col=c("black", "blue", "green", "red"), lwd=2)
dev.off()

pdf(file="qAssoc_ECDF_Zoom.pdf", height=4, width=12)
par(mfrow=c(1, 3), mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
plot(
    ecdf(yld.nc$P),
    xlab="P-value",
    ylab="Proportion of SNPs",
    main="\nYield",
    col="black",
    lwd=2,
    xlim=c(0, 0.05))
lines(ecdf(yld.syn$P), lwd=2, col="blue")
lines(ecdf(yld.nonsyn$P), lwd=2, col="green")
lines(ecdf(yld.del$P), lwd=2, col="red")

plot(
    ecdf(don.nc$P),
    xlab="P-value",
    ylab="Proportion of SNPs",
    main="\nDON",
    col="black",
    lwd=2,
    xlim=c(0, 0.05))
lines(ecdf(don.syn$P), lwd=2, col="blue")
lines(ecdf(don.nonsyn$P), lwd=2, col="green")
lines(ecdf(don.del$P), lwd=2, col="red")

plot(
    ecdf(pht.nc$P),
    xlab="P-value",
    ylab="Proportion of SNPs",
    main="\nHeight",
    col="black",
    lwd=2,
    xlim=c(0, 0.05))
lines(ecdf(pht.syn$P), lwd=2, col="blue")
lines(ecdf(pht.nonsyn$P), lwd=2, col="green")
lines(ecdf(pht.del$P), lwd=2, col="red")

legend("bottomright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"),
    col=c("black", "blue", "green", "red"), lwd=2)
dev.off()

pdf(file="qAssoc_Density.pdf", height=4, width=12)
par(mfrow=c(1, 3), mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
plot(
    density(yld.nc$P, na.rm=TRUE),
    xlab="P-value",
    ylab="Density",
    main="\nYield",
    col="black",
    lwd=2)
lines(density(yld.syn$P, na.rm=TRUE), lwd=2, col="blue")
lines(density(yld.nonsyn$P, na.rm=TRUE), lwd=2, col="green")
lines(density(yld.del$P, na.rm=TRUE), lwd=2, col="red")

plot(
    density(don.nc$P, na.rm=TRUE),
    xlab="P-value",
    ylab="Density",
    main="\nDON",
    col="black",
    lwd=2)
lines(density(don.syn$P, na.rm=TRUE), lwd=2, col="blue")
lines(density(don.nonsyn$P, na.rm=TRUE), lwd=2, col="green")
lines(density(don.del$P, na.rm=TRUE), lwd=2, col="red")

plot(
    density(pht.nc$P, na.rm=TRUE),
    xlab="P-value",
    ylab="Density",
    main="\nHeight",
    col="black",
    lwd=2)
lines(density(pht.syn$P, na.rm=TRUE), lwd=2, col="blue")
lines(density(pht.nonsyn$P, na.rm=TRUE), lwd=2, col="green")
lines(density(pht.del$P, na.rm=TRUE), lwd=2, col="red")

legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"),
    col=c("black", "blue", "green", "red"), lwd=2)
dev.off()
