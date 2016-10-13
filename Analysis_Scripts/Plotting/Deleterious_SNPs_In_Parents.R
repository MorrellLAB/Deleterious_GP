#   Script to produce a beeswarm plot of the number of deleterious SNPs carried
#   by each of the parents.

library(beeswarm)
dac <- read.table("/Volumes/Data_Disk/Dropbox/Projects/DM_GenomicPrediction/Results/Annotations/Derived_Allele_Counts.txt", header=T)

pdf(
    file="Parental_Deleterious_Count.pdf",
    height=6,
    width=3)
beeswarm(dac$Count, pwcol=dac$BP, pch=19, ylim=c(0, 300))
bxplot(dac$Count, add=TRUE)
legend("topright", col=1:3, pch=19, c("Busch Ag.", "UMN", "NDSU"))
dev.off()
