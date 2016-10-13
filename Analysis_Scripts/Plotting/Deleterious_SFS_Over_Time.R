#   Plot a derived site frequency spectrum for C1, C2, and C3 for various
#   partitions of SNPs.

library(ggplot2)
library(reshape)

setwd("/Volumes/Data_Disk/Dropbox/Projects/DM_GenomicPrediction/Results/DM_Selection")

make_sfs <- function(freqs, binsize, folded=TRUE) {
    #   Set the maximum frequency. If folded=TRUE, it is 0.5, else it is
    #   1.0
    if(folded) {
        maxfreq <- 0.5
    } else {
        maxfreq <- 1.0
    }
    #   Generate the bins
    bins <- seq(0.0, maxfreq, by=binsize)
    #   Divide the frequencies into the bins
    binned_freqs <- cut(freqs, breaks=bins, include.lowest=TRUE)
    #   Count how many are in each bin
    binned_freqs <- table(binned_freqs)
    #   And convert them into proportions
    binned_freqs <- binned_freqs / length(freqs)
    #   Return it.
    return(binned_freqs)
}

#   Read in the frequencies
syn.freq <- read.table("Syn_DAFs.txt", header=T)
tol.freq <- read.table("Tol_DAFs.txt", header=T)
del.freq <- read.table("Del_DAFs.txt", header=T)

syn.sfs <- apply(syn.freq, 2, make_sfs, 0.05, FALSE)
tol.sfs <- apply(tol.freq, 2, make_sfs, 0.05, FALSE)
del.sfs <- apply(del.freq, 2, make_sfs, 0.05, FALSE)

 #  Group them by cycle first, by functional class second
comb.sfs <- data.frame(
    Cycle=c(rep("C1", nrow(syn.sfs)), rep("C2", nrow(syn.sfs)), rep("C3", nrow(syn.sfs))),
    DAF=rep(row.names(syn.sfs), 3),
    Synonymous=c(syn.sfs[,"C1_DAF"], syn.sfs[,"C2_DAF"], syn.sfs[,"C3_DAF"]),
    Tolerated=c(tol.sfs[,"C1_DAF"], tol.sfs[,"C2_DAF"], tol.sfs[,"C3_DAF"]),
    Deleterious=c(del.sfs[,"C1_DAF"], del.sfs[,"C2_DAF"], del.sfs[,"C3_DAF"])
    )

#   Melt it for ggplot
comb.sfs <- melt(comb.sfs, id=c("Cycle", "DAF"))

#   Plot it!
pdf(file="DM_Selection.pdf", width=6, height=10)
ggplot(comb.sfs, aes(x=factor(DAF, levels=row.names(syn.sfs)), y=value, fill=factor(variable))) +
    geom_bar(stat="identity", position="dodge") +
    facet_grid(Cycle~.) +
    scale_fill_manual(
        name="",
        breaks=c("Synonymous", "Tolerated", "Deleterious"),
        labels=c("Synonymous", "Tolerated", "Deleterious"),
        values=c("black", "blue", "red")) +
    theme_bw() +
    theme(
        strip.background=element_blank(),
        strip.text.y=element_text(size=10, colour="black", angle=0),
        axis.text.x=element_text(size=8, colour="black", angle=90, hjust=1, vjust=0.5)) +
    labs(y="Proportion", x="Derived Allele Frequency")
dev.off()
