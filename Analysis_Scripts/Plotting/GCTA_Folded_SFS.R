#   Script to plot the folded site frequency spectra for the full population
#   used for the genomic prediction experiment. Will plot folded SFS for
#   All, Noncoding, Coding, Synonymous, BOPA SNPs.

#   Set working directory
setwd("/Users/tomkono/DataDisk/Dropbox/Projects/DM_GenomicPrediction/Data/Progeny_Genotypes/GCTA/MAFs")

#   Read the MAFs from the PLINK outputs
all.maf <- read.table("all.frq", header=T)
noncoding.maf <- read.table("noncoding.frq", header=T)
coding.maf <- read.table("coding.frq", header=T)
nonsyn.maf <- read.table("nonsyn.frq", header=T)
del.maf <- read.table("deleterious.frq", header=T)

#   Define a function that returns a SFS, given a vector of frequencies, a
#   frequency into which to bin, and whether or not the SFS is folded.
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

#   Generate SFS for all the partitions
all.sfs <- make_sfs(all.maf$MAF, 0.05, folded=TRUE)
noncoding.sfs <- make_sfs(noncoding.maf$MAF, 0.05, folded=TRUE)
coding.sfs <- make_sfs(coding.maf$MAF, 0.05, folded=TRUE)
nonsyn.sfs <- make_sfs(nonsyn.maf$MAF, 0.05, folded=TRUE)
del.sfs <- make_sfs(del.maf$MAF, 0.05, folded=TRUE)

#   And put them into a data frame for plotting
plot.sfs.data <- as.data.frame(
    cbind(
        all.sfs,
        noncoding.sfs,
        coding.sfs,
        nonsyn.sfs,
        del.sfs
        )
    )

#   Make a plot
pdf(file="Imputed_MAFs.pdf", 8, 6)
plt <- barplot(
    t(plot.sfs.data),
    ylim=c(0, 0.4),
    beside=TRUE,
    axisnames=F,
    xlab="Minor Allele Frequency",
    ylab="Proportion",
    col=c("black", "grey70", "blue", "green", "red")
    )
#   Get the labels of the frequency classes as the names of the output from
#   make_sfs()
labels <- names(all.sfs)
#   Calculate where the labels should be, as the midpoint (mean) of the bar
#   group positions.
at <- apply(plt, 2, mean)
#   And place the labels
axis(
    side=1,
    at=at,
    labels=labels,
    font=1,
    cex.axis=0.75
    )
#   Draw a legend
legend(
  "topright",
  c("All", "Noncoding", "Coding", "Nonsynonymous", "Deleterious"),
  fill=c("black", "grey70", "blue", "green", "red"),
  cex=1
  )
dev.off()
