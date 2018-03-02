# Make a plot of the number of deleterious SNPs per sequenced codon in 1Mb
# windows across the genome.

# read in the data
dat <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Prop_Deleterious_Per_Codon.txt.gz", header=TRUE)

# The chromosomes to plot
chroms <- c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H")
# And the centromeres
centromeres <- data.frame(
    Chr=c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H"),
    Start=c(NA, 136904248, 115794456, 58804724, 60207685, 180782470, 202235859),
    End=c(NA, 528069469, 438235281, 467767377, 368181264, 383275592, 449173785))

# Make a PDF
pdf(file="Prop_Del_Per_Codon.pdf", height=10.5, width=8)
par(mfrow=c(7, 1), mar=c(4, 4, 1, 2))
for(chrom in chroms) {
    cdat <- dat[dat$chr == chrom,]
    midp <- (cdat$starting + cdat$ending)/2
    cent <- centromeres[centromeres$Chr == chrom,]
    # Convert from bp to Mbp
    midp <- midp/1000000
    cent$Start <- cent$Start/1000000
    cent$End <- cent$End/1000000
    plot(
        cdat$delSNP_codonNb ~ midp,
        xlab="Position (Mb)",
        ylab="DSNPs per Codon",
        main=chrom,
        ylim=c(0, 0.01),
        xlim=c(0, 800),
        cex=0.8,
        axes=F)
    rect(cent$Start, 0, cent$End, 1, density=NA, col=rgb(0, 0, 0, alpha=0.2))
    axis(side=2)
    axis(side=1, at=seq(0, 800, by=50), labels=seq(0, 800, by=50))
}
dev.off()
