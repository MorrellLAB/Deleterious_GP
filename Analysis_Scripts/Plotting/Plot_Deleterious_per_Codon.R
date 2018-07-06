# Make a plot of the number of deleterious SNPs per sequenced codon in 1Mb
# windows across the genome.

# read in the data
dat <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Prop_Deleterious_Per_Codon.txt", header=TRUE)

# The chromosomes to plot
chroms <- c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H")
# And the pericentromeric regions. These are estimated from Figure 2 of
# Munoz et al. 2015 on the basis of recombination rate.
pericentromeres <- data.frame(
    Chr=c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H"),
    Start=c(70000000, 140000000, 140000000, 50000000, 50000000, 80000000, 175000000),
    End=c(240000000, 425000000, 315000000, 325000000, 260000000, 325000000, 375000000))
# The actual centromeres are from Mascher et al. 2017
centromeres <- data.frame(
    Chr=c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H"),
    Start=c(213600001, 333600001, 288000001, 286400001, 217600001, 257600001, 339200001),
    End=c(221600000, 341600000, 299200000, 296000000, 225600000, 265600000, 345600000))

# Make a PDF
pdf(file="Prop_Del_Per_Codon.pdf", height=10.5, width=8)
par(mfrow=c(7, 1), mar=c(4, 4, 1, 2))
for(chrom in chroms) {
    cdat <- dat[dat$chr == chrom,]
    midp <- (cdat$starting + cdat$ending)/2
    cent <- centromeres[centromeres$Chr == chrom,]
    pcent <- pericentromeres[pericentromeres$Chr == chrom,]
    # Convert from bp to Mbp
    midp <- midp/1000000
    pcent$Start <- pcent$Start/1000000
    pcent$End <- pcent$End/1000000
    cent$Start <- cent$Start/1000000
    cent$End <- cent$End/1000000
    plot(
        cdat$delSNP_codonNb ~ midp,
        xlab="Position (Mb)",
        ylab="dSNPs per Codon",
        main=chrom,
        ylim=c(0, 0.125),
        xlim=c(0, 800),
        axes=F,
        type="n")
    rect(pcent$Start, 0, pcent$End, 1, density=NA, col=rgb(0, 0, 0, alpha=0.2))
    rect(cent$Start, 0, cent$End, 1, density=NA, col=rgb(0, 0, 0, alpha=0.3))
    points(cdat$delSNP_codonNb ~ midp, pch=19)
    axis(side=2)
    axis(side=1, at=seq(0, 800, by=50), labels=seq(0, 800, by=50))
}
dev.off()
