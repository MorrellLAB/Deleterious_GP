# Make a plot of the number of deleterious SNPs per sequenced codon in 1Mb
# windows across the genome.

# read in the data
dsnp <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/del_per_codon.txt", header=TRUE)
nssnp <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/nonsyn_per_codon.txt", header=TRUE)

summary(dsnp)
summary(nssnp)

ns_col <- '#fdae61'
del_col <- '#d7191c'

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
pdf(file="Prop_dSNP_nsSNP_Per_Codon.pdf", height=10.5, width=8)
par(mfrow=c(7, 1), mar=c(4, 4, 1, 2))
for(chrom in chroms) {
    cdat <- dsnp[dsnp$chr == chrom,]
    ndat <- nssnp[nssnp$chr == chrom,]
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
        ylab="SNPs per Codon",
        main=chrom,
        ylim=c(1e-5, 0.02),
        xlim=c(0, 800),
        axes=F,
        type="n",
        log="y")
    rect(pcent$Start, 1e-5, pcent$End, 1, density=NA, col=rgb(0, 0, 0, alpha=0.2))
    rect(cent$Start, 1e-5, cent$End, 1, density=NA, col=rgb(0, 0, 0, alpha=0.3))
    ndat$delSNP_codonNb[ndat$delSNP_codonNb < 1e-5] <- 1e-5
    cdat$delSNP_codonNb[cdat$delSNP_codonNb < 1e-5] <- 1e-5
    points(ndat$delSNP_codonNb ~ midp, pch=19, col=ns_col)
    points(cdat$delSNP_codonNb ~ midp, pch=19, col=del_col)
    axis(side=2)
    axis(side=1, at=seq(0, 800, by=50), labels=seq(0, 800, by=50))
    if(chrom == "chr1H") {
        legend("topright", c("Tolerated", "Deleterious"), col=c(ns_col, del_col), pch=19, ncol=2)
    }
}
dev.off()
