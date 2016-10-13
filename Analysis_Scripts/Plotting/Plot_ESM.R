#   Plot the ESM test statistic across each of the seven barley chromosomes

library(ggplot2)
setwd("/Volumes/Data_Disk/tmp/Perms")

#   Read in each of the ESM output files. There is one per chromosome.
esm.1H <- read.table("1H/1H_ESM.txt", header=TRUE)
esm.2H <- read.table("2H/2H_ESM.txt", header=TRUE)
esm.3H <- read.table("3H/3H_ESM.txt", header=TRUE)
esm.4H <- read.table("4H/4H_ESM.txt", header=TRUE)
esm.5H <- read.table("5H/5H_ESM.txt", header=TRUE)
esm.6H <- read.table("6H/6H_ESM.txt", header=TRUE)
esm.7H <- read.table("7H/7H_ESM.txt", header=TRUE)

#   Make a large data frame for plotting
esm.data <- data.frame(
    Chromosome=c(
        rep("chr1H", nrow(esm.1H)),
        rep("chr2H", nrow(esm.2H)),
        rep("chr3H", nrow(esm.3H)),
        rep("chr4H", nrow(esm.4H)),
        rep("chr5H", nrow(esm.5H)),
        rep("chr6H", nrow(esm.6H)),
        rep("chr7H", nrow(esm.7H))),
    P.value=c(
        esm.1H$p.values,
        esm.2H$p.values,
        esm.3H$p.values,
        esm.4H$p.values,
        esm.5H$p.values,
        esm.6H$p.values,
        esm.7H$p.values),
    Position=c(
        esm.1H$loci.midpoint,
        esm.2H$loci.midpoint,
        esm.3H$loci.midpoint,
        esm.4H$loci.midpoint,
        esm.5H$loci.midpoint,
        esm.6H$loci.midpoint,
        esm.7H$loci.midpoint)
    )

#   How many tests were there, genome-wide?
ntests <- nrow(esm.data)
#   Set a significance threshold
sig <- 0.05/ntests
sig <- -log(sig, base=10)
print(c(ntests, sig))

#   And plot the data.
pdf(file="ESM_test.pdf", width=10.5, height=8)

ggplot(esm.data, aes(x=Position/1000000, y=-log(P.value, base=10))) +
    geom_line(size=0.25) +
    facet_grid(Chromosome~.) +
    theme_bw() +
    theme(
        strip.background=element_blank(),
        strip.text.y=element_text(size=10, colour="black", angle=0),
        axis.ticks.y=element_blank()) +
    scale_y_continuous(limits=c(0, 4), breaks=c(0, 2, 4)) +
    scale_x_continuous(limits=c(0, 540), breaks=seq(0, 540, by=50)) +
    labs(x="Physical Position (Mb)", y="-log(P value)", title="ESM Test")
dev.off()
