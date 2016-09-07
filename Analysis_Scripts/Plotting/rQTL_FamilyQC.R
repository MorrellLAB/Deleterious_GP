#   Uses r/QTL to check the progeny genotyping data by family. This script will
#   output diagnostic plots for missingness, and the genetic positions of the
#   typed markers for each family. It will eventually produce plots for
#   potential genotyping errors, too.

library(qtl)
args <- commandArgs(TRUE)

all_chr_names <- c("1H", "2H", "3H", "4H", "5H", "6H", "7H")
#   First, read in the cross object as an F2. We want to at least plot the
#   genotype distribution across the chromosomes.
pop_pre <- read.cross(
    format="csv",
    file=args[1])
genoplotname <- gsub(".csv", "_PreQC_PlotGeno.pdf", args[1])
pdf(
    file=genoplotname,
    height=15,
    width=15)
par(mfrow=c(4, 2))
par(mfrow=c(4, 2))
for(chr in all_chr_names) {
    if(chr %in% chrnames(pop_pre)) {
        if(length(pop_pre[["geno"]][[chr]][["map"]]) == 1) {
            plot(0, 0, type="n", main=paste("Chromosome ", chr))
        }
        else {
            plotGeno(pop_pre, chr=chr)
        }
    }
    else {
        plot(0, 0, type="n", main=paste("Chromosome ", chr))
    }
}
dev.off()


#   In our case, the cross type will be of 'bcsft' - meaning s generation of
#   backcrossing, and t generations of self fertilization. We use s=0 and t=3
#   for an F3.
pop <- read.cross(
    format="csv",
    file=args[1],
    BC.gen=0,
    F.gen=3)
#   Print out a nice summary of the family
summary(pop)

#   Plot the missingness for the family
m_plotname <- gsub(".csv", "_Missingness.pdf", args[1])
pdf(
    file=m_plotname,
    height=6,
    width=6)
plotMissing(pop)
dev.off()


#   Remove markers with high missingness (>50%)
non_missing <- ntyped(pop, "mar")
todrop_missing <- names(non_missing[non_missing < 0.5*nind(pop)])
pop <- drop.markers(pop, todrop_missing)
print(todrop_missing)

#   Remove markers with bad segregation patterns. These are largely
#   uninformative, or mis-informative. First, tabulate the genotypes for each
#   marker, and find those with distorted segregation patterns
geno_table <- geno.table(pop)
#   Drop them. We use 5%, corrected for multiple testing
todrop_segdist <- rownames(geno_table[geno_table$P.value < 0.05/totmar(pop),])
pop <- drop.markers(pop, todrop_segdist)
#   Print out which markers we should drop
print(todrop_segdist)

#   Identify markers with identical genotypes. Drop one from each pair
dup_mar <- findDupMarkers(pop, exact.only=FALSE)
pop <- drop.markers(pop, unlist(dup_mar))
print(dup_mar)


#   Then plot it after QC
genoplotname <- gsub(".csv", "_PostQC_PlotGeno.pdf", args[1])
pdf(
    file=genoplotname,
    height=15,
    width=15)
par(mfrow=c(4, 2))
for(chr in all_chr_names) {
    if(chr %in% chrnames(pop)) {
        if(length(pop[["geno"]][[chr]][["map"]]) == 1) {
            plot(0, 0, type="n", main=paste("Chromosome ", chr))
        }
        else {
            plotGeno(pop, chr=chr, include.xo=FALSE)
        }
    }
    else {
        plot(0, 0, type="n", main=paste("Chromosome ", chr))
    }
}
dev.off()


#   Plot the map pre- and post- QC
map_plotname <- gsub(".csv", "_Map.pdf", args[1])
pdf(
    file=map_plotname,
    height=6,
    width=6)
par(mfrow=c(1, 2))
plotMap(pop_pre)
plotMap(pop)
dev.off()
