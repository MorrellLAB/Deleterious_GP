#   Uses r/QTL to check the progeny genotyping data by family. This script will
#   output diagnostic plots for missingness, and the genetic positions of the
#   typed markers for each family. It will eventually produce plots for
#   potential genotyping errors, too.

library(qtl)
args <- commandArgs(TRUE)

#   First, read in the cross object. In our case, the cross type will be of
#   'bcsft' - meaning s generation of backcrossing, and t generations of self
#   fertilization. We use s=0 and t=3 for an F3.
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

#   Plot the map comparison
estmap_plotname <- gsub(".csv", "_EstMap.pdf", args[1])
pdf(
    file=estmap_plotname,
    height=6,
    width=6)
#   Remove markers with bad segregation patterns. These are largely
#   uninformative, or mis-informative. First, tabulate the genotypes for each
#   marker, and find those with distorted segregation patterns
geno_table <- geno.table(pop)
#   Drop them. We use 5%, corrected for multiple testing
todrop <- rownames(geno_table[geno_table$P.value < 0.05/totmar(pop),])
pop <- drop.markers(pop, todrop)
#   For now, we go with the Morgan map function
newmap <- est.map(pop, map.function="morgan")
plotMap(pop, newmap)
dev.off()

#   Plot the lrecombination fractions, and the LOD scores for association (the
#   LOD that r = 1/2 v linked)
rf_plotname <- gsub(".csv", "_RF.pdf", args[1])
pdf(
    file=rf_plotname,
    height=6,
    width=6)
pop <- est.rf(pop)
plotRF(pop)
dev.off()

#   Does it look like there are markers with switched alleles?
rflod_plotname <- gsub(".csv", "_RF_LOD.pdf", args[1])
pdf(
    file=rflod_plotname,
    height=6,
    width=6)
rf <- pull.rf(pop)
lod <- pull.rf(pop, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recomb. Frac.", ylab="LOD score")
dev.off()
