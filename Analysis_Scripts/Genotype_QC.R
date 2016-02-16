#   Script to create genotype quality and genetic map quality plots for the
#   three cycles of progeny genotypes.

#   To take arguments
args <- commandArgs(TRUE)

#   We need the r/qtl package for all our genotype quality checks
library(qtl)

#   Read in the three cycles' data. These are in 'csvr' format, which is
#   documented in the r/qtl manual here:
#       http://www.rqtl.org/sampledata/
#   The default cross type is an F2. Let's compare F2 to F3
cycle1.f2 <- read.cross(
    "csvr",
    file=args[1],
    genotypes=c(1, 2, 3))
cycle2.f2 <- read.cross(
    "csvr",
    file=args[2],
    genotypes=c(1, 2, 3))
cycle3.f2 <- read.cross(
    "csvr",
    file=args[3],
    genotypes=c(1, 2, 3))
#   Read in as F3s
cycle1.f3 <- read.cross(
    "csvr",
    file=args[1],
    genotypes=c(1, 2, 3),
    F.gen=3,
    BC.gen=0)
cycle2.f3 <- read.cross(
    "csvr",
    file=args[2],
    genotypes=c(1, 2, 3),
    F.gen=3,
    BC.gen=0)
cycle3.f3 <- read.cross(
    "csvr",
    file=args[3],
    genotypes=c(1, 2, 3),
    F.gen=3,
    BC.gen=0)

#   Plots of missingness should be the same for both cross types. Just use the
#   F2 object.
pdf(file="Cycle1_Missingness.pdf", 8, 8)
plotMissing(cycle1.f2)
dev.off()
pdf(file="Cycle2_Missingness.pdf", 8, 8)
plotMissing(cycle2.f2)
dev.off()
pdf(file="Cycle3_Missingness.pdf", 8, 8)
plotMissing(cycle3.f2)
dev.off()

#   Next, we will estimate maps from F2 and F3 and compare them
c1.f2.map <- est.map(cycle1.f2)
c1.f3.map <- est.map(cycle1.f3)
pdf(file="Cycle1_MapComparison.pdf", 8, 8)
plotMap(cycle1.f2, c1.f2.map, main="F2 and F3 Map Comparison for Cycle 1")
dev.off()
