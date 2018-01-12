# Plot the observed heterozygosity in the imputed SNPs from the AlphaPeel
# output.

# Take arguments
args <- commandArgs(TRUE)

# Make an output filename
outf <- gsub(".txt", "_Ho.pdf", args[1])

dat <- read.table(args[1], header=TRUE)
# Isolate the cycles
cycle1 <- grep("MS10", dat$ID)
cycle2 <- grep("MS11", dat$ID)
cycle3 <- grep("MS12", dat$ID)

cycle1 <- dat[cycle1,]
cycle2 <- dat[cycle2,]
cycle3 <- dat[cycle3,]

# Calculate the proportion of heterozygous sites relative to the number of
# polymorphic sites in each family.
c1.ho <- cycle1$NHet/cycle1$NSeg
c2.ho <- cycle2$NHet/cycle2$NSeg
c3.ho <- cycle3$NHet/cycle3$NSeg

# remove NAs
c1.ho <- c1.ho[!is.na(c1.ho)]
c2.ho <- c2.ho[!is.na(c2.ho)]
c3.ho <- c3.ho[!is.na(c3.ho)]

pdf(file=outf, height=4, width=6)
plot(
    density(c1.ho),
    col="black",
    lwd=2,
    xlab="Proportion of Heterozygous Sites",
    ylab="Density",
    main="Observed Heterozygosity at Segregating Sites")
lines(density(c2.ho), col="red", lwd=2)
lines(density(c3.ho), col="blue", lwd=2)
legend(
    "topright",
    c("Cycle 1", "Cycle 2", "Cycle 3"),
    col=c("black", "red", "blue"),
    lwd=2)
dev.off()
