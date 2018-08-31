# Make a plot of the phenotype (YLD and DON) versus the dosage of deleterious
# alleles for C1, C2, and C3 progeny.

nc_col <- '#2c7bb6'
syn_col <- '#abd9e9'
ns_col <- '#fdae61'
del_col <- '#d7191c'

# Read the phenotype data
phen <- read.csv("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Data/Phenotypic_Data/Adjusted_Phenotypic_Data_800Lines.csv", header=TRUE)
# Deleterious allele dosages
dosages <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Burden/Homozygous_Derived_Counts.txt.gz", header=TRUE)

# Merge the tables
merged <- merge(phen, dosages, by.x="line", by.y="LineID", all.x=TRUE, all.y=FALSE)


cor.test(merged$yld_kg, merged$Noncoding)
cor.test(merged$DON, merged$Noncoding)
cor.test(merged$height, merged$Noncoding)

cor.test(merged$yld_kg, merged$Synonymous)
cor.test(merged$DON, merged$Synonymous)
cor.test(merged$height, merged$Synonymous)

cor.test(merged$yld_kg, merged$Nonsynonymous)
cor.test(merged$DON, merged$Nonsynonymous)
cor.test(merged$height, merged$Nonsynonymous)

cor.test(merged$yld_kg, merged$Deleterious)
cor.test(merged$DON, merged$Deleterious)
cor.test(merged$height, merged$Deleterious)


# We will separate the cycles for visualization purposes
c0 <- merged[merged$cycle == "C0",]
c1 <- merged[merged$cycle == "C1",]
c2 <- merged[merged$cycle == "C2",]
c3 <- merged[merged$cycle == "C3",]


# Plot the distributions
pdf(file="Deleterious_Phenotype_Correlations.pdf", height=8, width=8)
par(mfrow=c(2, 2), mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
plot(
    c0$yld_kg ~ c0$Deleterious,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Homozygous Derived Deleterious SNPs",
    ylab="Yield BLUE (kg/ha)",
    main="\nYield",
    xlim=c(110, 300),
    ylim=c(3500, 8000))
points(c1$yld_kg ~ c1$Deleterious, pch=16, col="blue", cex=0.6)
points(c2$yld_kg ~ c2$Deleterious, pch=17, col="green", cex=0.6)
points(c3$yld_kg ~ c3$Deleterious, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))

plot(
    c0$DON ~ c0$Deleterious,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Homozygous Derived Deleterious SNPs",
    ylab="DON BLUE (ppm)",
    main="\nDON Concentration",
    xlim=c(110, 300),
    ylim=c(0, 35))
points(c1$DON ~ c1$Deleterious, pch=16, col="blue", cex=0.6)
points(c2$DON ~ c2$Deleterious, pch=17, col="green", cex=0.6)
points(c3$DON ~ c3$Deleterious, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))

plot(
    c0$height ~ c0$Deleterious,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Homozygous Derived Deleterious SNPs",
    ylab="Height BLUE (cm)",
    main="\nHeight",
    xlim=c(110, 300),
    ylim=c(50, 110))
points(c1$height ~ c1$Deleterious, pch=16, col="blue", cex=0.6)
points(c2$height ~ c2$Deleterious, pch=17, col="green", cex=0.6)
points(c3$height ~ c3$Deleterious, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))
dev.off()


pdf(file="Nonsynonymous_Phenotype_Correlations.pdf", height=8, width=8)
par(mfrow=c(2, 2), mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
plot(
    c0$yld_kg ~ c0$Nonsynonymous,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Homozygous Derived Nonsynonymous SNPs",
    ylab="Yield BLUE (kg/ha)",
    main="\nYield",
    xlim=c(11000, 15000),
    ylim=c(3500, 8000))
points(c1$yld_kg ~ c1$Nonsynonymous, pch=16, col="blue", cex=0.6)
points(c2$yld_kg ~ c2$Nonsynonymous, pch=17, col="green", cex=0.6)
points(c3$yld_kg ~ c3$Nonsynonymous, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))

plot(
    c0$DON ~ c0$Nonsynonymous,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Homozygous Derived Nonsynonymous SNPs",
    ylab="DON BLUE (ppm)",
    main="\nDON Concentration",
    xlim=c(11000, 15000),
    ylim=c(0, 35))
points(c1$DON ~ c1$Nonsynonymous, pch=16, col="blue", cex=0.6)
points(c2$DON ~ c2$Nonsynonymous, pch=17, col="green", cex=0.6)
points(c3$DON ~ c3$Nonsynonymous, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))

plot(
    c0$height ~ c0$Nonsynonymous,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Homozygous Derived Nonsynonymous SNPs",
    ylab="Height BLUE (cm)",
    main="\nHeight",
    xlim=c(11000, 15000),
    ylim=c(50, 110))
points(c1$height ~ c1$Nonsynonymous, pch=16, col="blue", cex=0.6)
points(c2$height ~ c2$Nonsynonymous, pch=17, col="green", cex=0.6)
points(c3$height ~ c3$Nonsynonymous, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))
dev.off()

pdf(file="Synonymous_Phenotype_Correlations.pdf", height=8, width=8)
par(mfrow=c(2, 2), mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
plot(
    c0$yld_kg ~ c0$Synonymous,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Homozygous Derived Synonymous SNPs",
    ylab="Yield BLUE (kg/ha)",
    main="\nYield",
    xlim=c(14000, 19000),
    ylim=c(3500, 8000))
points(c1$yld_kg ~ c1$Synonymous, pch=16, col="blue", cex=0.6)
points(c2$yld_kg ~ c2$Synonymous, pch=17, col="green", cex=0.6)
points(c3$yld_kg ~ c3$Synonymous, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))

plot(
    c0$DON ~ c0$Synonymous,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Homozygous Derived Synonymous SNPs",
    ylab="DON BLUE (ppm)",
    main="\nDON Concentration",
    xlim=c(14000, 19000),
    ylim=c(0, 35))
points(c1$DON ~ c1$Synonymous, pch=16, col="blue", cex=0.6)
points(c2$DON ~ c2$Synonymous, pch=17, col="green", cex=0.6)
points(c3$DON ~ c3$Synonymous, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))

plot(
    c0$height ~ c0$Synonymous,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Homozygous Derived Synonymous SNPs",
    ylab="Height BLUE (cm)",
    main="\nHeight",
    xlim=c(14000, 19000),
    ylim=c(50, 110))
points(c1$height ~ c1$Synonymous, pch=16, col="blue", cex=0.6)
points(c2$height ~ c2$Synonymous, pch=17, col="green", cex=0.6)
points(c3$height ~ c3$Synonymous, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))
dev.off()


pdf(file="Noncoding_Phenotype_Correlations.pdf", height=8, width=8)
par(mfrow=c(2, 2), mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
plot(
    c0$yld_kg ~ c0$Noncoding,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Homozygous Derived Noncoding SNPs",
    ylab="Yield BLUE (kg/ha)",
    main="\nYield",
    xlim=c(55000, 67000),
    ylim=c(3500, 8000))
points(c1$yld_kg ~ c1$Noncoding, pch=16, col="blue", cex=0.6)
points(c2$yld_kg ~ c2$Noncoding, pch=17, col="green", cex=0.6)
points(c3$yld_kg ~ c3$Noncoding, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))

plot(
    c0$DON ~ c0$Noncoding,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Homozygous Derived Noncoding SNPs",
    ylab="DON BLUE (ppm)",
    main="\nDON Concentration",
    xlim=c(55000, 67000),
    ylim=c(0, 35))
points(c1$DON ~ c1$Noncoding, pch=16, col="blue", cex=0.6)
points(c2$DON ~ c2$Noncoding, pch=17, col="green", cex=0.6)
points(c3$DON ~ c3$Noncoding, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))

plot(
    c0$height ~ c0$Noncoding,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Homozygous Derived Noncoding SNPs",
    ylab="Height BLUE (cm)",
    main="\nHeight",
    xlim=c(55000, 67000),
    ylim=c(50, 110))
points(c1$height ~ c1$Noncoding, pch=16, col="blue", cex=0.6)
points(c2$height ~ c2$Noncoding, pch=17, col="green", cex=0.6)
points(c3$height ~ c3$Noncoding, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))
dev.off()
