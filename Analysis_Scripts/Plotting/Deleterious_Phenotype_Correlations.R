# Make a plot of the phenotype (YLD and DON) versus the dosage of deleterious
# alleles for C1, C2, and C3 progeny.

nc_col <- '#2c7bb6'
syn_col <- '#abd9e9'
ns_col <- '#fdae61'
del_col <- '#d7191c'

# Read the phenotype data
phen <- read.csv("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Data/Phenotypic_Data/Adjusted_Phenotypic_Data.csv", header=TRUE)
# Deleterious allele dosages
dosages <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Summaries/GP_Deleterious_Homozygotes.txt.gz", header=TRUE)
dosages <- dosages[!dosages$line_name=="CELEBRATION",]
# Nonsynonymous dosages
nonsyn_dosages <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Summaries/GP_Nonsynonymous_Homozygotes.txt.gz", header=TRUE)
nonsyn_dosages <- nonsyn_dosages[!nonsyn_dosages$line_name == "CELEBRATION",]
syn_dosages <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Summaries/GP_Synonymous_Homozygotes.txt.gz", header=TRUE)
syn_dosages <- syn_dosages[!nonsyn_dosages$line_name == "CELEBRATION",]
# Pedigree file that gives sel/ran
sel <- read.csv("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Data/Pedigrees/Population_Pedigrees.csv", header=TRUE)

# Trim down the sel file to just line name and selection status
sel <- sel[, c("Line", "Type")]
# And change the names to match the other files
names(sel) <- c("line_name", "Type")

# Merge the tables
m_tmp <- merge(phen, dosages, by="line_name")
merged <- merge(m_tmp, sel, by="line_name", all.x=TRUE)

ns_m_tmp <- merge(phen, nonsyn_dosages, by="line_name")
ns_merged <- merge(ns_m_tmp, sel, by="line_name", all.x=TRUE)

s_m_tmp <- merge(phen, syn_dosages, by="line_name")
s_merged <- merge(s_m_tmp, sel, by="line_name", all.x=TRUE)

cor.test(merged$yld.BLUE_mv, merged$Homozygotes)
cor.test(merged$DON.BLUE_mv, merged$Homozygotes)
cor.test(merged$height.BLUE_raw, merged$Homozygotes)

cor.test(ns_merged$yld.BLUE_mv, ns_merged$Homozygotes)
cor.test(ns_merged$DON.BLUE_mv, ns_merged$Homozygotes)
cor.test(ns_merged$height.BLUE_raw, ns_merged$Homozygotes)

cor.test(s_merged$yld.BLUE_mv, s_merged$Homozygotes)
cor.test(s_merged$DON.BLUE_mv, s_merged$Homozygotes)
cor.test(s_merged$height.BLUE_raw, s_merged$Homozygotes)

# We will separate the cycles for visualization purposes
c0 <- merged[merged$cycle == "C0",]
c1 <- merged[merged$cycle == "C1",]
c2 <- merged[merged$cycle == "C2",]
c3 <- merged[merged$cycle == "C3",]


ns_c0 <- ns_merged[ns_merged$cycle == "C0",]
ns_c1 <- ns_merged[ns_merged$cycle == "C1",]
ns_c2 <- ns_merged[ns_merged$cycle == "C2",]
ns_c3 <- ns_merged[ns_merged$cycle == "C3",]

s_c0 <- s_merged[s_merged$cycle == "C0",]
s_c1 <- s_merged[s_merged$cycle == "C1",]
s_c2 <- s_merged[s_merged$cycle == "C2",]
s_c3 <- s_merged[s_merged$cycle == "C3",]

# Plot the distributions
pdf(file="Deleterious_Phenotype_Correlations.pdf", height=8, width=8)
par(mfrow=c(2, 2), mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
plot(
    c0$yld.BLUE_mv ~ c0$Homozygotes,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Number of Homozygous Deleterious SNPs",
    ylab="Yield BLUE (kg/ha)",
    main="\nYield",
    xlim=c(0, 300),
    ylim=c(4750, 7800))
points(c1$yld.BLUE_mv ~ c1$Homozygotes, pch=16, col="blue", cex=0.6)
points(c2$yld.BLUE_mv ~ c2$Homozygotes, pch=17, col="green", cex=0.6)
points(c3$yld.BLUE_mv ~ c3$Homozygotes, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))

plot(
    c0$DON.BLUE_mv ~ c0$Homozygotes,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Number of Homozygous Deleterious SNPs",
    ylab="DON BLUE (ppm)",
    main="\nDON Concentration",
    xlim=c(0, 300),
    ylim=c(0, 25))
points(c1$DON.BLUE_mv ~ c1$Homozygotes, pch=16, col="blue", cex=0.6)
points(c2$DON.BLUE_mv ~ c2$Homozygotes, pch=17, col="green", cex=0.6)
points(c3$DON.BLUE_mv ~ c3$Homozygotes, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))

plot(
    c0$height.BLUE_raw ~ c0$Homozygotes,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Number of Homozygous Deleterious SNPs",
    ylab="Height BLUE (cm)",
    main="\nHeight",
    xlim=c(0, 300),
    ylim=c(75, 105))
points(c1$height.BLUE_raw ~ c1$Homozygotes, pch=16, col="blue", cex=0.6)
points(c2$height.BLUE_raw ~ c2$Homozygotes, pch=17, col="green", cex=0.6)
points(c3$height.BLUE_raw ~ c3$Homozygotes, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))
dev.off()


pdf(file="Nonsynonymous_Phenotype_Correlations.pdf", height=8, width=8)
par(mfrow=c(2, 2), mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
plot(
    ns_c0$yld.BLUE_mv ~ ns_c0$Homozygotes,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Number of Homozygous Nonsynonymous SNPs",
    ylab="Yield BLUE (kg/ha)",
    main="\nYield",
    xlim=c(3500, 13500),
    ylim=c(4750, 7800))
points(ns_c1$yld.BLUE_mv ~ ns_c1$Homozygotes, pch=16, col="blue", cex=0.6)
points(ns_c2$yld.BLUE_mv ~ ns_c2$Homozygotes, pch=17, col="green", cex=0.6)
points(ns_c3$yld.BLUE_mv ~ ns_c3$Homozygotes, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))

plot(
    ns_c0$DON.BLUE_mv ~ ns_c0$Homozygotes,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Number of Homozygous Nonsynonymous SNPs",
    ylab="DON BLUE (ppm)",
    main="\nDON Concentration",
    xlim=c(3500, 13500),
    ylim=c(0, 25))
points(ns_c1$DON.BLUE_mv ~ ns_c1$Homozygotes, pch=16, col="blue", cex=0.6)
points(ns_c2$DON.BLUE_mv ~ ns_c2$Homozygotes, pch=17, col="green", cex=0.6)
points(ns_c3$DON.BLUE_mv ~ ns_c3$Homozygotes, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))

plot(
    ns_c0$height.BLUE_raw ~ ns_c0$Homozygotes,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Number of Homozygous Nonsynonymous SNPs",
    ylab="Height BLUE (cm)",
    main="\nHeight",
    xlim=c(3500, 13500),
    ylim=c(75, 105))
points(ns_c1$height.BLUE_raw ~ ns_c1$Homozygotes, pch=16, col="blue", cex=0.6)
points(ns_c2$height.BLUE_raw ~ ns_c2$Homozygotes, pch=17, col="green", cex=0.6)
points(ns_c3$height.BLUE_raw ~ ns_c3$Homozygotes, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))
dev.off()

pdf(file="Synynonymous_Phenotype_Correlations.pdf", height=8, width=8)
par(mfrow=c(2, 2), mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
plot(
    s_c0$yld.BLUE_mv ~ s_c0$Homozygotes,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Number of Homozygous Synonymous SNPs",
    ylab="Yield BLUE (kg/ha)",
    main="\nYield",
    xlim=c(3500, 18000),
    ylim=c(4750, 7800))
points(s_c1$yld.BLUE_mv ~ s_c1$Homozygotes, pch=16, col="blue", cex=0.6)
points(s_c2$yld.BLUE_mv ~ s_c2$Homozygotes, pch=17, col="green", cex=0.6)
points(s_c3$yld.BLUE_mv ~ s_c3$Homozygotes, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))

plot(
    s_c0$DON.BLUE_mv ~ s_c0$Homozygotes,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Number of Homozygous Synonymous SNPs",
    ylab="DON BLUE (ppm)",
    main="\nDON Concentration",
    xlim=c(3500, 18000),
    ylim=c(0, 25))
points(s_c1$DON.BLUE_mv ~ s_c1$Homozygotes, pch=16, col="blue", cex=0.6)
points(s_c2$DON.BLUE_mv ~ s_c2$Homozygotes, pch=17, col="green", cex=0.6)
points(s_c3$DON.BLUE_mv ~ s_c3$Homozygotes, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))

plot(
    s_c0$height.BLUE_raw ~ s_c0$Homozygotes,
    pch=15,
    cex=0.6,
    col="black",
    xlab="Number of Homozygous Synonymous SNPs",
    ylab="Height BLUE (cm)",
    main="\nHeight",
    xlim=c(3500, 18000),
    ylim=c(75, 105))
points(s_c1$height.BLUE_raw ~ s_c1$Homozygotes, pch=16, col="blue", cex=0.6)
points(s_c2$height.BLUE_raw ~ s_c2$Homozygotes, pch=17, col="green", cex=0.6)
points(s_c3$height.BLUE_raw ~ s_c3$Homozygotes, pch=19, col="red", cex=0.6)
legend("topright", c("C0", "C1", "C2", "C3"), col=c("black", "blue", "green", "red"), pch=c(15, 16, 17, 19))
dev.off()
