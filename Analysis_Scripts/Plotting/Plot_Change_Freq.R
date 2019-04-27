# Plot the difference in derived allele frequency for various functional
# classes of SNPs across the genome

library(beeswarm)

CHROMS <- c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H")
# nc_col <- "#2c7bb6"
nc_col <- rgb(44/255, 123/255, 182/255, 0.15)
# syn_col <- "#abd9e9"
syn_col <- rgb(171/255, 217/255, 233/255, 0.15)
# ns_col <- "#fdae61"
ns_col <- rgb(253/255, 174/255, 97/255, 0.15)
# del_col <- "#d7191c"
del_col <- rgb(215/255, 25/255, 28/255, 0.4)

# Read the data files
freqs <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Genotype_Freqs/DAF_By_Cycle.txt", header=TRUE)
nonc <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Noncoding.names.gz", header=FALSE)$V1)
syn <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Synonymous.names.gz", header=FALSE)$V1)
nonsyn <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Nonsynonymous.names.gz", header=FALSE)$V1)
del <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Deleterious.names.gz", header=FALSE)$V1)
pos <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Imputation/AlphaPeel/GP_AP.bim", header=FALSE)

# Set the names of the columns of the positions
names(pos) <- c("CHR", "SNP_ID", "GPOS", "BP", "A1", "A2")

# Merge the positions into the frequencies
posfreqs <- merge(freqs, pos, by="SNP_ID")
posfreqs <- posfreqs[order(posfreqs$CHR, posfreqs$BP),]
posfreqs$BP <- posfreqs$BP/1000000

# Calculate a new column for the change in frequency
posfreqs$ChFrq <- (posfreqs$C3_DAF - posfreqs$C1_DAF)/posfreqs$C1_DAF
# Add the centromeres
# centromeres <- data.frame(
#     Chr=c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H"),
#     Start=c(NA, 136904248, 115794456, 58804724, 60207685, 180782470, 202235859),
#     End=c(NA, 528069469, 438235281, 467767377, 368181264, 383275592, 449173785))
pericentromeres <- data.frame(
    Chr=c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H"),
    Start=c(70000000, 140000000, 140000000, 50000000, 50000000, 80000000, 175000000),
    End=c(240000000, 425000000, 315000000, 325000000, 260000000, 325000000, 375000000))
# The actual centromeres are from Mascher et al. 2017
centromeres <- data.frame(
    Chr=c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H"),
    Start=c(213600001, 333600001, 288000001, 286400001, 217600001, 257600001, 339200001),
    End=c(221600000, 341600000, 299200000, 296000000, 225600000, 265600000, 345600000))

# Separate the functional classes
df_nc <- posfreqs[posfreqs$SNP_ID %in% nonc,]
df_syn <- posfreqs[posfreqs$SNP_ID %in% syn,]
df_nonsyn <- posfreqs[(posfreqs$SNP_ID %in% nonsyn) & !(posfreqs$SNP_ID %in% del),]
df_del <- posfreqs[posfreqs$SNP_ID %in% del,]

df_nc <- df_nc[!is.na(df_nc$ChFrq),]
df_syn <- df_syn[!is.na(df_syn$ChFrq),]
df_nonsyn <- df_nonsyn[!is.na(df_nonsyn$ChFrq),]
df_del <- df_del[!is.na(df_del$ChFrq),]

# Make counts of increasing and decreasing
nc.inc <- sum(df_nc$ChFrq > 0) / nrow(df_nc)
nc.dec <- sum(df_nc$ChFrq < 0) / nrow(df_nc)
syn.inc <- sum(df_syn$ChFrq > 0) / nrow(df_syn)
syn.dec <- sum(df_syn$ChFrq < 0) / nrow(df_syn)
nonsyn.inc <- sum(df_nonsyn$ChFrq > 0) / nrow(df_nonsyn)
nonsyn.dec <- sum(df_nonsyn$ChFrq < 0) / nrow(df_nonsyn)
del.inc <- sum(df_del$ChFrq > 0) / nrow(df_del)
del.dec <- sum(df_del$ChFrq < 0) / nrow(df_del)

# Mean change in DAF

inc_dec <- data.frame(
    Class=c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"),
    MeanDeltaDAF=c(mean(df_nc$ChFrq), mean(df_syn$ChFrq), mean(df_nonsyn$ChFrq), mean(df_del$ChFrq)),
    MedianDeltaDAF=c(median(df_nc$ChFrq), median(df_syn$ChFrq), median(df_nonsyn$ChFrq), median(df_del$ChFrq)),
    Increasing=c(nc.inc, syn.inc, nonsyn.inc, del.inc),
    Decreasing=c(nc.dec, syn.dec, nonsyn.dec, del.dec))
print(inc_dec)

# # Make the 4H region expansion plot
# pdf(file="4H_ROI_DeltaDAF.pdf", height=3, width=6)
# nc_roi <- df_nc[df_nc$CHR == "chr4H" & df_nc$BP > 18 & df_nc$BP < 36,]
# syn_roi <- df_syn[df_syn$CHR == "chr4H" & df_syn$BP > 18 & df_syn$BP < 36,]
# nonsyn_roi <- df_nonsyn[df_nonsyn$CHR == "chr4H" & df_nonsyn$BP > 18 & df_nonsyn$BP < 36,]
# del_roi <- df_del[df_del$CHR == "chr4H" & df_del$BP > 18 & df_del$BP < 36,]
# plot(nc_roi$ChFrq ~ nc_roi$BP, type="l", lwd=1, col=nc_col, ylim=c(-1, 5), xlim=c(18, 36), xlab="Position (Mb)", ylab="Fold Change in DAF", main="4H Yield and DON Region", axes=F)
# axis(side=2)
# axis(side=1, at=seq(18, 36, by=2), labels=seq(18, 36, by=2))
# lines(syn_roi$ChFrq ~ syn_roi$BP, lwd=1, col=syn_col)
# lines(nonsyn_roi$ChFrq ~ nonsyn_roi$BP, lwd=1, col=ns_col)
# lines(del_roi$ChFrq ~ del_roi$BP, lwd=2, col=del_col)
# abline(h=0, lwd=2, col=nc_col, lty=3)
# dev.off()

# Plot it
png(file="Delta_DAF.png", res=300, height=10.5*300, width=8*300)
par(mfrow=c(7, 1), mar=c(4, 4, 1, 2))
for(chrom in CHROMS) {
    nc_toplot <- df_nc[df_nc$CHR == chrom,]
    syn_toplot <- df_syn[df_syn$CHR == chrom,]
    nonsyn_toplot <- df_nonsyn[df_nonsyn$CHR == chrom,]
    del_toplot <- df_del[df_del$CHR == chrom,]
    cent <- centromeres[centromeres$Chr == chrom,]
    peri <- pericentromeres[pericentromeres$Chr == chrom,]
    plot(nc_toplot$ChFrq ~ nc_toplot$BP, type="n", ylim=c(-1, 5), xlim=c(0, 800), xlab="Position (Mb)", ylab="Fold Change in DAF", main=chrom, axes=F)
    rect(peri$Start/1000000, -2, peri$End/1000000, 5, density=NA, col=rgb(0, 0, 0, alpha=0.2))
    rect(cent$Start/1000000, -2, cent$End/1000000, 5, density=NA, col=rgb(0, 0, 0, alpha=0.3))
    points(nc_toplot$ChFrq ~ nc_toplot$BP, type="p", cex=0.5, col=nc_col)
    axis(side=2)
    axis(side=1, at=seq(0, 800, by=50), labels=seq(0, 800, by=50))
    points(syn_toplot$ChFrq ~ syn_toplot$BP, pch=19, cex=0.5, col=syn_col)
    points(nonsyn_toplot$ChFrq ~ nonsyn_toplot$BP, pch=19, cex=0.5, col=ns_col)
    points(del_toplot$ChFrq ~ del_toplot$BP, pch=19, col=del_col, cex=0.75)
    abline(h=0, lwd=2, col="black", lty=3)
}
dev.off()

# pdf(file="Delta_DAF_Boxplot.pdf", 6, 6)
png(file="Delta_DAF_Beeswarm.png", res=300, height=1800, width=1800)
toplot <- data.frame(
    DDAF=c(df_nc$ChFrq, df_syn$ChFrq, df_nonsyn$ChFrq, df_del$ChFrq),
    Class=c(rep("Noncoding", nrow(df_nc)), rep("Synonymous", nrow(df_syn)), rep("Nonsynonymous", nrow(df_nonsyn)), rep("Deleterious", nrow(df_del))))
toplot$Class <- factor(toplot$Class, levels=c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"), ordered=TRUE)
beeswarm(toplot$DDAF ~ toplot$Class,
    xlab="Functional Class",
    ylab="Fold Change in DAF",
    main="Fold Change in DAF Over Three Cycles of Selection",
    col=c(nc_col, syn_col, ns_col, del_col),
    pch=19,
    cex=0.1,
    method="hex",
    corral="random",
    cex.axis=0.8,
    ylim=c(-1, 5))
bxplot(
    toplot$DDAF ~ toplot$Class,
    add=TRUE)
abline(h=0, lwd=1, col="black", lty=3)
dev.off()

# Fit an ANOVA
daf.aov <- aov(DDAF ~ Class, data=toplot)
summary(daf.aov)
TukeyHSD(daf.aov)
