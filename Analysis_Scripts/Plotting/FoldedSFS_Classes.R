# Generate a folded site frequency spectrum for the classes of SNPs, using the
# imputed genotypes.

nc_col <- '#2c7bb6'
syn_col <- '#abd9e9'
ns_col <- '#fdae61'
del_col <- '#d7191c'

# Read the data files
c0.freq <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Genotype_Freqs/C0_MAF.frq.gz", header=TRUE)
c1.freq <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Genotype_Freqs/C1_MAF.frq.gz", header=TRUE)
c2.freq <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Genotype_Freqs/C2_MAF.frq.gz", header=TRUE)
c3.freq <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Genotype_Freqs/C3_MAF.frq.gz", header=TRUE)

nonc <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Noncoding.names", header=FALSE)$V1)
syn <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Synonymous.names", header=FALSE)$V1)
nonsyn <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Nonsynonymous.names", header=FALSE)$V1)
del <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Deleterious.names", header=FALSE)$V1)

# slice the frequency table up by the classes
c0.nonc.freq <- c0.freq$MAF[c0.freq$SNP %in% nonc]
c0.syn.freq <- c0.freq$MAF[c0.freq$SNP %in% syn]
c0.nonsyn.freq <- c0.freq$MAF[c0.freq$SNP %in% nonsyn]
c0.del.freq <- c0.freq$MAF[c0.freq$SNP %in% del]

c1.nonc.freq <- c1.freq$MAF[c1.freq$SNP %in% nonc]
c1.syn.freq <- c1.freq$MAF[c1.freq$SNP %in% syn]
c1.nonsyn.freq <- c1.freq$MAF[c1.freq$SNP %in% nonsyn]
c1.del.freq <- c1.freq$MAF[c1.freq$SNP %in% del]

c2.nonc.freq <- c2.freq$MAF[c2.freq$SNP %in% nonc]
c2.syn.freq <- c2.freq$MAF[c2.freq$SNP %in% syn]
c2.nonsyn.freq <- c2.freq$MAF[c2.freq$SNP %in% nonsyn]
c2.del.freq <- c2.freq$MAF[c2.freq$SNP %in% del]

c3.nonc.freq <- c3.freq$MAF[c3.freq$SNP %in% nonc]
c3.syn.freq <- c3.freq$MAF[c3.freq$SNP %in% syn]
c3.nonsyn.freq <- c3.freq$MAF[c3.freq$SNP %in% nonsyn]
c3.del.freq <- c3.freq$MAF[c3.freq$SNP %in% del]

# Bin them up in 2.5% frequency bins
breaks <- seq(0, 0.5, by=0.025)
c0.nonc.sfs <- table(cut(c0.nonc.freq, breaks=breaks, include.lowest=TRUE))/length(c0.nonc.freq)
c0.syn.sfs <- table(cut(c0.syn.freq, breaks=breaks, include.lowest=TRUE))/length(c0.syn.freq)
c0.nonsyn.sfs <- table(cut(c0.nonsyn.freq, breaks=breaks, include.lowest=TRUE))/length(c0.nonsyn.freq)
c0.del.sfs <- table(cut(c0.del.freq, breaks=breaks, include.lowest=TRUE))/length(c0.del.freq)

c1.nonc.sfs <- table(cut(c1.nonc.freq, breaks=breaks, include.lowest=TRUE))/length(c1.nonc.freq)
c1.syn.sfs <- table(cut(c1.syn.freq, breaks=breaks, include.lowest=TRUE))/length(c1.syn.freq)
c1.nonsyn.sfs <- table(cut(c1.nonsyn.freq, breaks=breaks, include.lowest=TRUE))/length(c1.nonsyn.freq)
c1.del.sfs <- table(cut(c1.del.freq, breaks=breaks, include.lowest=TRUE))/length(c1.del.freq)

c2.nonc.sfs <- table(cut(c2.nonc.freq, breaks=breaks, include.lowest=TRUE))/length(c2.nonc.freq)
c2.syn.sfs <- table(cut(c2.syn.freq, breaks=breaks, include.lowest=TRUE))/length(c2.syn.freq)
c2.nonsyn.sfs <- table(cut(c2.nonsyn.freq, breaks=breaks, include.lowest=TRUE))/length(c2.nonsyn.freq)
c2.del.sfs <- table(cut(c2.del.freq, breaks=breaks, include.lowest=TRUE))/length(c2.del.freq)

c3.nonc.sfs <- table(cut(c3.nonc.freq, breaks=breaks, include.lowest=TRUE))/length(c3.nonc.freq)
c3.syn.sfs <- table(cut(c3.syn.freq, breaks=breaks, include.lowest=TRUE))/length(c3.syn.freq)
c3.nonsyn.sfs <- table(cut(c3.nonsyn.freq, breaks=breaks, include.lowest=TRUE))/length(c3.nonsyn.freq)
c3.del.sfs <- table(cut(c3.del.freq, breaks=breaks, include.lowest=TRUE))/length(c3.del.freq)

# Make a matrix of it
c0.toplot <- cbind(c0.nonc.sfs, c0.syn.sfs, c0.nonsyn.sfs, c0.del.sfs)
c1.toplot <- cbind(c1.nonc.sfs, c1.syn.sfs, c1.nonsyn.sfs, c1.del.sfs)
c2.toplot <- cbind(c2.nonc.sfs, c2.syn.sfs, c2.nonsyn.sfs, c2.del.sfs)
c3.toplot <- cbind(c3.nonc.sfs, c3.syn.sfs, c3.nonsyn.sfs, c3.del.sfs)
pdf(file="Imputed_Folded_SFS.pdf", width=8.5, height=11)
par(mfrow=c(4, 1))
barplot(t(c0.toplot),
    beside=TRUE,
    col=c(nc_col, syn_col, ns_col, del_col),
    xlab="Minor Allele Frequency",
    ylab="Proportion",
    main="Folded SFS in Parents")
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"), fill=c(nc_col, syn_col, ns_col, del_col))
barplot(t(c1.toplot),
    beside=TRUE,
    col=c(nc_col, syn_col, ns_col, del_col),
    xlab="Minor Allele Frequency",
    ylab="Proportion",
    main="Folded SFS in Cycle 1")
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"), fill=c(nc_col, syn_col, ns_col, del_col))
barplot(t(c2.toplot),
    beside=TRUE,
    col=c(nc_col, syn_col, ns_col, del_col),
    xlab="Minor Allele Frequency",
    ylab="Proportion",
    main="Folded SFS in Cycle 2")
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"), fill=c(nc_col, syn_col, ns_col, del_col))
barplot(t(c3.toplot),
    beside=TRUE,
    col=c(nc_col, syn_col, ns_col, del_col),
    xlab="Minor Allele Frequency",
    ylab="Proportion",
    main="Folded SFS in Cycle 3")
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"), fill=c(nc_col, syn_col, ns_col, del_col))
dev.off()
