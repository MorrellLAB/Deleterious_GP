# Plot the unfolded site frequency spectra for each cycle, separated by the
# functional class.

nc_col <- '#2c7bb6'
syn_col <- '#abd9e9'
ns_col <- '#fdae61'
del_col <- '#d7191c'

freqs <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Genotype_Freqs/DAF_By_Cycle.txt.gz", header=TRUE)
nonc <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Noncoding.names.gz", header=FALSE)$V1)
syn <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Synonymous.names.gz", header=FALSE)$V1)
nonsyn <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Nonsynonymous.names.gz", header=FALSE)$V1)
del <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Deleterious.names.gz", header=FALSE)$V1)

# Separate the functional classes
nc_freqs <- freqs[freqs$SNP_ID %in% nonc,]
s_freqs <- freqs[freqs$SNP_ID %in% syn,]
ns_freqs <- freqs[(freqs$SNP_ID %in% nonsyn) & !(freqs$SNP_ID %in% del),]
del_freqs <- freqs[freqs$SNP_ID %in% del,]

# Chop them up by frequency
bins <- seq(0, 1, by=0.05)

c1_nc_sfs <- table(cut(nc_freqs$C1_DAF, breaks=bins, include.lowest=TRUE))
c2_nc_sfs <- table(cut(nc_freqs$C2_DAF, breaks=bins, include.lowest=TRUE))
c3_nc_sfs <- table(cut(nc_freqs$C3_DAF, breaks=bins, include.lowest=TRUE))

c1_s_sfs <- table(cut(s_freqs$C1_DAF, breaks=bins, include.lowest=TRUE))
c2_s_sfs <- table(cut(s_freqs$C2_DAF, breaks=bins, include.lowest=TRUE))
c3_s_sfs <- table(cut(s_freqs$C3_DAF, breaks=bins, include.lowest=TRUE))

c1_ns_sfs <- table(cut(ns_freqs$C1_DAF, breaks=bins, include.lowest=TRUE))
c2_ns_sfs <- table(cut(ns_freqs$C2_DAF, breaks=bins, include.lowest=TRUE))
c3_ns_sfs <- table(cut(ns_freqs$C3_DAF, breaks=bins, include.lowest=TRUE))

c1_del_sfs <- table(cut(del_freqs$C1_DAF, breaks=bins, include.lowest=TRUE))
c2_del_sfs <- table(cut(del_freqs$C2_DAF, breaks=bins, include.lowest=TRUE))
c3_del_sfs <- table(cut(del_freqs$C3_DAF, breaks=bins, include.lowest=TRUE))

# Turn them into proportions
c1_nc_sfs <- c1_nc_sfs/sum(c1_nc_sfs)
c2_nc_sfs <- c2_nc_sfs/sum(c2_nc_sfs)
c3_nc_sfs <- c3_nc_sfs/sum(c3_nc_sfs)

c1_s_sfs <- c1_s_sfs/sum(c1_s_sfs)
c2_s_sfs <- c2_s_sfs/sum(c2_s_sfs)
c3_s_sfs <- c3_s_sfs/sum(c3_s_sfs)

c1_ns_sfs <- c1_ns_sfs/sum(c1_ns_sfs)
c2_ns_sfs <- c2_ns_sfs/sum(c2_ns_sfs)
c3_ns_sfs <- c3_ns_sfs/sum(c3_ns_sfs)

c1_del_sfs <- c1_del_sfs/sum(c1_del_sfs)
c2_del_sfs <- c2_del_sfs/sum(c2_del_sfs)
c3_del_sfs <- c3_del_sfs/sum(c3_del_sfs)

# Make matrices for plotting
c1 <- cbind(c1_nc_sfs, c1_s_sfs, c1_ns_sfs, c1_del_sfs)
c2 <- cbind(c2_nc_sfs, c2_s_sfs, c2_ns_sfs, c2_del_sfs)
c3 <- cbind(c3_nc_sfs, c3_s_sfs, c3_ns_sfs, c3_del_sfs)

# And plot them
pdf(file="Derived_SFS_By_Class.pdf", height=10.5, width=8)
par(mfrow=c(3, 1), mar=c(4, 4, 1, 2))
at <- barplot(t(c1),
    beside=TRUE,
    col=c(nc_col, syn_col, ns_col, del_col),
    xlab="Derived Allele Frequency",
    ylab="Proportion",
    main="Cycle 1",
    axes=FALSE,
    ylim=c(0, 0.7))
axis(side=2)
axis(side=1, at=apply(at, 2, mean), labels=names(c1_nc_sfs))
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"), fill=c(nc_col, syn_col, ns_col, del_col))
at <- barplot(t(c2),
    beside=TRUE,
    col=c(nc_col, syn_col, ns_col, del_col),
    xlab="Derived Allele Frequency",
    ylab="Proportion",
    main="Cycle 2",
    axes=FALSE,
    ylim=c(0, 0.7))
axis(side=2)
axis(side=1, at=apply(at, 2, mean), labels=names(c1_nc_sfs))
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"), fill=c(nc_col, syn_col, ns_col, del_col))
at <- barplot(t(c3),
    beside=TRUE,
    col=c(nc_col, syn_col, ns_col, del_col),
    xlab="Derived Allele Frequency",
    ylab="Proportion",
    main="Cycle 3",
    axes=FALSE,
    ylim=c(0, 0.7))
axis(side=2)
axis(side=1, at=apply(at, 2, mean), labels=names(c1_nc_sfs))
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"), fill=c(nc_col, syn_col, ns_col, del_col))
dev.off()
