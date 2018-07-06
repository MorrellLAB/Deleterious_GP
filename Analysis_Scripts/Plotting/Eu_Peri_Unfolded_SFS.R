# Plot the unfolded site frequency spectra for each cycle, separated by the
# functional class and euchromatic and pericentromere.

freqs <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Genotype_Freqs/DAF_By_Cycle.txt.gz", header=TRUE)
peri <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Pericentromeric.txt", header=FALSE)$V1)
eu <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Euchromatic.txt", header=FALSE)$V1)
nonc <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Noncoding.txt", header=FALSE)$V1)
syn <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Synonymous.txt", header=FALSE)$V1)
nonsyn <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Nonsynonymous.txt", header=FALSE)$V1)
del <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Deleterious.txt", header=FALSE)$V1)

# Separate the functional classes
peri.nc_freqs <- freqs[(freqs$SNP_ID %in% nonc) & (freqs$SNP_ID %in% peri),]
peri.s_freqs <- freqs[(freqs$SNP_ID %in% syn) & (freqs$SNP_ID %in% peri),]
peri.ns_freqs <- freqs[(freqs$SNP_ID %in% nonsyn) & !(freqs$SNP_ID %in% del) & (freqs$SNP_ID %in% peri),]
peri.del_freqs <- freqs[(freqs$SNP_ID %in% del) & (freqs$SNP_ID %in% peri),]

eu.nc_freqs <- freqs[(freqs$SNP_ID %in% nonc) & (freqs$SNP_ID %in% eu),]
eu.s_freqs <- freqs[(freqs$SNP_ID %in% syn) & (freqs$SNP_ID %in% eu),]
eu.ns_freqs <- freqs[(freqs$SNP_ID %in% nonsyn) & !(freqs$SNP_ID %in% del) & (freqs$SNP_ID %in% eu),]
eu.del_freqs <- freqs[(freqs$SNP_ID %in% del) & (freqs$SNP_ID %in% eu),]

bins <- seq(0, 1, by=0.05)

peri.c1_nc_sfs <- table(cut(peri.nc_freqs$C1_DAF, breaks=bins, include.lowest=TRUE))
peri.c2_nc_sfs <- table(cut(peri.nc_freqs$C2_DAF, breaks=bins, include.lowest=TRUE))
peri.c3_nc_sfs <- table(cut(peri.nc_freqs$C3_DAF, breaks=bins, include.lowest=TRUE))

peri.c1_s_sfs <- table(cut(peri.s_freqs$C1_DAF, breaks=bins, include.lowest=TRUE))
peri.c2_s_sfs <- table(cut(peri.s_freqs$C2_DAF, breaks=bins, include.lowest=TRUE))
peri.c3_s_sfs <- table(cut(peri.s_freqs$C3_DAF, breaks=bins, include.lowest=TRUE))

peri.c1_ns_sfs <- table(cut(peri.ns_freqs$C1_DAF, breaks=bins, include.lowest=TRUE))
peri.c2_ns_sfs <- table(cut(peri.ns_freqs$C2_DAF, breaks=bins, include.lowest=TRUE))
peri.c3_ns_sfs <- table(cut(peri.ns_freqs$C3_DAF, breaks=bins, include.lowest=TRUE))

peri.c1_del_sfs <- table(cut(peri.del_freqs$C1_DAF, breaks=bins, include.lowest=TRUE))
peri.c2_del_sfs <- table(cut(peri.del_freqs$C2_DAF, breaks=bins, include.lowest=TRUE))
peri.c3_del_sfs <- table(cut(peri.del_freqs$C3_DAF, breaks=bins, include.lowest=TRUE))

# Turn them into proportions
peri.c1_nc_sfs <- peri.c1_nc_sfs/sum(peri.c1_nc_sfs)
peri.c2_nc_sfs <- peri.c2_nc_sfs/sum(peri.c2_nc_sfs)
peri.c3_nc_sfs <- peri.c3_nc_sfs/sum(peri.c3_nc_sfs)

peri.c1_s_sfs <- peri.c1_s_sfs/sum(peri.c1_s_sfs)
peri.c2_s_sfs <- peri.c2_s_sfs/sum(peri.c2_s_sfs)
peri.c3_s_sfs <- peri.c3_s_sfs/sum(peri.c3_s_sfs)

peri.c1_ns_sfs <- peri.c1_ns_sfs/sum(peri.c1_ns_sfs)
peri.c2_ns_sfs <- peri.c2_ns_sfs/sum(peri.c2_ns_sfs)
peri.c3_ns_sfs <- peri.c3_ns_sfs/sum(peri.c3_ns_sfs)

peri.c1_del_sfs <- peri.c1_del_sfs/sum(peri.c1_del_sfs)
peri.c2_del_sfs <- peri.c2_del_sfs/sum(peri.c2_del_sfs)
peri.c3_del_sfs <- peri.c3_del_sfs/sum(peri.c3_del_sfs)

# Make matrices for plotting
peri.c1 <- cbind(peri.c1_nc_sfs, peri.c1_s_sfs, peri.c1_ns_sfs, peri.c1_del_sfs)
peri.c2 <- cbind(peri.c2_nc_sfs, peri.c2_s_sfs, peri.c2_ns_sfs, peri.c2_del_sfs)
peri.c3 <- cbind(peri.c3_nc_sfs, peri.c3_s_sfs, peri.c3_ns_sfs, peri.c3_del_sfs)

# And plot them
pdf(file="Peri_Derived_SFS_By_Class.pdf", height=10.5, width=8)
par(mfrow=c(3, 1), mar=c(4, 4, 1, 2))
at <- barplot(t(peri.c1),
    beside=TRUE,
    col=c("blacK", "blue", "green", "red"),
    xlab="Derived Allele Frequency",
    ylab="Proportion",
    main="Cycle 1 Pericentromere",
    axes=FALSE,
    ylim=c(0, 0.7))
axis(side=2)
axis(side=1, at=apply(at, 2, mean), labels=names(peri.c1_nc_sfs))
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"), fill=c("black", "blue", "green", "red"))
at <- barplot(t(peri.c2),
    beside=TRUE,
    col=c("blacK", "blue", "green", "red"),
    xlab="Derived Allele Frequency",
    ylab="Proportion",
    main="Cycle 2 Pericentromere",
    axes=FALSE,
    ylim=c(0, 0.7))
axis(side=2)
axis(side=1, at=apply(at, 2, mean), labels=names(peri.c1_nc_sfs))
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"), fill=c("black", "blue", "green", "red"))
at <- barplot(t(peri.c3),
    beside=TRUE,
    col=c("blacK", "blue", "green", "red"),
    xlab="Derived Allele Frequency",
    ylab="Proportion",
    main="Cycle 3 Pericentromere",
    axes=FALSE,
    ylim=c(0, 0.7))
axis(side=2)
axis(side=1, at=apply(at, 2, mean), labels=names(peri.c1_nc_sfs))
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"), fill=c("black", "blue", "green", "red"))
dev.off()

eu.c1_nc_sfs <- table(cut(eu.nc_freqs$C1_DAF, breaks=bins, include.lowest=TRUE))
eu.c2_nc_sfs <- table(cut(eu.nc_freqs$C2_DAF, breaks=bins, include.lowest=TRUE))
eu.c3_nc_sfs <- table(cut(eu.nc_freqs$C3_DAF, breaks=bins, include.lowest=TRUE))

eu.c1_s_sfs <- table(cut(eu.s_freqs$C1_DAF, breaks=bins, include.lowest=TRUE))
eu.c2_s_sfs <- table(cut(eu.s_freqs$C2_DAF, breaks=bins, include.lowest=TRUE))
eu.c3_s_sfs <- table(cut(eu.s_freqs$C3_DAF, breaks=bins, include.lowest=TRUE))

eu.c1_ns_sfs <- table(cut(eu.ns_freqs$C1_DAF, breaks=bins, include.lowest=TRUE))
eu.c2_ns_sfs <- table(cut(eu.ns_freqs$C2_DAF, breaks=bins, include.lowest=TRUE))
eu.c3_ns_sfs <- table(cut(eu.ns_freqs$C3_DAF, breaks=bins, include.lowest=TRUE))

eu.c1_del_sfs <- table(cut(eu.del_freqs$C1_DAF, breaks=bins, include.lowest=TRUE))
eu.c2_del_sfs <- table(cut(eu.del_freqs$C2_DAF, breaks=bins, include.lowest=TRUE))
eu.c3_del_sfs <- table(cut(eu.del_freqs$C3_DAF, breaks=bins, include.lowest=TRUE))

# Turn them into proportions
eu.c1_nc_sfs <- eu.c1_nc_sfs/sum(eu.c1_nc_sfs)
eu.c2_nc_sfs <- eu.c2_nc_sfs/sum(eu.c2_nc_sfs)
eu.c3_nc_sfs <- eu.c3_nc_sfs/sum(eu.c3_nc_sfs)

eu.c1_s_sfs <- eu.c1_s_sfs/sum(eu.c1_s_sfs)
eu.c2_s_sfs <- eu.c2_s_sfs/sum(eu.c2_s_sfs)
eu.c3_s_sfs <- eu.c3_s_sfs/sum(eu.c3_s_sfs)

eu.c1_ns_sfs <- eu.c1_ns_sfs/sum(eu.c1_ns_sfs)
eu.c2_ns_sfs <- eu.c2_ns_sfs/sum(eu.c2_ns_sfs)
eu.c3_ns_sfs <- eu.c3_ns_sfs/sum(eu.c3_ns_sfs)

eu.c1_del_sfs <- eu.c1_del_sfs/sum(eu.c1_del_sfs)
eu.c2_del_sfs <- eu.c2_del_sfs/sum(eu.c2_del_sfs)
eu.c3_del_sfs <- eu.c3_del_sfs/sum(eu.c3_del_sfs)

# Make matrices for plotting
eu.c1 <- cbind(eu.c1_nc_sfs, eu.c1_s_sfs, eu.c1_ns_sfs, eu.c1_del_sfs)
eu.c2 <- cbind(eu.c2_nc_sfs, eu.c2_s_sfs, eu.c2_ns_sfs, eu.c2_del_sfs)
eu.c3 <- cbind(eu.c3_nc_sfs, eu.c3_s_sfs, eu.c3_ns_sfs, eu.c3_del_sfs)

# And plot them
pdf(file="Eu_Derived_SFS_By_Class.pdf", height=10.5, width=8)
par(mfrow=c(3, 1), mar=c(4, 4, 1, 2))
at <- barplot(t(eu.c1),
    beside=TRUE,
    col=c("blacK", "blue", "green", "red"),
    xlab="Derived Allele Frequency",
    ylab="Proportion",
    main="Cycle 1 Euchromatin",
    axes=FALSE,
    ylim=c(0, 0.7))
axis(side=2)
axis(side=1, at=apply(at, 2, mean), labels=names(eu.c1_nc_sfs))
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"), fill=c("black", "blue", "green", "red"))
at <- barplot(t(eu.c2),
    beside=TRUE,
    col=c("blacK", "blue", "green", "red"),
    xlab="Derived Allele Frequency",
    ylab="Proportion",
    main="Cycle 2 Euchromatin",
    axes=FALSE,
    ylim=c(0, 0.7))
axis(side=2)
axis(side=1, at=apply(at, 2, mean), labels=names(eu.c1_nc_sfs))
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"), fill=c("black", "blue", "green", "red"))
at <- barplot(t(eu.c3),
    beside=TRUE,
    col=c("blacK", "blue", "green", "red"),
    xlab="Derived Allele Frequency",
    ylab="Proportion",
    main="Cycle 3 Euchromatin",
    axes=FALSE,
    ylim=c(0, 0.7))
axis(side=2)
axis(side=1, at=apply(at, 2, mean), labels=names(eu.c1_nc_sfs))
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"), fill=c("black", "blue", "green", "red"))
dev.off()
