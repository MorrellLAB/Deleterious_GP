# Plot an unfolded SFS for the functional classes in the parental data.

freqs <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Ancestral_State/GP_Ancestral.txt.gz", header=TRUE)
nonc <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Noncoding.txt", header=FALSE)$V1)
syn <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Synonymous.txt", header=FALSE)$V1)
nonsyn <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Nonsynonymous.txt", header=FALSE)$V1)
del <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Deleterious.txt", header=FALSE)$V1)

# Separate the classes
nc_freqs <- freqs[freqs$SNPID %in% nonc,]
s_freqs <- freqs[freqs$SNPID %in% syn,]
ns_freqs <- freqs[(freqs$SNPID %in% nonsyn) & !(freqs$SNPID %in% del),]
del_freqs <- freqs[freqs$SNPID %in% del,]

# Chop them up by frequency
bins <- seq(0, 1, by=0.1)

nc_sfs <- table(cut(nc_freqs$DAF, breaks=bins, include.lowest=TRUE))
s_sfs <- table(cut(s_freqs$DAF, breaks=bins, include.lowest=TRUE))
ns_sfs <- table(cut(ns_freqs$DAF, breaks=bins, include.lowest=TRUE))
del_sfs <- table(cut(del_freqs$DAF, breaks=bins, include.lowest=TRUE))

nc_sfs <- nc_sfs/sum(nc_sfs)
s_sfs <- s_sfs/sum(s_sfs)
ns_sfs <- ns_sfs/sum(ns_sfs)
del_sfs <- del_sfs/sum(del_sfs)

toplot <- cbind(nc_sfs, s_sfs, ns_sfs, del_sfs)

pdf(file="Parental_Derived_SFS_By_Class.pdf", height=4, width=8)
at <- barplot(t(toplot),
    beside=TRUE,
    col=c("blacK", "blue", "green", "red"),
    xlab="Derived Allele Frequency",
    ylab="Proportion",
    main="Unfolded SFS in Parents",
    axes=FALSE,
    ylim=c(0, 0.7))
axis(side=2)
axis(side=1, at=apply(at, 2, mean), labels=names(nc_sfs))
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"), fill=c("black", "blue", "green", "red"))
dev.off()
