#   Script to plot the proportions of homozygous major/minor and het calls for
#   each minor allele count class

freebayes <- read.table("/Volumes/Data_Disk/Dropbox/GitHub/Deleterious_GP/Data/FreeBayes_GenoCounts.txt", header=T)
gatk <- read.table("/Volumes/Data_Disk/Dropbox/GitHub/Deleterious_GP/Data/GATK_GenoCounts.txt", header=T)

#   Make some plots of these. They are ready to be coerced into matrices and
#   fed (transposed) into barplot()
pdf(file="FreeBayes_Genotype_Classes.pdf", width=11, height=6)
toplot <- as.matrix(freebayes[2:4])
at <- barplot(
    t(toplot),
    beside=TRUE,
    ylim=c(0, 1.1),
    main="FreeBayes Genotype Class Summary on 7H",
    xlab="Minor Allele Count",
    ylab="Proportion of Genotypes",
    col=c("darkred", "purple", "darkblue")
    )
axis(
    side=1,
    at=apply(at, 2, mean),
    labels=freebayes$MAC
    )
legend(
    "topright",
    c("Hom. Min.", "Het.", "Hom. Maj."),
    fill=c("darkred", "purple", "darkblue")
    )
dev.off()

#   And do the same for GATK
pdf(file="GATK_Genotype_Classes.pdf", width=11, height=6)
toplot <- as.matrix(gatk[2:4])
at <- barplot(
    t(toplot),
    beside=TRUE,
    ylim=c(0, 1.1),
    main="GATK Genotype Class Summary on 7H",
    xlab="Minor Allele Count",
    ylab="Proportion of Genotypes",
    col=c("darkred", "purple", "darkblue")
    )
axis(
    side=1,
    at=apply(at, 2, mean),
    labels=gatk$MAC
    )
legend(
    "topright",
    c("Hom. Min.", "Het.", "Hom. Maj."),
    fill=c("darkred", "purple", "darkblue")
    )
dev.off()
