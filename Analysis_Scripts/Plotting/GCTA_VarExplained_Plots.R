#   Make a side-by-side boxplot to show the distribution of variance explained
#   by random subsets of 101 markers across various partitions of marker
#   classes.

#   Set working directory
#setwd("/Users/tomkono/DataDisk/Dropbox/Projects/DM_GenomicPrediction/Data/Progeny_Genotypes/GCTA")
#setwd("/home/tom/DataDisk/Dropbox/Projects/DM_GenomicPrediction/Data/Progeny_Genotypes/GCTA")
setwd("/Volumes/Data_Disk/VBox_Shared/PropPheno")
#   Load ggplot2
library(ggplot2)

#   Read the variance data. There is no header. The first column is the
#   proportion of variance explained, and the second column is the SE of the
#   proportion.
null_snps <- read.table("Null/VG_VP.txt", header=F)
all_snps <- read.table("AllSNPs/VG_VP.txt", header=F)
noncoding <- read.table("Noncoding/VG_VP.txt", header=F)
coding <- read.table("Coding/VG_VP.txt", header=F)
nonsyn <- read.table("Nonsyn/VG_VP.txt", header=F)
del <- read.table("Deleterious/VG_VP.txt", header=F)

#   Then, put them into a dataframe for plotting
plot_data <- data.frame(
    PropVar=c(null_snps$V1, all_snps$V1, noncoding$V1, coding$V1, nonsyn$V1, del$V1),
    Partition=factor(c(
        rep("Null", length(null_snps$V1)),
        rep("All", length(all_snps$V1)),
        rep("Noncoding", length(noncoding$V1)),
        rep("Coding", length(coding$V1)),
        rep("Nonsynonymous", length(nonsyn$V1)),
        rep("Deleterious", length(del$V1))), levels=c("Null", "All", "Noncoding", "Coding", "Nonsynonymous", "Deleterious"))
    )

#   Make a plot
pdf(
    file="Proportion_Var_Explained.pdf",
    width=6,
    height=6)
plt <- ggplot(
    plot_data,
    aes(x=factor(Partition), y=PropVar, fill=factor(Partition)))
plt + geom_boxplot() +
    scale_fill_manual(
        name="Partition",
        breaks=c("Null", "All", "Noncoding", "Coding", "Nonsynonymous", "Deleteirous"),
        values=c("#ffffff", "#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000")) +
    theme_bw() +
    theme(
        axis.text.x = element_text(colour="black",size=10),
        axis.text.y = element_text(colour="black",size=10),
        axis.title.y = element_text(colour="black", size=12, face="bold"),
        axis.title.x = element_text(color="black", size=12, face="bold"),
        legend.position="none"
        ) +
    labs(x="Partition", y="Proportion of Phenotypic Variance Explained")
dev.off()
