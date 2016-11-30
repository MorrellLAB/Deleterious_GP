#   Calculate the correlation between the number of deleterious SNPs and the
#   phenotype.

#   Read in the data files
n_del <- read.table("/Volumes/Data_Disk/Dropbox/Projects/DM_GenomicPrediction/Results/DM_Selection/Deleterious_Counts_By_Line.txt", header=T)
phen <- read.csv("/Volumes/Data_Disk/Dropbox/Projects/DM_GenomicPrediction/Data/Phenotypic_Data/Adjusted_Phenotypic_Data.csv", header=T)

#   Make the Line ID variable consistent for merging
names(n_del) <- c("Cycle", "line_name", "Count")
#   And merge them
dm_phen <- merge(n_del, phen, by="line_name", all.x=F)
#   Remove lines without phenotypic data
dm_phen <- dm_phen[!is.na(dm_phen$yld.BLUE_mv), ]

#   correlate the number of deleterious SNPs with yield
r_yld <- cor(dm_phen$yld.BLUE_mv, dm_phen$Count)
#   do a two-sided t-test for significance
cor.test(dm_phen$yld.BLUE_mv, dm_phen$Count, alternative="two.sided", conf.level=0.99)

#   Make a plot of the phenotype against the number of deleterious SNPs
pdf(file="DM_Yield.pdf", 6, 6)
plot(
    dm_phen$yld.BLUE_mv ~ dm_phen$Count,
    pch=19,
    cex=0.75,
    xlab="Number of Putatively Deleterious SNPs",
    ylab="Yield BLUE (kg/ha)")
dev.off()

#   And the same with DON
cor(dm_phen$DON.BLUE_mv, dm_phen$Count)
cor.test(dm_phen$DON.BLUE_mv, dm_phen$Count, alternative="two.sided", conf.level=0.99)
pdf(file="DM_DON.pdf", 6, 6)
plot(
    dm_phen$DON.BLUE_mv ~ dm_phen$Count,
    pch=19,
    cex=0.75,
    xlab="Number of Putatively Deleterious SNPs",
    ylab="DON BLUE (ppm)")
dev.off()
