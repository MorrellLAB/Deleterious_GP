#   Test if there are any differences in mean crossover rate in families that
#   each of the 21 genomic selection parents contributed to.

#   Read in the population XO summary
popsum <- read.table("/Volumes/Data_Disk/Dropbox/GitHub/Deleterious_GP/Data/Genotyping_Data/rQTL_Inputs/Genomic_Prediction_Family_Genetic_Summary.txt", header=T)

#   Change WD to the contribution directory
setwd("/Volumes/Data_Disk/Dropbox/Projects/DM_GenomicPrediction/Data/Progeny_Genotypes/Parental_Conribution")

#   Get all the filenames here. They contain the families that each parent contributed to.
par_filename <- list.files(".", pattern="*.txt", full.names=FALSE)

#   If we remove the '_Families.txt' bit, then we get the parental ID
fam_ids <- sapply(par_filename, strsplit, "_")
fam_ids <- as.vector(unlist(lapply(fam_ids, "[[", 1)))

#   Read in all the files
par_fams <- sapply(
    fam_ids,
    function(x) {
        y <- read.table(paste(x, "_Families.txt", sep=""), header=F)
        y <- as.character(y$V1)
        return(y)
        }
    )

#   Then, over each parent, get the meanXO observed in each of its descendant families
par_fam_XO <- lapply(par_fams,
    function(x) {
        return(popsum[popsum$FamilyID %in% x, "MeanXO"])
    })
#   Build a vector of family names to associate with mean XO numbers
fam_names <- as.character(unlist(sapply(names(par_fam_XO),
    function(x) {
        return(rep(x, length(par_fam_XO[[x]])))
    })))
#   And flatten the list of XO numbers
xo_fam <- unlist(par_fam_XO, use.names=F)

#   Drop them into a data frame
xo_data <- data.frame(
    ParName=fam_names,
    XOPerFam=xo_fam)

#   Make a plot
pdf("Parental_XO.pdf", 10, 6)
boxplot(XOPerFam ~ ParName, data=xo_data, main="Mean Number of XO Observed Per Family", ylab="Mean XO Detected Per Family", xlab="", las=2)
dev.off()

#   And do an ANOVA
xo_aov <- anova(lm(XOPerFam~ParName, data=xo_data))
print(xo_aov)
