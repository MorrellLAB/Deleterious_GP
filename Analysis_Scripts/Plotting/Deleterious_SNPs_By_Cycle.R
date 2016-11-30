#   Plot the number of deleterious SNPs in each line, separated by cycle.
#   We want to make beeswarm plots
library(beeswarm)
#   Read in the data file with the raw counts
dm_sel <- read.table("/Volumes/Data_Disk/Dropbox/Projects/DM_GenomicPrediction/Results/DM_Selection/Deleterious_Counts_By_Line.txt", header=TRUE)
#   Read in the lines that are selected
selected <- read.table("/Volumes/Data_Disk/Dropbox/Projects/DM_GenomicPrediction/Results/DM_Selection/Selected_Lines.txt", header=FALSE)
#   And the lines that were random
rand <- read.table("/Volumes/Data_Disk/Dropbox/Projects/DM_GenomicPrediction/Results/DM_Selection/Random_Lines.txt", header=FALSE)

#   Assign a color to plot the random, selected, and 'none' lines
#   They will be black, red, and grey, respectively.
type <- sapply(
    dm_sel$LineID,
    function(x) {
        if(as.character(x) %in% selected$V1) {
            return("Sel")
        }
        else if(as.character(x) %in% rand$V1) {
            return("Ran")
        }
        else {
            return("Non")
        }
    }
    )

dm_sel$Type <- type

#   Separate the data now, so we can clearly show the random and selected lines
ran <- dm_sel[dm_sel$Type == "Ran",]
sel <- dm_sel[dm_sel$Type == "Sel",]
non <- dm_sel[dm_sel$Type == "Non",]

pdf(file="DM_By_Cycle.pdf", height=6, width=8)
#   Plot the not selected nor random lines
beeswarm(
    non$Count ~ non$Cycle,
    col="#cccccc",
    pch=19,
    cex=0.35,
    method="hex",
    ylim=c(600, 750),
    xlab="Cycle",
    ylab="Number of Deleterious SNPs",
    main="Deleterious SNPs Over Time",
    axes=F)
beeswarm(
    ran$Count ~ ran$Cycle,
    col="#333333",
    pch=19,
    cex=0.5,
    method="hex",
    ylim=c(600, 750),
    add=TRUE,
    side=-1,
    axes=F)
beeswarm(
    sel$Count ~ sel$Cycle,
    col="#aa0000",
    pch=19,
    cex=0.5,
    method="hex",
    ylim=c(600, 750),
    add=TRUE,
    side=1,
    axes=F)
boxplot(
    Count~Cycle + Type,
    data=droplevels(rbind(ran, sel)),
    at=c(1.75, 2.75, 3.75, 2.25, 3.25, 4.25),
    boxwex=0.2,
    lwd=1,
    border=c("black", "black", "black", "red", "red", "red"),
    axes=F,
    add=TRUE)
legend("topright", pch=19, col=c("black", "red"), legend=c("Random", "Selected"))
axis(
    side=1,
    at=c(1, 2, 3, 4),
    labels=c("Parents", "C1", "C2", "C3"))
axis(
    side=2)
dev.off()
