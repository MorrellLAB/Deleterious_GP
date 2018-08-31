#   Plot the number of deleterious SNPs in each line, separated by cycle.
#   We want to make beeswarm plots
library(beeswarm)

# Read the dosages and homozygous genotype counts
hom_counts <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Burden/Homozygous_Derived_Counts.txt.gz", header=TRUE)
del_dosages <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Burden/Full_T3_Deleterious_Burden.profile.gz", header=TRUE)
# Random lines
rand <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Data/Pedigrees/Random_Lines.txt", header=FALSE)$V1)
# Selected lines
selected <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Data/Pedigrees/Selected_Lines.txt", header=FALSE)$V1)


#   Assign a color to plot the random, selected, and 'none' lines
#   They will be black, red, and grey, respectively.
hom_counts$Type <- sapply(
    hom_counts$LineID,
    function(x) {
        if(as.character(x) %in% selected) {
            return("Sel")
        }
        else if(as.character(x) %in% rand) {
            return("Ran")
        }
        else {
            return("Non")
        }
    }
    )

# Assign a cycle
cyc <- sapply(
    hom_counts$LineID,
    function(x) {
        if(grepl("MS10", x)) {
            return("C1")
        }
        else if(grepl("MS11", x)) {
            return("C2")
        }
        else if(grepl("MS12", x)) {
            return("C3")
        }
        else {
            return("C0")
        }
    }
    )
hom_counts$Cycle <- factor(cyc, levels=c("C0", "C1", "C2", "C3"))

# Do the same for dosages
del_dosages$Type <- sapply(
    del_dosages$IID,
    function(x) {
        if(as.character(x) %in% selected) {
            return("Sel")
        }
        else if(as.character(x) %in% rand) {
            return("Ran")
        }
        else {
            return("Non")
        }
    }
    )
cyc <- sapply(
    del_dosages$IID,
    function(x) {
        if(grepl("MS10", x)) {
            return("C1")
        }
        else if(grepl("MS11", x)) {
            return("C2")
        }
        else if(grepl("MS12", x)) {
            return("C3")
        }
        else {
            return("C0")
        }
    }
    )
del_dosages$Cycle <- factor(cyc, levels=c("C0", "C1", "C2", "C3"))

#   Separate the data now, so we can clearly show the random and selected lines
ran <- hom_counts[hom_counts$Type == "Ran",]
sel <- hom_counts[hom_counts$Type == "Sel",]
non <- hom_counts[hom_counts$Type == "Non",]

pdf(file="DM_By_Cycle.pdf", height=6, width=4)
par(mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
#   Plot the not selected nor random lines
beeswarm(
    non$Deleterious ~ non$Cycle,
    col="#cccccc",
    pch=19,
    cex=0.25,
    method="hex",
    ylim=c(110, 300),
    corral="wrap",
    xlab="Cycle",
    ylab="Number of Homozygous Deleterious SNPs",
    main="",
    axes=F)
beeswarm(
    ran$Deleterious ~ ran$Cycle,
    col="#333333",
    pch=19,
    cex=0.5,
    method="hex",
    ylim=c(110, 300),
    corral="wrap",
    add=TRUE,
    side=-1,
    axes=F)
beeswarm(
    sel$Deleterious ~ sel$Cycle,
    col="#990000",
    pch=19,
    cex=0.5,
    method="swarm",
    ylim=c(110, 300),
    add=TRUE,
    side=1,
    axes=F)
boxplot(
    Deleterious~Cycle + Type,
    data=droplevels(rbind(ran, sel)),
    at=c(1.75, 2.75, 3.75, 2.25, 3.25, 4.25),
    boxwex=0.2,
    lwd=2,
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
box()
dev.off()

# Do the same for average burden
ran <- del_dosages[del_dosages$Type == "Ran",]
sel <- del_dosages[del_dosages$Type == "Sel",]
non <- del_dosages[del_dosages$Type == "Non",]

pdf(file="DM_Burden_By_Cycle.pdf", height=6, width=4)
par(mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
#   Plot the not selected nor random lines
beeswarm(
    non$SCORE ~ non$Cycle,
    col="#cccccc",
    pch=19,
    cex=0.25,
    method="hex",
    ylim=c(0.1, 0.25),
    corral="wrap",
    xlab="Cycle",
    ylab="Average Burden of Deleterious Alleles",
    main="",
    axes=F)
beeswarm(
    ran$SCORE ~ ran$Cycle,
    col="#333333",
    pch=19,
    cex=0.5,
    method="hex",
    ylim=c(0.1, 0.25),
    corral="wrap",
    add=TRUE,
    side=-1,
    axes=F)
beeswarm(
    sel$SCORE ~ sel$Cycle,
    col="#990000",
    pch=19,
    cex=0.5,
    method="swarm",
    ylim=c(0.1, 0.25),
    add=TRUE,
    side=1,
    axes=F)
boxplot(
    SCORE~Cycle + Type,
    data=droplevels(rbind(ran, sel)),
    at=c(1.75, 2.75, 3.75, 2.25, 3.25, 4.25),
    boxwex=0.2,
    lwd=2,
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
box()
dev.off()
