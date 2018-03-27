#   Plot the number of deleterious SNPs in each line, separated by cycle.
#   We want to make beeswarm plots
library(beeswarm)

# Read the dosages
dos <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Summaries/GP_Deleterious_Homozygotes.txt", header=TRUE)
# Random lines
rand <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Data/Pedigrees/Random_Lines.txt", header=FALSE)$V1)
# Selected lines
selected <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Data/Pedigrees/Selected_Lines.txt", header=FALSE)$V1)


#   Assign a color to plot the random, selected, and 'none' lines
#   They will be black, red, and grey, respectively.
type <- sapply(
    dos$line_name,
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

dos$Type <- type

# Assign a cycle
cyc <- sapply(
    dos$line_name,
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

dos$Cycle <- factor(cyc, levels=c("C0", "C1", "C2", "C3"))


#   Separate the data now, so we can clearly show the random and selected lines
ran <- dos[dos$Type == "Ran",]
sel <- dos[dos$Type == "Sel",]
non <- dos[dos$Type == "Non",]

pdf(file="DM_By_Cycle.pdf", height=6, width=4)
par(mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
#   Plot the not selected nor random lines
beeswarm(
    non$Homozygotes ~ non$Cycle,
    col="#cccccc",
    pch=19,
    cex=0.25,
    method="hex",
    ylim=c(0, 300),
    corral="wrap",
    xlab="Cycle",
    ylab="Number of Homozygous Deleterious SNPs",
    main="",
    axes=F)
beeswarm(
    ran$Homozygotes ~ ran$Cycle,
    col="#333333",
    pch=19,
    cex=0.5,
    method="hex",
    ylim=c(0, 300),
    corral="wrap",
    add=TRUE,
    side=-1,
    axes=F)
beeswarm(
    sel$Homozygotes ~ sel$Cycle,
    col="#990000",
    pch=19,
    cex=0.5,
    method="swarm",
    ylim=c(0, 300),
    add=TRUE,
    side=1,
    axes=F)
boxplot(
    Homozygotes~Cycle + Type,
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
