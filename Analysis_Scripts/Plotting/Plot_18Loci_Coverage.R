#   Script to plot coverage values for each of the "18 Loci"

#   Color palette, from https://www.r-bloggers.com/the-paul-tol-21-color-salute/
tol21rainbow <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

#   Read in the data table
coverage <- read.table("/Volumes/Data_Disk/Dropbox/GitHub/Deleterious_GP/Data/18Loci_Coverage.txt", header=T)

genes <- as.character(levels(coverage$Gene))

#   Make a plot for each gene
pdf(file="18Loci_Coverage.pdf", width=20, height=20)
par(mfrow=c(5, 4))
sapply(
    genes,
    function(x) {
        genecov <- coverage[coverage$Gene == x,]
        samplenames <- colnames(coverage)[3:21]
        #   Get the highest observed coverage for the gene at any point - this
        #   is the boundary of our plot
        maxcov <- max(as.matrix(genecov[,3:21]))
        #   And get the length of the gene, add for legend padding
        maxpos <- max(genecov$Pos) + 200
        #   Open a blank plot. Again, this is hardcoded and should change, but
        #   it is easier for now
        plot(c(0, 0), ylim=c(0, maxcov), xlim=c(0, maxpos), type="n", xlab="Gene Position (bp)", ylab="Number of Reads", main=paste("Coverage Over", x))
        sapply(
            seq_along(samplenames),
            function(y) {
                lines(genecov[,samplenames[y]] ~ genecov$Pos, col=tol21rainbow[y], lwd=1)
            }
            )
        legend(
            "topright",
            gsub("X", "", samplenames),
            col=tol21rainbow,
            lwd=1,
            cex=0.55)
    }
    )
dev.off()
