#   Create a plot that summarizes the total number of reads and those that
#   are soft-masked by BWA as a stacked bar chart.

args <- commandArgs(TRUE)

#   Read the text table
counts <- read.table(args[1], header=TRUE)

#   Calculate the number of reads that were soft-masked
counts$Soft <- counts$All - counts$Nonmasked

#   Then make the data for a stacked barplot
plotdata <- as.matrix(counts[, c("Nonmasked", "Soft")]/1000000)
rownames(plotdata) <- counts$SampleName

#   Then plot it
pdf(file="Clipping_Summary.pdf", width=10, height=12)
par(mfrow=c(2,1))
at <- barplot(
    t(plotdata),
    col=c("white", "grey"),
    legend=c("Non-masked", "Softmasked"),
    xaxt="n",
    ylab="Read Counts (Millions)",
    main="Cleaned Read Clipping Summary\nTotal Counts")
axis(
    side=1,
    at=at,
    lab=counts$SampleName,
    las=2,
    cex.axis=0.75)

plotdata <- as.matrix((counts[, c("Nonmasked", "Soft")]/counts$All)*100)
rownames(plotdata) <- counts$SampleName
at <- barplot(
    t(plotdata),
    col=c("white", "grey"),
    xaxt="n",
    ylab="Percent Mapped Reads",
    main="Percentage")
axis(
    side=1,
    at=at,
    lab=counts$SampleName,
    las=2,
    cex.axis=0.75)
dev.off()
