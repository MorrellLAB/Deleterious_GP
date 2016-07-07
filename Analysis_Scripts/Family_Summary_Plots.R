#   Read the family summary file and visualize the information

#   Set some colors for the chromosomes. This is a 'qualitative' set from
#   Rcolorbrewer
colors <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
            "#b3de69")

library(reshape2)
famsum <- read.table("Genomic_Prediction_Family_Genetic_Summary.txt", header=TRUE)

#   Make a data frame only for only the number of markers per chromosome
markers <- famsum[, c("Cycle", "NMarkers_1H", "NMarkers_2H", "NMarkers_3H",
                      "NMarkers_4H", "NMarkers_5H", "NMarkers_6H", "NMarkers_7H")]
#   Melt it for side-by-side plots
markers <- melt(markers, id="Cycle")

#   And plot it!
pdf(file="Genomic_Prediction_NMarkers.pdf", width=10.5, height=8)
boxplot(
    value~variable+Cycle,
    data=markers,
    at=c(1:7, 10:16, 19:25),
    names=rep(c("1H", "2H", "3H", "4H", "5H", "6H", "7H"), 3),
    col=colors,
    ylab="Number of Polymorphic Markers",
    cex.axis=0.85)
axis(
    side=3,
    at=c(mean(1:7), mean(10:16), mean(19:25)),
    labels=c("Cycle 1", "Cycle 2", "Cycle 3"))
dev.off()

#   Now do the mean XO plots
xo <- famsum[, c("Cycle", "MeanXO_1H", "MeanXO_2H", "MeanXO_3H", "MeanXO_4H",
                 "MeanXO_5H", "MeanXO_6H", "MeanXO_7H")]
xo <- melt(xo, id="Cycle")
pdf(file="Genomic_Prediction_XO.pdf", width=10.5, height=8)
boxplot(
    value~variable+Cycle,
    data=xo,
    at=c(1:7, 10:16, 19:25),
    names=rep(c("1H", "2H", "3H", "4H", "5H", "6H", "7H"), 3),
    col=colors,
    ylab="Mean Number of Detectable Crossovers",
    cex.axis=0.85)
axis(
    side=3,
    at=c(mean(1:7), mean(10:16), mean(19:25)),
    labels=c("Cycle1", "Cycle2", "Cycle3"))
dev.off()

#   And the marker spacing plot
spacing <- famsum[, c("Cycle", "MarkerSpacing_1H", "MarkerSpacing_2H",
                      "MarkerSpacing_3H", "MarkerSpacing_4H", "MarkerSpacing_5H",
                      "MarkerSpacing_6H", "MarkerSpacing_7H")]
spacing <- melt(spacing, id="Cycle")
pdf(file="Genomic_Prediction_MarkerSpacing.pdf", width=10.5, height=8)
boxplot(
    value~variable+Cycle,
    data=spacing,
    at=c(1:7, 10:16, 19:25),
    names=rep(c("1H", "2H", "3H", "4H", "5H", "6H", "7H"), 3),
    col=colors,
    ylab="Mean Marker Spacing (cM)",
    cex.axis=0.85)
axis(
    side=3,
    at=c(mean(1:7), mean(10:16), mean(19:25)),
    labels=c("Cycle 1", "Cycle 2", "Cycle 3"))
dev.off()
