#   Make side-by-side boxplots for selected and random lines, for each
#   cycle of genomewide selection.

library(ggplot2)
library(reshape)

setwd("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Data/Phenotypic_Data")
#   Read the DON and yield data
adjusted_phenotypes <- read.csv("Adjusted_Phenotypic_Data.csv", header=T)
don_data <- adjusted_phenotypes[,c("line_name", "cycle", "type", "DON.BLUE_mv")]
yield_data <- adjusted_phenotypes[,c("line_name", "cycle", "type", "yld.BLUE_mv")]
height_data <- adjusted_phenotypes[,c("line_name", "cycle", "type", "height.BLUE_raw")]

#   Remove NA values
don_data <- don_data[!is.na(don_data$type),]
yield_data <- yield_data[!is.na(yield_data$type),]
height_data <- height_data[!is.na(height_data$type),]

#   Put the factors in order for plotting: Check, C0, C1, C2, C3
don_data$cycle <- factor(don_data$cycle, levels=c("chk", "C0", "C1", "C2", "C3"))
yield_data$cycle <- factor(yield_data$cycle, levels=c("chk", "C0", "C1", "C2", "C3"))
height_data$cycle <- factor(height_data$cycle, levels=c("chk", "C0", "C1", "C2", "C3"))

#   Plot the DON data
pdf(
    file="Adjusted_DON_Boxplots_BLUE.pdf",
    width=2.75,
    height=4)
plt <- ggplot(
    don_data,
    aes(x=factor(cycle), y=DON.BLUE_mv, fill=factor(type)))
plt + geom_boxplot(position=position_dodge(width=0.75)) +
    scale_fill_manual(
        values=c("white", "white", "white", "grey", "grey"),
        name="type") +
    theme_bw() +
    scale_color_manual(
        values=c("black", "black", "black", "red", "red"),
        name="type") +
    scale_x_discrete(
        breaks=c("chk", "C0", "C1", "C2", "C3", "C3"),
        labels=c("Chk", "Par", "C1", "C2", "C3", "C3")) +
    theme(
        legend.position="none",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()
        ) +
    labs(x="Cycle", y="DON Concentration BLUE (ppm)")
dev.off()

#   Plot the yield data
pdf(
    file="Adjusted_Yield_Boxplots_BLUE.pdf",
    width=2.75,
    height=4)
plt <- ggplot(
    yield_data,
    aes(x=factor(cycle), y=yld.BLUE_mv, fill=factor(type)))
plt + geom_boxplot(position=position_dodge(width=0.75)) +
    scale_fill_manual(
        values=c("white", "white", "white", "grey", "grey"),
        name="type") +
    theme_bw() +
    scale_color_manual(
        values=c("black", "black", "black", "red", "red"),
        name="type") +
    scale_x_discrete(
        breaks=c("chk", "C0", "C1", "C2", "C3"),
        labels=c("Chk", "Par", "C1", "C2", "C3")) +
    theme(
        legend.position="none",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()
        ) +
    labs(x="Cycle", y="Yield BLUE (kg/ha)")
dev.off()

pdf(
    file="Adjusted_Height_Boxplots_BLUE.pdf",
    width=2.75,
    height=4)
plt <- ggplot(
    height_data,
    aes(x=factor(cycle), y=height.BLUE_raw, fill=factor(type)))
plt + geom_boxplot(position=position_dodge(width=0.75)) +
    scale_fill_manual(
        values=c("white", "white", "white", "grey", "grey"),
        name="type") +
    theme_bw() +
    scale_color_manual(
        values=c("black", "black", "black", "red", "red"),
        name="type") +
    scale_x_discrete(
        breaks=c("chk", "C0", "C1", "C2", "C3"),
        labels=c("Chk", "Par", "C1", "C2", "C3")) +
    theme(
        legend.position="none",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()
        ) +
    labs(x="Cycle", y="Height BLUE (cm)")
dev.off()
