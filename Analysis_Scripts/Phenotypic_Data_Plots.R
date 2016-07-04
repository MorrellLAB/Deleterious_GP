#   Make side-by-side boxplots for selected and random lines, for each
#   cycle of genomewide selection.

library(ggplot2)
library(reshape)

setwd("/Users/tomkono/DataDisk/Dropbox/Projects/DM_GenomicPrediction/Data/Phenotypic_Data")
#   Read the DON and yield data
don_data <- read.csv("Raw_DON_Data.csv", header=TRUE)
yield_data <- read.csv("Raw_Yield_Data.csv", header=TRUE)

#   Get "selected" and "random" lines
don_sel <- don_data[don_data$Type == "Sel",]
don_ran <- don_data[don_data$Type == "Ran",]
yield_sel <- yield_data[don_data$Type == "Sel",]
yield_ran <- yield_data[yield_data$Type == "Ran",]

#   Then, we want to separate them by cycle
don_sel_c1 <- don_sel[don_sel$Cycle == "C1",]
don_sel_c2 <- don_sel[don_sel$Cycle == "C2",]
don_sel_c3 <- don_sel[don_sel$Cycle == "C3",]
don_ran_c1 <- don_ran[don_ran$Cycle == "C1",]
don_ran_c2 <- don_ran[don_ran$Cycle == "C2",]
don_ran_c3 <- don_ran[don_ran$Cycle == "C3",]

yield_sel_c1 <- yield_sel[yield_sel$Cycle == "C1",]
yield_sel_c2 <- yield_sel[yield_sel$Cycle == "C2",]
yield_sel_c3 <- yield_sel[yield_sel$Cycle == "C3",]
yield_ran_c1 <- yield_ran[yield_ran$Cycle == "C1",]
yield_ran_c2 <- yield_ran[yield_ran$Cycle == "C2",]
yield_ran_c3 <- yield_ran[yield_ran$Cycle == "C3",]

#   Next, we want to average the values for each unique line name
don_sel_c1_avg <- sapply(
    unique(as.character(don_sel_c1$line_name)),
    function(x) {
        return(
            mean(don_sel_c1$DON[don_sel_c1$line_name == x], na.rm=TRUE)
            )
        })
don_sel_c2_avg <- sapply(
    unique(as.character(don_sel_c2$line_name)),
    function(x) {
        return(
            mean(don_sel_c2$DON[don_sel_c2$line_name == x], na.rm=TRUE)
            )
        })
don_sel_c3_avg <- sapply(
    unique(as.character(don_sel_c3$line_name)),
    function(x) {
        return(
            mean(don_sel_c3$DON[don_sel_c3$line_name == x], na.rm=TRUE)
            )
        })
don_ran_c1_avg <- sapply(
    unique(as.character(don_ran_c1$line_name)),
    function(x) {
        return(
            mean(don_ran_c1$DON[don_ran_c1$line_name == x], na.rm=TRUE)
            )
        })
don_ran_c2_avg <- sapply(
    unique(as.character(don_ran_c2$line_name)),
    function(x) {
        return(
            mean(don_ran_c2$DON[don_ran_c2$line_name == x], na.rm=TRUE)
            )
        })
don_ran_c3_avg <- sapply(
    unique(as.character(don_ran_c3$line_name)),
    function(x) {
        return(
            mean(don_ran_c3$DON[don_ran_c3$line_name == x], na.rm=TRUE)
            )
        })

# #   Remove missing values
# don_sel_c1_avg <- don_sel_c1_avg[!is.na(don_sel_c1_avg)]
# don_sel_c2_avg <- don_sel_c2_avg[!is.na(don_sel_c2_avg)]
# don_sel_c3_avg <- don_sel_c3_avg[!is.na(don_sel_c3_avg)]
# don_ran_c1_avg <- don_ran_c1_avg[!is.na(don_ran_c1_avg)]
# don_ran_c2_avg <- don_ran_c2_avg[!is.na(don_ran_c2_avg)]
# don_ran_c3_avg <- don_ran_c3_avg[!is.na(don_ran_c3_avg)]


#   Next, we want to average the yield values for each unique line name
yield_sel_c1_avg <- sapply(
    unique(as.character(yield_sel_c1$line_name)),
    function(x) {
        return(
            mean(yield_sel_c1$yld_bu[yield_sel_c1$line_name == x], na.rm=TRUE)
            )
        })
yield_sel_c2_avg <- sapply(
    unique(as.character(yield_sel_c2$line_name)),
    function(x) {
        return(
            mean(yield_sel_c2$yld_bu[yield_sel_c2$line_name == x], na.rm=TRUE)
            )
        })
yield_sel_c3_avg <- sapply(
    unique(as.character(yield_sel_c3$line_name)),
    function(x) {
        return(
            mean(yield_sel_c3$yld_bu[yield_sel_c3$line_name == x], na.rm=TRUE)
            )
        })
yield_ran_c1_avg <- sapply(
    unique(as.character(yield_ran_c1$line_name)),
    function(x) {
        return(
            mean(yield_ran_c1$yld_bu[yield_ran_c1$line_name == x], na.rm=TRUE)
            )
        })
yield_ran_c2_avg <- sapply(
    unique(as.character(yield_ran_c2$line_name)),
    function(x) {
        return(
            mean(yield_ran_c2$yld_bu[yield_ran_c2$line_name == x], na.rm=TRUE)
            )
        })
yield_ran_c3_avg <- sapply(
    unique(as.character(yield_ran_c3$line_name)),
    function(x) {
        return(
            mean(yield_ran_c3$yld_bu[yield_ran_c3$line_name == x], na.rm=TRUE)
            )
        })

# #   Remove missing values
# yield_sel_c1_avg <- yield_sel_c1_avg[!is.na(yield_sel_c1_avg)]
# yield_sel_c2_avg <- yield_sel_c2_avg[!is.na(yield_sel_c2_avg)]
# yield_sel_c3_avg <- yield_sel_c3_avg[!is.na(yield_sel_c3_avg)]
# yield_ran_c1_avg <- yield_ran_c1_avg[!is.na(yield_ran_c1_avg)]
# yield_ran_c2_avg <- yield_ran_c2_avg[!is.na(yield_ran_c2_avg)]
# yield_ran_c3_avg <- yield_ran_c3_avg[!is.na(yield_ran_c3_avg)]

#   Put DON and yield data into data frames. We have to make sure that we
#   access them in the proper order.
don_sel_c1_line_names <- as.character(unique(don_sel_c1$line_name))
don_sel_c1_line_names <- don_sel_c1_line_names[!is.na(don_sel_c1_line_names)]
don_ran_c1_line_names <- as.character(unique(don_ran_c1$line_name))
don_ran_c1_line_names <- don_ran_c1_line_names[!is.na(don_ran_c1_line_names)]

don_sel_c2_line_names <- as.character(unique(don_sel_c2$line_name))
don_sel_c2_line_names <- don_sel_c2_line_names[!is.na(don_sel_c2_line_names)]
don_ran_c2_line_names <- as.character(unique(don_ran_c2$line_name))
don_ran_c2_line_names <- don_ran_c2_line_names[!is.na(don_ran_c2_line_names)]

don_sel_c3_line_names <- as.character(unique(don_sel_c3$line_name))
don_sel_c3_line_names <- don_sel_c3_line_names[!is.na(don_sel_c3_line_names)]
don_ran_c3_line_names <- as.character(unique(don_ran_c3$line_name))
don_ran_c3_line_names <- don_ran_c3_line_names[!is.na(don_ran_c3_line_names)]

#   For some reason, there are 51 selected lines in Cycle3, so we append an NA
#   to selected C1 and C2.
don_plot_data <- data.frame(
    Type=c(rep("Selected", length(don_sel_c1_line_names)+1), rep("Random", length(don_ran_c1_line_names))),
    C1=c(don_sel_c1_avg[don_sel_c1_line_names], NA, don_ran_c1_avg[don_ran_c1_line_names]),
    C2=c(don_sel_c2_avg[don_sel_c2_line_names], NA, don_ran_c2_avg[don_ran_c2_line_names]),
    C3=c(don_sel_c3_avg[don_sel_c3_line_names], don_ran_c3_avg[don_ran_c3_line_names])
    )
#   Melt it for ggplot
don_plot_data <- melt(don_plot_data)

#   We want to do the same for yield
yield_sel_c1_line_names <- as.character(unique(yield_sel_c1$line_name))
yield_sel_c1_line_names <- yield_sel_c1_line_names[!is.na(yield_sel_c1_line_names)]
yield_ran_c1_line_names <- as.character(unique(yield_ran_c1$line_name))
yield_ran_c1_line_names <- yield_ran_c1_line_names[!is.na(yield_ran_c1_line_names)]

yield_sel_c2_line_names <- as.character(unique(yield_sel_c2$line_name))
yield_sel_c2_line_names <- yield_sel_c2_line_names[!is.na(yield_sel_c2_line_names)]
yield_ran_c2_line_names <- as.character(unique(yield_ran_c2$line_name))
yield_ran_c2_line_names <- yield_ran_c2_line_names[!is.na(yield_ran_c2_line_names)]

yield_sel_c3_line_names <- as.character(unique(yield_sel_c3$line_name))
yield_sel_c3_line_names <- yield_sel_c3_line_names[!is.na(yield_sel_c3_line_names)]
yield_ran_c3_line_names <- as.character(unique(yield_ran_c3$line_name))
yield_ran_c3_line_names <- yield_ran_c3_line_names[!is.na(yield_ran_c3_line_names)]

#   Similarly, the maximum number selected was 96, and the maxiumum number
#   of random was 51. We pad with NAs accordingly.
yield_plot_data <- data.frame(
    Type=c(rep("Selected", 96), rep("Random", 51)),
    C1=c(yield_sel_c1_avg[yield_sel_c1_line_names], NA, NA, yield_ran_c1_avg[yield_ran_c1_line_names], NA),
    C2=c(yield_sel_c2_avg[yield_sel_c2_line_names], yield_ran_c2_avg[yield_ran_c2_line_names], NA),
    C3=c(yield_sel_c3_avg[yield_sel_c3_line_names], yield_ran_c3_avg[yield_ran_c3_line_names])
    )
yield_plot_data <- melt(yield_plot_data)

#   Then create a side-by-side boxplot for DON, and one for yield.
pdf(
    file="DON_Boxplots.pdf",
    width=4.75,
    height=7)
plt <- ggplot(
    don_plot_data,
    aes(x=factor(variable), y=value))
plt + geom_boxplot(aes(fill=factor(Type))) +
    scale_fill_manual(
        values=c("white", "red"),
        name="Type") +
    theme_bw() +
    theme(
        axis.text.x = element_text(colour="black",size=14,face="bold"),
        axis.text.y = element_text(colour="black",size=14,face="bold"),
        axis.title.y = element_text(colour="black", size=16, face="bold"),
        axis.title.x = element_text(color="black", size=16, face="bold"),
        legend.position="none"
        ) +
    labs(x="Cycle", y="DON Concentration (ppm)")
dev.off()

pdf(
    file="Yield_Boxplots.pdf",
    width=4.75,
    height=7)
plt <- ggplot(
    yield_plot_data,
    aes(x=factor(variable), y=value*53.7996))
plt + geom_boxplot(aes(fill=factor(Type))) +
    scale_fill_manual(
        values=c("white", "red"),
        name="Type") +
    theme_bw() +
    theme(
        axis.text.x = element_text(colour="black",size=14,face="bold"),
        axis.text.y = element_text(colour="black",size=14,face="bold"),
        axis.title.y = element_text(colour="black", size=16, face="bold"),
        axis.title.x = element_text(color="black", size=16, face="bold"),
        legend.position="none"
        ) +
    labs(x="Cycle", y="Yield (kg/ha)")
dev.off()
