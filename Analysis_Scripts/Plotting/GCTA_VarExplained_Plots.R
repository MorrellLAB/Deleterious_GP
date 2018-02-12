#   Make a side-by-side boxplot to show the distribution of variance explained
#   by random subsets of 101 markers across various partitions of marker
#   classes.

#   Set working directory
setwd("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/GCTA")

# Read the data files. There are a lot.
yield.null.250 <- read.table("Analysis/Yield_GCTA/250/null/VG_VP.txt", header=FALSE)$V1
yield.nonc.250 <- read.table("Analysis/Yield_GCTA/250/noncoding/VG_VP.txt", header=FALSE)$V1
yield.syn.250 <- read.table("Analysis/Yield_GCTA/250/synonymous/VG_VP.txt", header=FALSE)$V1
yield.nonsyn.250 <- read.table("Analysis/Yield_GCTA/250/nonsynonymous/VG_VP.txt", header=FALSE)$V1
yield.del.250 <- read.table("Analysis/Yield_GCTA/250/deleterious/VG_VP.txt", header=FALSE)$V1

yield.null.500 <- read.table("Analysis/Yield_GCTA/500/null/VG_VP.txt", header=FALSE)$V1
yield.nonc.500 <- read.table("Analysis/Yield_GCTA/500/noncoding/VG_VP.txt", header=FALSE)$V1
yield.syn.500 <- read.table("Analysis/Yield_GCTA/500/synonymous/VG_VP.txt", header=FALSE)$V1
yield.nonsyn.500 <- read.table("Analysis/Yield_GCTA/500/nonsynonymous/VG_VP.txt", header=FALSE)$V1
yield.del.500 <- read.table("Analysis/Yield_GCTA/500/deleterious/VG_VP.txt", header=FALSE)$V1

yield.null.750 <- read.table("Analysis/Yield_GCTA/750/null/VG_VP.txt", header=FALSE)$V1
yield.nonc.750 <- read.table("Analysis/Yield_GCTA/750/noncoding/VG_VP.txt", header=FALSE)$V1
yield.syn.750 <- read.table("Analysis/Yield_GCTA/750/synonymous/VG_VP.txt", header=FALSE)$V1
yield.nonsyn.750 <- read.table("Analysis/Yield_GCTA/750/nonsynonymous/VG_VP.txt", header=FALSE)$V1
yield.del.750 <- read.table("Analysis/Yield_GCTA/750/deleterious/VG_VP.txt", header=FALSE)$V1

yield.null.1000 <- read.table("Analysis/Yield_GCTA/1000/null/VG_VP.txt", header=FALSE)$V1
yield.nonc.1000 <- read.table("Analysis/Yield_GCTA/1000/noncoding/VG_VP.txt", header=FALSE)$V1
yield.syn.1000 <- read.table("Analysis/Yield_GCTA/1000/synonymous/VG_VP.txt", header=FALSE)$V1
yield.nonsyn.1000 <- read.table("Analysis/Yield_GCTA/1000/nonsynonymous/VG_VP.txt", header=FALSE)$V1
yield.del.1000 <- read.table("Analysis/Yield_GCTA/1000/deleterious/VG_VP.txt", header=FALSE)$V1

yield.null.1250 <- read.table("Analysis/Yield_GCTA/1250/null/VG_VP.txt", header=FALSE)$V1
yield.nonc.1250 <- read.table("Analysis/Yield_GCTA/1250/noncoding/VG_VP.txt", header=FALSE)$V1
yield.syn.1250 <- read.table("Analysis/Yield_GCTA/1250/synonymous/VG_VP.txt", header=FALSE)$V1
yield.nonsyn.1250 <- read.table("Analysis/Yield_GCTA/1250/nonsynonymous/VG_VP.txt", header=FALSE)$V1
yield.del.1250 <- read.table("Analysis/Yield_GCTA/1250/deleterious/VG_VP.txt", header=FALSE)$V1

yield.null.1500 <- read.table("Analysis/Yield_GCTA/1500/null/VG_VP.txt", header=FALSE)$V1
yield.nonc.1500 <- read.table("Analysis/Yield_GCTA/1500/noncoding/VG_VP.txt", header=FALSE)$V1
yield.syn.1500 <- read.table("Analysis/Yield_GCTA/1500/synonymous/VG_VP.txt", header=FALSE)$V1
yield.nonsyn.1500 <- read.table("Analysis/Yield_GCTA/1500/nonsynonymous/VG_VP.txt", header=FALSE)$V1
yield.del.1500 <- read.table("Analysis/Yield_GCTA/1500/deleterious/VG_VP.txt", header=FALSE)$V1

don.null.250 <- read.table("Analysis/DON_GCTA/250/null/VG_VP.txt", header=FALSE)$V1
don.nonc.250 <- read.table("Analysis/DON_GCTA/250/noncoding/VG_VP.txt", header=FALSE)$V1
don.syn.250 <- read.table("Analysis/DON_GCTA/250/synonymous/VG_VP.txt", header=FALSE)$V1
don.nonsyn.250 <- read.table("Analysis/DON_GCTA/250/nonsynonymous/VG_VP.txt", header=FALSE)$V1
don.del.250 <- read.table("Analysis/DON_GCTA/250/deleterious/VG_VP.txt", header=FALSE)$V1

don.null.500 <- read.table("Analysis/DON_GCTA/500/null/VG_VP.txt", header=FALSE)$V1
don.nonc.500 <- read.table("Analysis/DON_GCTA/500/noncoding/VG_VP.txt", header=FALSE)$V1
don.syn.500 <- read.table("Analysis/DON_GCTA/500/synonymous/VG_VP.txt", header=FALSE)$V1
don.nonsyn.500 <- read.table("Analysis/DON_GCTA/500/nonsynonymous/VG_VP.txt", header=FALSE)$V1
don.del.500 <- read.table("Analysis/DON_GCTA/500/deleterious/VG_VP.txt", header=FALSE)$V1

don.null.750 <- read.table("Analysis/DON_GCTA/750/null/VG_VP.txt", header=FALSE)$V1
don.nonc.750 <- read.table("Analysis/DON_GCTA/750/noncoding/VG_VP.txt", header=FALSE)$V1
don.syn.750 <- read.table("Analysis/DON_GCTA/750/synonymous/VG_VP.txt", header=FALSE)$V1
don.nonsyn.750 <- read.table("Analysis/DON_GCTA/750/nonsynonymous/VG_VP.txt", header=FALSE)$V1
don.del.750 <- read.table("Analysis/DON_GCTA/750/deleterious/VG_VP.txt", header=FALSE)$V1

don.null.1000 <- read.table("Analysis/DON_GCTA/1000/null/VG_VP.txt", header=FALSE)$V1
don.nonc.1000 <- read.table("Analysis/DON_GCTA/1000/noncoding/VG_VP.txt", header=FALSE)$V1
don.syn.1000 <- read.table("Analysis/DON_GCTA/1000/synonymous/VG_VP.txt", header=FALSE)$V1
don.nonsyn.1000 <- read.table("Analysis/DON_GCTA/1000/nonsynonymous/VG_VP.txt", header=FALSE)$V1
don.del.1000 <- read.table("Analysis/DON_GCTA/1000/deleterious/VG_VP.txt", header=FALSE)$V1

don.null.1250 <- read.table("Analysis/DON_GCTA/1250/null/VG_VP.txt", header=FALSE)$V1
don.nonc.1250 <- read.table("Analysis/DON_GCTA/1250/noncoding/VG_VP.txt", header=FALSE)$V1
don.syn.1250 <- read.table("Analysis/DON_GCTA/1250/synonymous/VG_VP.txt", header=FALSE)$V1
don.nonsyn.1250 <- read.table("Analysis/DON_GCTA/1250/nonsynonymous/VG_VP.txt", header=FALSE)$V1
don.del.1250 <- read.table("Analysis/DON_GCTA/1250/deleterious/VG_VP.txt", header=FALSE)$V1

don.null.1500 <- read.table("Analysis/DON_GCTA/1500/null/VG_VP.txt", header=FALSE)$V1
don.nonc.1500 <- read.table("Analysis/DON_GCTA/1500/noncoding/VG_VP.txt", header=FALSE)$V1
don.syn.1500 <- read.table("Analysis/DON_GCTA/1500/synonymous/VG_VP.txt", header=FALSE)$V1
don.nonsyn.1500 <- read.table("Analysis/DON_GCTA/1500/nonsynonymous/VG_VP.txt", header=FALSE)$V1
don.del.1500 <- read.table("Analysis/DON_GCTA/1500/deleterious/VG_VP.txt", header=FALSE)$V1


# Define a function to "sanitize" the values. Values that are greater than
# 0.999 or less than 0 will be trimmed.
clean <- function(x) {
    cleaned <- x[x>0 & x<=0.999]
    return(cleaned)
}

# Clean the data
yield.null.250 <- clean(yield.null.250)
yield.nonc.250 <- clean(yield.nonc.250)
yield.syn.250 <- clean(yield.syn.250)
yield.nonsyn.250 <- clean(yield.nonsyn.250)
yield.del.250 <- clean(yield.del.250)
yield.null.500 <- clean(yield.null.500)
yield.nonc.500 <- clean(yield.nonc.500)
yield.syn.500 <- clean(yield.syn.500)
yield.nonsyn.500 <- clean(yield.nonsyn.500)
yield.del.500 <- clean(yield.del.500)
yield.null.750 <- clean(yield.null.750)
yield.nonc.750 <- clean(yield.nonc.750)
yield.syn.750 <- clean(yield.syn.750)
yield.nonsyn.750 <- clean(yield.nonsyn.750)
yield.del.750 <- clean(yield.del.750)
yield.null.1000 <- clean(yield.null.1000)
yield.nonc.1000 <- clean(yield.nonc.1000)
yield.syn.1000 <- clean(yield.syn.1000)
yield.nonsyn.1000 <- clean(yield.nonsyn.1000)
yield.del.1000 <- clean(yield.del.1000)
yield.null.1250 <- clean(yield.null.1250)
yield.nonc.1250 <- clean(yield.nonc.1250)
yield.syn.1250 <- clean(yield.syn.1250)
yield.nonsyn.1250 <- clean(yield.nonsyn.1250)
yield.del.1250 <- clean(yield.del.1250)
yield.null.1500 <- clean(yield.null.1500)
yield.nonc.1500 <- clean(yield.nonc.1500)
yield.syn.1500 <- clean(yield.syn.1500)
yield.nonsyn.1500 <- clean(yield.nonsyn.1500)
yield.del.1500 <- clean(yield.del.1500)
don.null.250 <- clean(don.null.250)
don.nonc.250 <- clean(don.nonc.250)
don.syn.250 <- clean(don.syn.250)
don.nonsyn.250 <- clean(don.nonsyn.250)
don.del.250 <- clean(don.del.250)
don.null.500 <- clean(don.null.500)
don.nonc.500 <- clean(don.nonc.500)
don.syn.500 <- clean(don.syn.500)
don.nonsyn.500 <- clean(don.nonsyn.500)
don.del.500 <- clean(don.del.500)
don.null.750 <- clean(don.null.750)
don.nonc.750 <- clean(don.nonc.750)
don.syn.750 <- clean(don.syn.750)
don.nonsyn.750 <- clean(don.nonsyn.750)
don.del.750 <- clean(don.del.750)
don.null.1000 <- clean(don.null.1000)
don.nonc.1000 <- clean(don.nonc.1000)
don.syn.1000 <- clean(don.syn.1000)
don.nonsyn.1000 <- clean(don.nonsyn.1000)
don.del.1000 <- clean(don.del.1000)
don.null.1250 <- clean(don.null.1250)
don.nonc.1250 <- clean(don.nonc.1250)
don.syn.1250 <- clean(don.syn.1250)
don.nonsyn.1250 <- clean(don.nonsyn.1250)
don.del.1250 <- clean(don.del.1250)
don.null.1500 <- clean(don.null.1500)
don.nonc.1500 <- clean(don.nonc.1500)
don.syn.1500 <- clean(don.syn.1500)
don.nonsyn.1500 <- clean(don.nonsyn.1500)
don.del.1500 <- clean(don.del.1500)

# Put the trait-intensity samples into a data frame for plotting. We will make
# side-by-side boxplots for prop. var. explained for each partition, separated
# by trait and resampling intensity.

yield.250 <- data.frame(
    Value=c(
        yield.null.250,
        yield.nonc.250,
        yield.syn.250,
        yield.nonsyn.250,
        yield.del.250),
    Label=c(
        rep("N", length(yield.null.250)),
        rep("NC", length(yield.nonc.250)),
        rep("S", length(yield.syn.250)),
        rep("NS", length(yield.nonsyn.250)),
        rep("D", length(yield.del.250)))
    )

yield.500 <- data.frame(
    Value=c(
        yield.null.500,
        yield.nonc.500,
        yield.syn.500,
        yield.nonsyn.500,
        yield.del.500),
    Label=c(
        rep("N", length(yield.null.500)),
        rep("NC", length(yield.nonc.500)),
        rep("S", length(yield.syn.500)),
        rep("NS", length(yield.nonsyn.500)),
        rep("D", length(yield.del.500)))
    )

yield.750 <- data.frame(
    Value=c(
        yield.null.750,
        yield.nonc.750,
        yield.syn.750,
        yield.nonsyn.750,
        yield.del.750),
    Label=c(
        rep("N", length(yield.null.750)),
        rep("NC", length(yield.nonc.750)),
        rep("S", length(yield.syn.750)),
        rep("NS", length(yield.nonsyn.750)),
        rep("D", length(yield.del.750)))
    )

yield.1000 <- data.frame(
    Value=c(
        yield.null.1000,
        yield.nonc.1000,
        yield.syn.1000,
        yield.nonsyn.1000,
        yield.del.1000),
    Label=c(
        rep("N", length(yield.null.1000)),
        rep("NC", length(yield.nonc.1000)),
        rep("S", length(yield.syn.1000)),
        rep("NS", length(yield.nonsyn.1000)),
        rep("D", length(yield.del.1000)))
    )

yield.1250 <- data.frame(
    Value=c(
        yield.null.1250,
        yield.nonc.1250,
        yield.syn.1250,
        yield.nonsyn.1250,
        yield.del.1250),
    Label=c(
        rep("N", length(yield.null.1250)),
        rep("NC", length(yield.nonc.1250)),
        rep("S", length(yield.syn.1250)),
        rep("NS", length(yield.nonsyn.1250)),
        rep("D", length(yield.del.1250)))
    )

yield.1500 <- data.frame(
    Value=c(
        yield.null.1500,
        yield.nonc.1500,
        yield.syn.1500,
        yield.nonsyn.1500,
        yield.del.1500),
    Label=c(
        rep("N", length(yield.null.1500)),
        rep("NC", length(yield.nonc.1500)),
        rep("S", length(yield.syn.1500)),
        rep("NS", length(yield.nonsyn.1500)),
        rep("D", length(yield.del.1500)))
    )

don.250 <- data.frame(
    Value=c(
        don.null.250,
        don.nonc.250,
        don.syn.250,
        don.nonsyn.250,
        don.del.250),
    Label=c(
        rep("N", length(don.null.250)),
        rep("NC", length(don.nonc.250)),
        rep("S", length(don.syn.250)),
        rep("NS", length(don.nonsyn.250)),
        rep("D", length(don.del.250)))
    )

don.500 <- data.frame(
    Value=c(
        don.null.500,
        don.nonc.500,
        don.syn.500,
        don.nonsyn.500,
        don.del.500),
    Label=c(
        rep("N", length(don.null.500)),
        rep("NC", length(don.nonc.500)),
        rep("S", length(don.syn.500)),
        rep("NS", length(don.nonsyn.500)),
        rep("D", length(don.del.500)))
    )

don.750 <- data.frame(
    Value=c(
        don.null.750,
        don.nonc.750,
        don.syn.750,
        don.nonsyn.750,
        don.del.750),
    Label=c(
        rep("N", length(don.null.750)),
        rep("NC", length(don.nonc.750)),
        rep("S", length(don.syn.750)),
        rep("NS", length(don.nonsyn.750)),
        rep("D", length(don.del.750)))
    )

don.1000 <- data.frame(
    Value=c(
        don.null.1000,
        don.nonc.1000,
        don.syn.1000,
        don.nonsyn.1000,
        don.del.1000),
    Label=c(
        rep("N", length(don.null.1000)),
        rep("NC", length(don.nonc.1000)),
        rep("S", length(don.syn.1000)),
        rep("NS", length(don.nonsyn.1000)),
        rep("D", length(don.del.1000)))
    )

don.1250 <- data.frame(
    Value=c(
        don.null.1250,
        don.nonc.1250,
        don.syn.1250,
        don.nonsyn.1250,
        don.del.1250),
    Label=c(
        rep("N", length(don.null.1250)),
        rep("NC", length(don.nonc.1250)),
        rep("S", length(don.syn.1250)),
        rep("NS", length(don.nonsyn.1250)),
        rep("D", length(don.del.1250)))
    )

don.1500 <- data.frame(
    Value=c(
        don.null.1500,
        don.nonc.1500,
        don.syn.1500,
        don.nonsyn.1500,
        don.del.1500),
    Label=c(
        rep("N", length(don.null.1500)),
        rep("NC", length(don.nonc.1500)),
        rep("S", length(don.syn.1500)),
        rep("NS", length(don.nonsyn.1500)),
        rep("D", length(don.del.1500)))
    )

yield.250$Label <- factor(yield.250$Label, levels=c("N", "NC", "S", "NS", "D"))
yield.500$Label <- factor(yield.500$Label, levels=c("N", "NC", "S", "NS", "D"))
yield.750$Label <- factor(yield.750$Label, levels=c("N", "NC", "S", "NS", "D"))
yield.1000$Label <- factor(yield.1000$Label, levels=c("N", "NC", "S", "NS", "D"))
yield.1250$Label <- factor(yield.1250$Label, levels=c("N", "NC", "S", "NS", "D"))
yield.1500$Label <- factor(yield.1500$Label, levels=c("N", "NC", "S", "NS", "D"))
don.250$Label <- factor(don.250$Label, levels=c("N", "NC", "S", "NS", "D"))
don.500$Label <- factor(don.500$Label, levels=c("N", "NC", "S", "NS", "D"))
don.750$Label <- factor(don.750$Label, levels=c("N", "NC", "S", "NS", "D"))
don.1000$Label <- factor(don.1000$Label, levels=c("N", "NC", "S", "NS", "D"))
don.1250$Label <- factor(don.1250$Label, levels=c("N", "NC", "S", "NS", "D"))
don.1500$Label <- factor(don.1500$Label, levels=c("N", "NC", "S", "NS", "D"))

# Make plots
pdf(file="GCTA_Yield.pdf", width=11, height=8.5)
par(mfrow=c(2, 3))
boxplot(yield.250$Value ~ yield.250$Label, main="250", xlab="Partition", ylab="VG/VP", ylim=c(0, 0.15))
boxplot(yield.500$Value ~ yield.500$Label, main="Yield\n500", xlab="Partition", ylab="VG/VP", ylim=c(0, 0.15))
boxplot(yield.750$Value ~ yield.750$Label, main="750", xlab="Partition", ylab="VG/VP", ylim=c(0, 0.15))
boxplot(yield.1000$Value ~ yield.1000$Label, main="1000", xlab="Partition", ylab="VG/VP", ylim=c(0, 0.15))
boxplot(yield.1250$Value ~ yield.1250$Label, main="1250", xlab="Partition", ylab="VG/VP", ylim=c(0, 0.15))
boxplot(yield.1500$Value ~ yield.1500$Label, main="1500", xlab="Partition", ylab="VG/VP", ylim=c(0, 0.15))
dev.off()

pdf(file="GCTA_DON.pdf", width=11, height=8.5)
par(mfrow=c(2, 3))
boxplot(don.250$Value ~ don.250$Label, main="250", xlab="Partition", ylab="VG/VP", ylim=c(0, 0.15))
boxplot(don.500$Value ~ don.500$Label, main="DON\n500", xlab="Partition", ylab="VG/VP", ylim=c(0, 0.15))
boxplot(don.750$Value ~ don.750$Label, main="750", xlab="Partition", ylab="VG/VP", ylim=c(0, 0.15))
boxplot(don.1000$Value ~ don.1000$Label, main="1000", xlab="Partition", ylab="VG/VP", ylim=c(0, 0.15))
boxplot(don.1250$Value ~ don.1250$Label, main="1250", xlab="Partition", ylab="VG/VP", ylim=c(0, 0.15))
boxplot(don.1500$Value ~ don.1500$Label, main="1500", xlab="Partition", ylab="VG/VP", ylim=c(0, 0.15))
dev.off()
