# Summarize the pairwise diversity among separate cycles at the BOPA markers

# Read all the chromosomes in
chr1H <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/AlphaPeel_Expanded/chr1H/AlphaPeel_BOPA_chr1H.txt", header=FALSE)
chr2H <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/AlphaPeel_Expanded/chr2H/AlphaPeel_BOPA_chr2H.txt", header=FALSE)
chr3H <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/AlphaPeel_Expanded/chr3H/AlphaPeel_BOPA_chr3H.txt", header=FALSE)
chr4H <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/AlphaPeel_Expanded/chr4H/AlphaPeel_BOPA_chr4H.txt", header=FALSE)
chr5H <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/AlphaPeel_Expanded/chr5H/AlphaPeel_BOPA_chr5H.txt", header=FALSE)
chr6H <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/AlphaPeel_Expanded/chr6H/AlphaPeel_BOPA_chr6H.txt", header=FALSE)
chr7H <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/AlphaPeel_Expanded/chr7H/AlphaPeel_BOPA_chr7H.txt", header=FALSE)

# Read the summary of the segregating sites
seg <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Imputation/AlphaPeel/BOPA_Segregation_Summary.txt", header=TRUE)

# Merge chromosomes by the first column
dat <- as.matrix(cbind(
    chr1H,
    chr2H[,-1],
    chr3H[,-1],
    chr4H[,-1],
    chr5H[,-1],
    chr6H[,-1],
    chr7H[,-1]))

colnames(dat) <- NULL

# Summarize missing data
prop_missing <- function(x) {
    pmiss <- sum(x == 9)/length(x)
    return(pmiss)
}

# Define a function that will calculate pairwise similarity
pairwise_div <- function(x) {
    # Missing data is coded as 9
    missing <- c(9)
    calls <- x[!is.element(x, missing)]
    # Count the alleles
    a1 <- (sum(calls == 0) * 2) + sum(calls == 1)
    a2 <- (sum(calls == 2) * 2) + sum(calls == 1)
    # pairwise similarity is
    #   [(A1Count choose 2) + (A2Count choose 2)] / (NChr choose 2)
    # Basically, the number of A1-A1 comparisons + A2-A2 comparisons, divided
    # by the total number of comparisons
    sim <- (choose(a1, 2) + choose(a2, 2)) / choose(a1+a2, 2)
    # Diversity is 1-similarity
    return(1-sim)
    }

# Define a function that will calculate observed heterozygosity
obs_het <- function(x) {
    missing <- c(9)
    calls <- x[!is.element(x, missing)]
    het <- sum(calls == 1)/length(calls)
    return(het)
}

# Partition the sample by cycle
c0 <- 1:21
c1 <- grep("MS10", dat[,1])
c2 <- grep("MS11", dat[,1])
c3 <- grep("MS12", dat[,1])

# Do the same partitions, but in the segregating summary. The summary table
# gives them by individual, and we want by family. Repeating the numbers will
# mess up the variance estimates, so we just take the first individual from each
# family.
c1.seg <- grep("MS10.+-001", seg$ID, perl=TRUE)
c2.seg <- grep("MS11.+-001", seg$ID, perl=TRUE)
c3.seg <- grep("MS12.+-001", seg$ID, perl=TRUE)


# Make a data frame to store these in
result <- data.frame(
    Cycle=c("C0", "C1", "C2", "C3"),
    NInd=c(length(c0), length(c1), length(c2), length(c3)),
    NFam=c(NA, length(c1.seg), length(c2.seg), length(c3.seg)),
    AvgPropMissing=c(
        mean(apply(dat[c0, -1], 1, prop_missing)),
        mean(apply(dat[c1, -1], 1, prop_missing)),
        mean(apply(dat[c2, -1], 1, prop_missing)),
        mean(apply(dat[c3, -1], 1, prop_missing))),
    AvgPairwiseDiv=c(
        mean(apply(dat[c0, -1], 2, pairwise_div)),
        mean(apply(dat[c1, -1], 2, pairwise_div)),
        mean(apply(dat[c2, -1], 2, pairwise_div)),
        mean(apply(dat[c3, -1], 2, pairwise_div))),
    AvgObsHet=c(
        mean(apply(dat[c0, -1], 1, obs_het)),
        mean(apply(dat[c1, -1], 1, obs_het)),
        mean(apply(dat[c2, -1], 1, obs_het)),
        mean(apply(dat[c3, -1], 1, obs_het))),
    AvgSegFam=c(
        NA,
        mean(seg[c1.seg, "NSeg"], na.rm=TRUE),
        mean(seg[c2.seg, "NSeg"], na.rm=TRUE),
        mean(seg[c3.seg, "NSeg"], na.rm=TRUE)),
    SdSegFam=c(
        NA,
        sd(seg[c1.seg, "NSeg"], na.rm=TRUE),
        sd(seg[c2.seg, "NSeg"], na.rm=TRUE),
        sd(seg[c3.seg, "NSeg"], na.rm=TRUE))
    )

# Write the result
write.csv(
    result,
    file="GenomicPrediction_BOPA_Summary.csv",
    row.names=FALSE,
    quote=FALSE)

# Plot the number of segregating markers per family
pdf(file="Seg_Markers_Per_Fam.pdf", 6, 6)
toplot <- data.frame(
    Cycle=c(
        rep("C1", length(c1.seg)),
        rep("C2", length(c2.seg)),
        rep("C3", length(c3.seg))),
    Value=c(
        seg[c1.seg, "NSeg"],
        seg[c2.seg, "NSeg"],
        seg[c3.seg, "NSeg"]))
boxplot(
    toplot$Value ~ toplot$Cycle,
    xlab="Cycle",
    ylab="Number of Polymorphic Markers Per Family",
    lwd=2)
dev.off()
