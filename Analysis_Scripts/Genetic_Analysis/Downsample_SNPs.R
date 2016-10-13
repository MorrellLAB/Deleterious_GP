#   Script to downsample SNPs from one MAF distribution to match another. This
#   was written generate samples from two classes of SNPs that have very
#   different MAF distributions that have the same MAF distribution. This is
#   accomplished in a somewhat simplistic way:
#       1) The source distribution (the one that is used as the "first") is
#          cut into bins, with the count of the SNPs in each bin representing
#          the probability it is chosen.
#       2) For each bin in the target distribution (the one to be downsampled),
#          the number of SNPs from the previous sampling step are randomly
#          sampled.
#       3) The IDs of the subsampled SNPs are returned.

#   Required input files:
#       1) The source minor allele frequencies (from PLINK)
#       2) The target minor allele frequencies (from PLINK)

source_freq_dist <- function(maf, binsize) {
    #   Generate the breakpoints
    bins <- seq(0, 0.5, by=binsize)
    #   Count the number of SNPs in each bin
    counts <- table(cut(maf, breaks=bins, include.lowest=TRUE))
    return(counts)
}

downsample <- function(maf, binsize, weights) {
    #   To downsample, we will decide which bin the target SNP lands in, then
    #   assign the weights of the source distribution to it.
    bins <- seq(0, 0.5, by=binsize)
    #   Build a data frame that gives the frequency bin points
    intervals <- data.frame(
        Start=bins[1:length(bins)-1],
        End=bins[2:length(bins)]
        )
    #   Then, for each chosen interval, sample the SNPs in that frequency class
    sampled <- unlist(sapply(
        seq_along(weights),
        function(x) {
            interval_start <- intervals[x, "Start"]
            interval_end <- intervals[x, "End"]
            nsamp <- weights[x]
            #   Get all the SNPs that are within the frequency class
            in_int <- as.character(maf[maf$MAF > interval_start & maf$MAF <= interval_end, "SNP"])
            #   Randomly sample the IDs
            return(sample(in_int, nsamp, replace=F))
        }))
    #   Then, randomly prune the set down to the subsampled size
    sampled <- sample(sampled, 500, replace=FALSE)
    return(sampled)
}

args <- commandArgs(TRUE)
source_maf <- read.table(args[1], header=T)
target_maf <- read.table(args[2], header=T)

source.dist <- source_freq_dist(source_maf$MAF, 0.01)
target.dist <- downsample(target_maf, 0.01, source.dist)

#   Density plots to check that the distributions are equal
#plot(density(source_maf[sample(1:nrow(source_maf)), "MAF"]))
#plot(density(target_maf[target_maf$SNP %in% target.dist, "MAF"]))

#   Print the SNPs to stdout
write(target.dist, file="", ncolumns=1)
