#   Calculate the number of exome capture targets per Mb in the pseudomolecule
#   assembly. Calculates in non-overlapping windows. Counts partial overlaps of
#   capture targets with the window as a full target (may double-count a few).

#   Define a function to return the number of overlapping intervals, given an
#   interval size.
num_targets <- function(chr_len, win_size, targets) {
    #   Get the break points of the windows
    win_points <- seq(0, chr_len, by=win_size)
    #   If the window size divides evenly into the chromosome length, we are
    #   fine. Otherwise, put the length on as the final endpoint.
    if(chr_len %% win_size != 0) {
        win_points <- c(win_points, chr_len)
    }
    #   Then, count how many exome capture targets are in each window
    num_in_win <- sapply(
        2:length(win_points),
        function(x) {
            in_win <- (targets$V2 >= win_points[x-1]) & (targets$V3 <= win_points[x])
            return(sum(in_win))
        }
    )
    return(num_in_win)
}

#   Function to transform a vector of numbers into a list of intervals
make_intervals <- function(endpoints) {
    starts <- endpoints[1:length(endpoints)-1]
    ends <- endpoints[2:length(endpoints)]
    return(data.frame(Start=starts, End=ends))
}

setwd("/Volumes/Data_Disk/Dropbox/Projects/DM_GenomicPrediction/Results/Annotations/")
capture <- read.table("50x_Capture.bed", header=F)
chr_sizes <- read.table("Chrom_Lengths.txt", header=T)

#   Get the chromosome names
chromosomes <- as.character(unique(capture$V1))

#   Get the number of targets in each window on each chromosome
target_num <- sapply(
    chromosomes,
    function(x) {
        chr_len <- as.integer(chr_sizes[chr_sizes$Chromosome == x, "Length"])
        chr_targets <- capture[capture$V1 == x,]
        return(num_targets(chr_len, 1000000, chr_targets))
    }
    )

chr_windows <- sapply(
    chromosomes,
    function(c) {
        len <- as.integer(chr_sizes[chr_sizes$Chromosome == c, "Length"])
        win_points <- seq(0, len, by=1000000)
        if(len %% 1000000 != 0) {
            win_points <- c(win_points, len)
        }
        wins <- make_intervals(win_points)
        return(wins)
    })

num_target_df <- data.frame(
    Chromosome=c(
        rep("chr1H", length(chr_windows["Start", "chr1H"][[1]])),
        rep("chr2H", length(chr_windows["Start", "chr2H"][[1]])),
        rep("chr3H", length(chr_windows["Start", "chr3H"][[1]])),
        rep("chr4H", length(chr_windows["Start", "chr4H"][[1]])),
        rep("chr5H", length(chr_windows["Start", "chr5H"][[1]])),
        rep("chr6H", length(chr_windows["Start", "chr6H"][[1]])),
        rep("chr7H", length(chr_windows["Start", "chr7H"][[1]])),
        rep("chrUn", length(chr_windows["Start", "chrUn"][[1]]))
            ),
    Start=as.vector(unlist(chr_windows["Start",])),
    End=as.vector(unlist(chr_windows["End", ])),
    NExCap=as.vector(unlist(target_num))
    )

write.table(
    num_target_df,
    file="ExomeCaptureTargets_per_Mb.txt",
    sep="\t",
    quote=F,
    row.names=F)
