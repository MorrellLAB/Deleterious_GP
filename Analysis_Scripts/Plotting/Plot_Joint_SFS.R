#   Generate and plot a joint site frequency spectrum, given two files that
#   have vectors of derived allele frequencies. Takes three arguments:
#       1) File with frequencies in partition 1
#       2) File with frequencies in partition 2
#       3) output filename

library(fields)
#   Take arguments. 
args <- commandArgs(TRUE)

part1 <- args[1]
part2 <- args[2]
output <- args[3]

#   Define a function that returns the matrix cell assignment for a single SNP
#   based on its frequencies in partition 1 and partition 2.
jsfs_cell <- function(freqs, binsize, folded=TRUE) {
    #   Set the maximum frequency, if folded=TRUE, it is 1.0. Else 0.5.
    if(folded) {
        maxfreq <- 0.5
    } else {
        maxfreq <- 1.0
    }
    #   Generate the bins
    bins <- seq(0.0, maxfreq, by=binsize)
    #   Generate the names of the frequency classes for matching
    classes <- names(table(cut(bins, breaks=bins, include.lowest=TRUE)))
    #   Then, get which bins the two frequencies fall in
    cell <- cut(freqs, breaks=bins, include.lowest=TRUE)
    cell <- match(cell, classes)
    return(cell)
}

#   Define a function that generates a matrix object and fills it with data
#   from the joint SFS assignments
jsfs_matrix <- function(cells, binsize, folded=TRUE) {
    #   Set the maximum frequency based on folded
    if(folded) {
        maxfreq <- 0.5
    } else {
        maxfreq <- 1.0
    }
    #   Generate the bins
    bins <- seq(0.0, maxfreq, by=binsize)
    #   Build the matrix
    #       Get the names of the frequency classes. This is somewhat ugly, but
    #       it works.
    freq.names <- names(table(cut(bins, breaks=bins, include.lowest=TRUE)))
    #       Build an empty matrix to populate
    freq.matrix <- matrix(0, length(freq.names)**2, nrow=length(freq.names), ncol=length(freq.names))
    #       Then populate it
    for(i in 1:ncol(cells)) {
            #   Increment the correct cell.
            mat_row <- cells[2, i]
            mat_col <- cells[1, i]
            freq.matrix[mat_row, mat_col] <- freq.matrix[mat_row, mat_col] + 1
        }
    #   Name the dimensions
    row.names(freq.matrix) <- freq.names
    colnames(freq.matrix) <- freq.names
    #   Convert to proportions
    #freq.matrix <- freq.matrix/sum(freq.matrix)
    #   Make a list with the components for image()
    hmp <- list(
        x=(bins[1:length(bins)-1] + bins[2: length(bins)])/2,
        y=(bins[1:length(bins)-1] + bins[2: length(bins)])/2,
        z=log(freq.matrix+1)
        )
    return(hmp)
}

#   Define the heatmap colors. This is from the 9-level YlOrRd pallette from
#   RColorBrewer
heatmap_colors <- c(
    "#fff7ec",
    "#fee8c8",
    "#fdd49e",
    "#fdbb84",
    "#fc8d59",
    "#ef6548",
    "#d7301f",
    "#b30000",
    "#7f0000")

#   Read in the frequencies and make them numeric vectors
part1.freq <- read.table(part1, header=F)
part1.freq <- as.numeric(part1.freq$V1)
part2.freq <- read.table(part2, header=F)
part2.freq <- as.numeric(part2.freq$V1)

#   Get the cell assignments for each SNP
joint_sfs <- apply(
    cbind(part1.freq, part2.freq),
    1,
    jsfs_cell,
    binsize=0.05,
    folded=FALSE
    )

#   Then build the joint SFS
joint_sfs_mat <- jsfs_matrix(joint_sfs, binsize=0.05, folded=FALSE)
print(joint_sfs_mat)
#   And plot it!
#   Then we start building the image
#       This first call builds the base of the heatmap, with axis labels
#       and title string
brk <- c(0, 0.0001, 0.00025, 0.0005, 0.001, 0.0025, 0.005, 0.01, 0.1, 0.25, 0.3, 0.4, 0.6)
pdf(file=output, 6, 6)
image.plot(joint_sfs_mat,
    col=rev(heat.colors(12)),
    breaks=c(seq(0, 5, by=0.5), 6, 7),
    lab.breaks=c(seq(0, 5, by=0.5), 6, 7),
    xlab="Cycle 3 Derived Allele Frequency",
    ylab="Cycle 1 Derived Allele Frequency",
    main="Random Panel",
    cex.axis=1.2)
#   Then we add our own axes
#       x-axis. number of rows in our matrix corresponds to the read length
#axis(1, at=seq(0, 1, length.out=nrow(joint_sfs_mat$z)), labels=row.names(joint_sfs_mat$z), cex.axis=1.5)
#       y-axis. We use the number of quality scores and the character vector
#       of the scores themselves to build the axis
#axis(2, at=seq(0, 1, length.out=ncol(joint_sfs_mat$z)), labels=colnames(joint_sfs_mat$z), cex.axis=1.5)
#   Put a box on it
box()
#   Draw the diagonal
abline(0, 1)
dev.off()
