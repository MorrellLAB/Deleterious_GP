#   Generate and plot a joint site frequency spectrum, given two files that
#   have vectors of derived allele frequencies.

library(fields)

freqs <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Genotype_Freqs/DAF_By_Cycle.txt.gz", header=TRUE)

nonc <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Noncoding.txt.gz", header=FALSE)$V1)
syn <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Synonymous.txt.gz", header=FALSE)$V1)
nonsyn <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Nonsynonymous.txt.gz", header=FALSE)$V1)
del <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Deleterious.txt.gz", header=FALSE)$V1)

nc.freqs <- freqs[freqs$SNP_ID %in% nonc,]
syn.freqs <- freqs[freqs$SNP_ID %in% syn,]
nonsyn.freqs <- freqs[freqs$SNP_ID %in% nonsyn,]
del.freqs <- freqs[freqs$SNP_ID %in% del,]

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


#   Get the cell assignments for each SNP
nc.joint_sfs <- apply(nc.freqs[,c("C1_DAF", "C2_DAF")], 1, jsfs_cell, binsize=0.05, folded=FALSE)
syn.joint_sfs <- apply(syn.freqs[,c("C1_DAF", "C2_DAF")], 1, jsfs_cell, binsize=0.05, folded=FALSE)
nonsyn.joint_sfs <- apply(nonsyn.freqs[,c("C1_DAF", "C2_DAF")], 1, jsfs_cell, binsize=0.05, folded=FALSE)
del.joint_sfs <- apply(del.freqs[,c("C1_DAF", "C2_DAF")], 1, jsfs_cell, binsize=0.05, folded=FALSE)

#   Then build the joint SFS
nc.joint_sfs_mat <- jsfs_matrix(nc.joint_sfs, binsize=0.05, folded=FALSE)
syn.joint_sfs_mat <- jsfs_matrix(syn.joint_sfs, binsize=0.05, folded=FALSE)
nonsyn.joint_sfs_mat <- jsfs_matrix(nonsyn.joint_sfs, binsize=0.05, folded=FALSE)
del.joint_sfs_mat <- jsfs_matrix(del.joint_sfs, binsize=0.05, folded=FALSE)

pdf(file="Noncoding_Joint_SFS.pdf", 6, 6)
image.plot(nc.joint_sfs_mat,
    col=rev(heat.colors(12)),
    breaks=c(seq(0, 5, by=0.5), 7.5, 10),
    lab.breaks=c(seq(0, 5, by=0.5), 7.5, 10),
    xlab="Cycle 3 Derived Allele Frequency",
    ylab="Cycle 1 Derived Allele Frequency",
    main="Noncoding",
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

pdf(file="Synonymous_Joint_SFS.pdf", 6, 6)
image.plot(syn.joint_sfs_mat,
    col=rev(heat.colors(12)),
    breaks=c(seq(0, 5, by=0.5), 7.5, 10),
    lab.breaks=c(seq(0, 5, by=0.5), 7.5, 10),
    xlab="Cycle 3 Derived Allele Frequency",
    ylab="Cycle 1 Derived Allele Frequency",
    main="Synonymous",
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

pdf(file="Nonsynonymous_Joint_SFS.pdf", 6, 6)
image.plot(nonsyn.joint_sfs_mat,
    col=rev(heat.colors(12)),
    breaks=c(seq(0, 5, by=0.5), 7.5, 10),
    lab.breaks=c(seq(0, 5, by=0.5), 7.5, 10),
    xlab="Cycle 3 Derived Allele Frequency",
    ylab="Cycle 1 Derived Allele Frequency",
    main="Nonsynonymous",
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

pdf(file="Deleterious_Joint_SFS.pdf", 6, 6)
image.plot(del.joint_sfs_mat,
    col=rev(heat.colors(12)),
    breaks=c(seq(0, 5, by=0.5), 7.5, 10),
    lab.breaks=c(seq(0, 5, by=0.5), 7.5, 10),
    xlab="Cycle 3 Derived Allele Frequency",
    ylab="Cycle 1 Derived Allele Frequency",
    main="Deleterious",
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
