#   A script to print out summary information for all the families used in the
#   genomic prediction experiment. Writes a tab-delimited file with the
#   following information:
#       Cycle
#       Family ID
#       Number of Markers
#       N. 1H
#       N. 2H
#       N. 3H
#       N. 4H
#       N. 5H
#       N. 6H
#       N. 7H
#       Mean XO per gen
#       Mean XO 1H per gen
#       Mean XO 2H per gen
#       Mean XO 3H per gen
#       Mean XO 4H per gen
#       Mean XO 5H per gen
#       Mean XO 6H per gen
#       Mean XO 7H per gen

#   We use r/QTL for this job
library(qtl)

#   Define a function to extract family IDs
get_famid <- function(path) {
    without_dir <- unlist(strsplit(path, '/'))[2]
    without_ext <- unlist(strsplit(without_dir, '[.]'))[1]
    return(without_ext)
}

#   Define a function for counting the number of markers on each chromosome
get_chr_nmarkers <- function(cross) {
    #   These are the chromosomes in barley
    chroms <- c("1H", "2H", "3H", "4H", "5H", "6H", "7H")
    #   Then count the number of markers on each chromosome. In the case of
    #   a chromosome have 0 markers, we return 0 (since ncol returns NULL)
    counts  <- sapply(
        chroms,
        function(x) {
            n <- ncol(cross$geno[[x]]$data)
            if(is.null(n)) {
                n <- 0
            }
            return(n)
        }
        )
    return(counts)
}

#   Define a function to calculate the mean number of crossovers per chromosome
get_chr_xo <- function(cross) {
    chroms <- c("1H", "2H", "3H", "4H", "5H", "6H", "7H")
    counts <- sapply(
        chroms,
        function(x) {
            #   If the chromosome is not in the cross object, then we cannot
            #   estimate the crossover rate. We return NA
            if(!x %in% chrnames(cross)) {
                return(NA)
            }
            else {
                return(mean(countXO(cross, chr=x)))
            }
        }
    )
    return(counts)
}

#   Define a function to calculate the average marker spacing over all
#   chromosomes
get_marker_spacing <- function(cross) {
    #   Define a sub-function to calculate the average interval
    intervals <- function(distances) {
        sp <- sapply(
            seq_along(distances)[1:length(distances)-1],
            function(d) {
                return(distances[d+1] - distances[d])
            }
        )
    }
    chroms <- c("1H", "2H", "3H", "4H", "5H", "6H", "7H")
    spacing <- sapply(
        chroms,
        function(x) {
            #   If the chromosome is not in the cross object, return NA
            if(!x %in% chrnames(cross)) {
                return(NA)
            }
            #   Else if the number of markers on the chromosome is 1, then
            #   return NA
            else if(length(cross$geno[[x]]$map) == 1) {
                return(NA)
            }
            #   otherwise calculate the intervals
            else {
                return(mean(intervals(cross$geno[[x]]$map)))
            }
        }
    )
    return(spacing)
}

#   Read in all the crosses
c1_fnames <- list.files("Cycle_1", pattern="*.csv", full.names=TRUE)
c2_fnames <- list.files("Cycle_2", pattern="*.csv", full.names=TRUE)
c3_fnames <- list.files("Cycle_3", pattern="*.csv", full.names=TRUE)

c1_crosses <- lapply(
    seq_along(c1_fnames),
    function(x) {
        read.cross(format="csv", file=c1_fnames[[x]])
        }
    )
c2_crosses <- lapply(
    seq_along(c2_fnames),
    function(x) {
        read.cross(format="csv", file=c2_fnames[[x]])
        }
    )
c3_crosses <- lapply(
    seq_along(c3_fnames),
    function(x) {
        read.cross(format="csv", file=c3_fnames[[x]])
        }
    )

#   For writing the data later, we generate a vector of family names
c1_famids <- sapply(c1_fnames, get_famid)
c2_famids <- sapply(c2_fnames, get_famid)
c3_famids <- sapply(c3_fnames, get_famid)

#   For each family, get the number of markers
c1_nmarkers <- unlist(lapply(c1_crosses, totmar))
c2_nmarkers <- unlist(lapply(c2_crosses, totmar))
c3_nmarkers <- unlist(lapply(c3_crosses, totmar))

#   Then get the chromosome-specific counts
c1_chr_nmarkers <- lapply(c1_crosses, get_chr_nmarkers)
c2_chr_nmarkers <- lapply(c2_crosses, get_chr_nmarkers)
c3_chr_nmarkers <- lapply(c3_crosses, get_chr_nmarkers)

#   Then calculate the mean number of crossovers on all chromosomes
c1_tot_xo <- unlist(lapply(c1_crosses, function(x) mean(countXO(x))))
c2_tot_xo <- unlist(lapply(c2_crosses, function(x) mean(countXO(x))))
c3_tot_xo <- unlist(lapply(c3_crosses, function(x) mean(countXO(x))))

#   And chromosome-by-chromosome
c1_chr_xo <- lapply(c1_crosses, get_chr_xo)
c2_chr_xo <- lapply(c2_crosses, get_chr_xo)
c3_chr_xo <- lapply(c3_crosses, get_chr_xo)

#   Then calculate the marker spacing
c1_marker_spacing <- lapply(c1_crosses, get_marker_spacing)
c2_marker_spacing <- lapply(c2_crosses, get_marker_spacing)
c3_marker_spacing <- lapply(c3_crosses, get_marker_spacing)

#   Finally, put it all together into a data frame
pop_summary <- data.frame(
    Cycle=c(rep(1, length(c1_fnames)), rep(2, length(c2_fnames)), rep(3, length(c3_fnames))),
    FamilyID=c(c1_famids, c2_famids, c3_famids),
    NMarkers=c(c1_nmarkers, c2_nmarkers, c3_nmarkers),
    NMarkers_1H=c(
        unlist(lapply(c1_chr_nmarkers, function(x) return(x[["1H"]]))),
        unlist(lapply(c2_chr_nmarkers, function(x) return(x[["1H"]]))),
        unlist(lapply(c3_chr_nmarkers, function(x) return(x[["1H"]])))),
    NMarkers_2H=c(
        unlist(lapply(c1_chr_nmarkers, function(x) return(x[["2H"]]))),
        unlist(lapply(c2_chr_nmarkers, function(x) return(x[["2H"]]))),
        unlist(lapply(c3_chr_nmarkers, function(x) return(x[["2H"]])))),
    NMarkers_3H=c(
        unlist(lapply(c1_chr_nmarkers, function(x) return(x[["3H"]]))),
        unlist(lapply(c2_chr_nmarkers, function(x) return(x[["3H"]]))),
        unlist(lapply(c3_chr_nmarkers, function(x) return(x[["3H"]])))),
    NMarkers_4H=c(
        unlist(lapply(c1_chr_nmarkers, function(x) return(x[["4H"]]))),
        unlist(lapply(c2_chr_nmarkers, function(x) return(x[["4H"]]))),
        unlist(lapply(c3_chr_nmarkers, function(x) return(x[["4H"]])))),
    NMarkers_5H=c(
        unlist(lapply(c1_chr_nmarkers, function(x) return(x[["5H"]]))),
        unlist(lapply(c2_chr_nmarkers, function(x) return(x[["5H"]]))),
        unlist(lapply(c3_chr_nmarkers, function(x) return(x[["5H"]])))),
    NMarkers_6H=c(
        unlist(lapply(c1_chr_nmarkers, function(x) return(x[["6H"]]))),
        unlist(lapply(c2_chr_nmarkers, function(x) return(x[["6H"]]))),
        unlist(lapply(c3_chr_nmarkers, function(x) return(x[["6H"]])))),
    NMarkers_7H=c(
        unlist(lapply(c1_chr_nmarkers, function(x) return(x[["7H"]]))),
        unlist(lapply(c2_chr_nmarkers, function(x) return(x[["7H"]]))),
        unlist(lapply(c3_chr_nmarkers, function(x) return(x[["7H"]])))),
    MarkerSpacing_1H=c(
        unlist(lapply(c1_marker_spacing, function(x) return(x[["1H"]]))),
        unlist(lapply(c2_marker_spacing, function(x) return(x[["1H"]]))),
        unlist(lapply(c3_marker_spacing, function(x) return(x[["1H"]])))),
    MarkerSpacing_2H=c(
        unlist(lapply(c1_marker_spacing, function(x) return(x[["2H"]]))),
        unlist(lapply(c2_marker_spacing, function(x) return(x[["2H"]]))),
        unlist(lapply(c3_marker_spacing, function(x) return(x[["2H"]])))),
    MarkerSpacing_3H=c(
        unlist(lapply(c1_marker_spacing, function(x) return(x[["3H"]]))),
        unlist(lapply(c2_marker_spacing, function(x) return(x[["3H"]]))),
        unlist(lapply(c3_marker_spacing, function(x) return(x[["3H"]])))),
    MarkerSpacing_4H=c(
        unlist(lapply(c1_marker_spacing, function(x) return(x[["4H"]]))),
        unlist(lapply(c2_marker_spacing, function(x) return(x[["4H"]]))),
        unlist(lapply(c3_marker_spacing, function(x) return(x[["4H"]])))),
    MarkerSpacing_5H=c(
        unlist(lapply(c1_marker_spacing, function(x) return(x[["5H"]]))),
        unlist(lapply(c2_marker_spacing, function(x) return(x[["5H"]]))),
        unlist(lapply(c3_marker_spacing, function(x) return(x[["5H"]])))),
    MarkerSpacing_6H=c(
        unlist(lapply(c1_marker_spacing, function(x) return(x[["6H"]]))),
        unlist(lapply(c2_marker_spacing, function(x) return(x[["6H"]]))),
        unlist(lapply(c3_marker_spacing, function(x) return(x[["6H"]])))),
    MarkerSpacing_7H=c(
        unlist(lapply(c1_marker_spacing, function(x) return(x[["7H"]]))),
        unlist(lapply(c2_marker_spacing, function(x) return(x[["7H"]]))),
        unlist(lapply(c3_marker_spacing, function(x) return(x[["7H"]])))),
    MeanXO=c(c1_tot_xo, c2_tot_xo, c3_tot_xo),
    MeanXO_1H=c(
        unlist(lapply(c1_chr_xo, function(x) return(x[["1H"]]))),
        unlist(lapply(c2_chr_xo, function(x) return(x[["1H"]]))),
        unlist(lapply(c3_chr_xo, function(x) return(x[["1H"]])))),
    MeanXO_2H=c(
        unlist(lapply(c1_chr_xo, function(x) return(x[["2H"]]))),
        unlist(lapply(c2_chr_xo, function(x) return(x[["2H"]]))),
        unlist(lapply(c3_chr_xo, function(x) return(x[["2H"]])))),
    MeanXO_3H=c(
        unlist(lapply(c1_chr_xo, function(x) return(x[["3H"]]))),
        unlist(lapply(c2_chr_xo, function(x) return(x[["3H"]]))),
        unlist(lapply(c3_chr_xo, function(x) return(x[["3H"]])))),
    MeanXO_4H=c(
        unlist(lapply(c1_chr_xo, function(x) return(x[["4H"]]))),
        unlist(lapply(c2_chr_xo, function(x) return(x[["4H"]]))),
        unlist(lapply(c3_chr_xo, function(x) return(x[["4H"]])))),
    MeanXO_5H=c(
        unlist(lapply(c1_chr_xo, function(x) return(x[["5H"]]))),
        unlist(lapply(c2_chr_xo, function(x) return(x[["5H"]]))),
        unlist(lapply(c3_chr_xo, function(x) return(x[["5H"]])))),
    MeanXO_6H=c(
        unlist(lapply(c1_chr_xo, function(x) return(x[["6H"]]))),
        unlist(lapply(c2_chr_xo, function(x) return(x[["6H"]]))),
        unlist(lapply(c3_chr_xo, function(x) return(x[["6H"]])))),
    MeanXO_7H=c(
        unlist(lapply(c1_chr_xo, function(x) return(x[["7H"]]))),
        unlist(lapply(c2_chr_xo, function(x) return(x[["7H"]]))),
        unlist(lapply(c3_chr_xo, function(x) return(x[["7H"]]))))
    )

#   And write it
write.table(
    pop_summary,
    file="Genomic_Prediction_Family_Genetic_Summary.txt",
    row.names=FALSE,
    quote=FALSE,
    sep="\t")
