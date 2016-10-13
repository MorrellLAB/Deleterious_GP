#   Script to plot the heatmaps for phenotypic data from the genomic prediction
#   experiment.

#   Define a function to build a matrix of phenotypic data
phenotype_map <- function(yl, dat, adj, variable) {
    #   Trim the data down to just the year-location we are interested in
    year_loc <- dat[dat$trial == yl, c("row", "column", "line_name", "Type")]
    #   Start an empty matrix to fill with data
    pheno_mat <- matrix(0, nrow=max(year_loc$row), ncol=max(year_loc$column))
    #   Start empty vectors for check data
    check_row <- c()
    check_col <- c()
    check_names <- c()
    #   Fill in the matrices
    for(i in 1:nrow(year_loc)) {
        mat.r <- as.numeric(year_loc[i, "row"])
        mat.c <- as.numeric(year_loc[i, "column"])
        chk <- as.character(year_loc[i, "Type"])
        name <- as.character(year_loc[i, "line_name"])
        adj_row <- adj$line_name == name
        yld <- as.numeric(adj[adj_row, variable])
        #   Change NA type to something else
        if(is.na(chk)) {
            chk <- "None"
        }
        pheno_mat[mat.r, mat.c] <- yld
        #   If the type is a check line, then put that into the checks dataframe
        if(chk == "Chk") {
            check_row <- c(check_row, mat.r)
            check_col <- c(check_col, mat.c)
            check_names <- c(check_names, name)
        }
    }
    #   Convert it to a z-score so we can easily compare them
    pheno_mat <- (pheno_mat - mean(pheno_mat, na.rm=T))/sd(pheno_mat, na.rm=T)
    #   Put the checks into a data frame
    checks <- data.frame(Row=check_row, Col=check_col, Name=check_names)
    return(list(Mat=pheno_mat, Checks=checks))
}

#   Define a function to make a plot. Simple.
plot_hmp <- function(yl, phenmaps, variable) {
    #   Pull out the proper heatmap that we want to plot
    to_plot <- phenmaps[[yl]][["Mat"]]
    checks <- phenmaps[[yl]][["Checks"]]
    #   Get the minimum and maximum defined values, as they will serve for the
    #   breakpoints later.
    def_min <- min(to_plot, na.rm=TRUE)
    def_max <- max(to_plot, na.rm=TRUE)
    #   We assign NA values to something ridiculous, so that we can fill those
    #   cells with grey.
    to_plot[is.na(to_plot)] <- -999
    #   Define the colors. We will use the blue-white-red colors from
    #   Rcolorbrewer, with 11 levels. These go from dark red to white to dark
    #   blue. Should correspond to high-avg-low
    colors <- c(
        "#67001f",
        "#b2182b",
        "#d6604d",
        "#f4a582",
        "#fddbc7",
        "#f7f7f7",
        "#d1e5f0",
        "#92c5de",
        "#4393c3",
        "#2166ac",
        "#053061",
        "#aaaaaa")
    #   Make the plot
    pdf(
        file=paste(yl, "_Adjusted", variable, ".pdf", sep=""),
        width=10,
        height=8)
    #   Plot the heatmap
    image(
        x=1:ncol(to_plot),
        y=1:nrow(to_plot),
        z=t(to_plot),
        col=rev(colors),
        xlab="Column",
        ylab="Row",
        main=paste(yl, "Adjusted", variable, "Data", sep=" "),
        breaks=c(-999, seq(def_min, def_max, length.out=12)),
        axes=FALSE)
    #   Put axes on it
    axis(side=1, at=1:ncol(to_plot), labels=1:ncol(to_plot))
    axis(side=2, at=1:nrow(to_plot), labels=1:nrow(to_plot))
    #   Draw lines for the blocks
    abline(v=seq(0, 40, by=5) + 0.5, col="black", lwd=5)
    abline(h=seq(0, 9, by=3) + 0.5, col="black", lwd=5)
    #   Draw segments for the checks
    for(check in 1:nrow(checks)) {
        x <- checks[check, "Col"]
        y <- checks[check, "Row"]
        #   Get the points of the line segments
        xs <- c(x-0.5, x-0.5, x-0.5, x+0.5)
        ys <- c(y-0.5, y-0.5, y+0.5, y+0.5)
        xe <- c(x+0.5, x-0.5, x+0.5, x+0.5)
        ye <- c(y-0.5, y+0.5, y+0.5, y-0.5)
        segments(xs, ys, xe, ye, col="#ffff00", lwd=5)
        text(x, y, checks[check, "Name"], cex=0.75, srt=90)
    }
    dev.off()
}

#   Read in the data files
raw_yield <- "/Volumes/Data_Disk/Dropbox/Projects/DM_GenomicPrediction/Data/Phenotypic_Data/Raw_Yield_Data.csv"
raw_don <- "/Volumes/Data_Disk/Dropbox/Projects/DM_GenomicPrediction/Data/Phenotypic_Data/Raw_DON_Data.csv"
adjusted <- "/Volumes/Data_Disk/Dropbox/Projects/DM_GenomicPrediction/Data/Phenotypic_Data/Adjusted_Phenotypic_Data.csv"

raw_yield <- read.csv(raw_yield, header=TRUE)
raw_don <- read.csv(raw_don, header=TRUE)
adjusted <- read.csv(adjusted, header=TRUE)

#   Get the trial names
yld_trials <- as.character(unique(raw_yield$trial))
don_trials <- as.character(unique(raw_don$trial))

#   Convert the columnar data to matrices for plotting
yield_maps <- lapply(yld_trials, phenotype_map, raw_yield, adjusted, "yld.BLUE_mv")
names(yield_maps) <- yld_trials
don_maps <- lapply(don_trials, phenotype_map, raw_don, adjusted, "DON.BLUE_mv")
names(don_maps) <- don_trials

#   Make plots of the heatmaps
lapply(yld_trials, plot_hmp, yield_maps, "Yield")
lapply(don_trials, plot_hmp, don_maps, "DON")
