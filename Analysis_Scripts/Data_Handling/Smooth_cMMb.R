#   Script to calculate Lowess-smoothed estimates of cM/Mb in barley

#   define a function to return the smoothed values
smooth <- function(chrom, mapdata, f, delta) {
    map_subset <- mapdata[mapdata$Chromosome == chrom,]
    midpoints <- (map_subset$LeftBP + map_subset$RightBP)/2
    smoothed <- lowess(map_subset$cMMb~midpoints, f=f, delta=delta)
    smoothed$y[smoothed$y < 0] <- 0
    return(smoothed)
}

setwd("/Volumes/Data_Disk/Dropbox/Projects/DM_GenomicPrediction/Data")
cmmb <- read.table("BOPA_cMMb.txt", header=T)

chromosomes <- unique(cmmb$Chromosome)
smoothed_map <- sapply(chromosomes, smooth, cmmb, f=0.02, delta=3000000)

#   Associate the smoothed values with the original map intervals
new_intervals <- sapply(
    seq_along(chromosomes),
    function(chrom) {
        map_subset <- cmmb[cmmb$Chromosome==chromosomes[chrom], ]
        map_subset$Smoothed_cMMb <- smoothed_map["y", chrom]$y
        return(map_subset)
    })

df_new_intervals <- data.frame(
    Chromosome=as.character(unlist(new_intervals["Chromosome",])),
    LeftBP=as.numeric(unlist(new_intervals["LeftBP",])),
    RightBP=as.numeric(unlist(new_intervals["RightBP",])),
    Smoothed_cMMb=as.numeric(unlist(new_intervals["Smoothed_cMMb",]))
    )

#   Save the new file!
write.table(
    df_new_intervals,
    file="BOPA_cMMb_Smoothed.txt",
    quote=F,
    row.names=F,
    sep="\t")
