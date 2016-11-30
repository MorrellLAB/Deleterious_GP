#   Make a fancy plot of -log(p) from the LRT and recombination rate

library(ggplot2)

#   Read the recombination rate table and the predictions
rec_rate <- read.table("/Volumes/Data_Disk/Dropbox/Projects/DM_GenomicPrediction/Data/BOPA_cMMb_Smoothed.txt", header=T)
effects <- read.table("/Volumes/Data_Disk/Dropbox/Projects/DM_GenomicPrediction/Results/Annotations/DMII_Effects_PPH2_PROVEAN_DAF_BM.txt", header=T)
excap <- read.table("/Volumes/Data_Disk/Dropbox/Projects/DM_GenomicPrediction/Results/Annotations/ExomeCaptureTargets_per_Mb.txt", header=T)

#   Drop the unmapped chromosome
rec_rate <- rec_rate[rec_rate$Chromosome != "chrUn",]
effects <- effects[effects$Chromosome != "chrUn",]
excap <- excap[excap$Chromosome != "chrUn",]

#   Get get the deleterious SNPs
#   PROVEAN parameters determined by power analysis, LL
provean_cutoff <- -4.1528
ncodons <- sum(effects$Silent == "No")
sig <- 0.05/ncodons
max_const <- 0.91347390
minseq <- 10

pph.del <- effects$PPH2 == "deleterious"
pro.del <- as.numeric(effects$PROVEAN) <= provean_cutoff
lrt.del <- (effects$MaskedP.value <= sig) & (effects$SeqCount >= minseq) & (effects$MaskedConstraint <= max_const)

deleterious <- effects[pph.del&pro.del&lrt.del, c("Chromosome", "Position", "MaskedP.value")]

get_recomb_rate <- function(snp) {
    chrom <- as.character(snp["Chromosome"])
    pos <- as.numeric(snp["Position"])
    sel_row <- which((rec_rate$Chromosome == chrom) & (pos > rec_rate$LeftBP) & (pos < rec_rate$RightBP))
    rrate <- rec_rate[sel_row, "Smoothed_cMMb"]
    #   Set recombination rates that are too high to NA. Same with those that
    #   do not match any given interval
    if(length(rrate) == 0) {
        return(NA)
    }
    if(is.na(rrate)) {
        return(NA)
    }
    if(rrate > 10) {
        return(NA)
    }
    return(rrate)
}

get_num_excap <- function(snp) {
    chrom <- as.character(snp["Chromosome"])
    pos <- as.numeric(snp["Position"])
    sel_row <- which((excap$Chromosome == chrom) & (pos > excap$Start) & (pos <= excap$End))
    n_excap <- excap[sel_row, "NExCap"]
    if(length(n_excap) == 0) {
        return(NA)
    }
    return(n_excap)
}

#   Set 0 values to very small
deleterious$MaskedP.value[deleterious$MaskedP.value==0] <- 1e-18
#   Set really high cM/Mb values to NA
rec_rate$Smoothed_cMMb[rec_rate$Smoothed_cMMb > 10] <- NA
#   Get rid of rows with NA for position (shouldn't happen, but)
effects$Position <- as.numeric(effects$Position)
effects <- effects[!is.na(effects$Position),]

pdf(file="RecRate_AllSNPs.pdf", 10, 6)
ggplot(effects) +
    geom_vline(aes(xintercept=Position/1000000), size=0.05, alpha=0.1, color="#a6cee3") +
    geom_line(aes(x=(Start+End)/2000000, y=NExCap/10), data=excap, size=0.75, color="#1f78b4", alpha=0.7) +
    geom_line(aes(x= (LeftBP+RightBP)/2000000, y=Smoothed_cMMb), data=rec_rate, color="#329f2a", size=1.1, alpha=0.7) +
    facet_grid(Chromosome~.) +
    scale_y_continuous(limits=c(0, 10), breaks=c(0, 5, 10)) +
    scale_x_continuous(limits=c(0, 725), breaks=seq(0, 725, by=50)) +
    theme_bw() +
    theme(
        strip.background=element_blank(),
        strip.text.y=element_text(size=10, colour="black", angle=0),
        axis.ticks.y=element_blank()) +
    labs(y="", x="Physical Position (Mb)")

dev.off()

pdf(file="RecRate_LRT_DelOnly2.pdf", 10, 6)
ggplot(deleterious, aes(x=Position/1000000, y=-log(MaskedP.value, base=100))) +
    geom_line(aes(x=(Start+End)/2000000, y=NExCap/15), data=excap, size=0.75, color="#1f78b4", alpha=0.9) +
    geom_line(aes(x= (LeftBP+RightBP)/2000000, y=Smoothed_cMMb), data=rec_rate, color="#329f2a", size=1.1, alpha=0.7) +
    geom_point(pch=19, size=1.5, alpha=0.9, color="#ef8a62") +
    facet_grid(Chromosome~.) +
    scale_y_continuous(limits=c(0, 10), breaks=c(0, 5, 10)) +
    scale_x_continuous(limits=c(0, 725), breaks=seq(0, 725, by=50)) +
    theme_bw() +
    theme(
        strip.background=element_blank(),
        strip.text.y=element_text(size=10, colour="black", angle=0),
        axis.ticks.y=element_blank()) +
    labs(y="", x="Physical Position (Mb)")
dev.off()

# pdf(file="RecRate_PValue_DelOnly.pdf", 6, 6)
# y <- -log(deleterious$MaskedP.value, base=100)
# x <- deleterious$Recomb
# plot(y~x, pch=19, cex=0.5, main="LRT P-value and Recombination Rate", xlab="Recombination Rate (cM/Mb)", ylab="-log100(P-value)")
# dev.off()
