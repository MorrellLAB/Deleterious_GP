#   Make a fancy plot of -log(p) from the LRT and recombination rate

library(ggplot2)

#   Read the recombination rate table and the predictions
rec_rate <- read.table("/Volumes/DataDisk/Dropbox/GitHub/Deleterious_GP/Data/SNP_Positions/BOPA_cMMb_Smoothed.txt", header=T)
effects <- read.table("/Volumes/DataDisk/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GenomicPrediction_Effcts_PROVEAN_PPH_BM.txt", header=T)
excap <- read.table("/Volumes/DataDisk/Dropbox/GitHub/Deleterious_GP/Data/Resequencing_Summaries/ExomeCaptureTargets_per_Mb.txt", header=T)
bopa <- read.table("/Volumes/DataDisk/Dropbox/GitHub/Deleterious_GP/Data/SNP_Positions/384-Capture_Positions.txt", header=TRUE)

#   Drop the unmapped chromosome
chroms <- c("chr1H", "chr2H", "chr3H", "chr4H", "chr5H", "chr6H", "chr7H")
rec_rate <- rec_rate[rec_rate$Chromosome %in% chroms,]
effects <- effects[effects$Chromosome %in% chroms,]
excap <- excap[excap$Chromosome %in% chroms,]
bopa <- bopa[bopa$Chromosome %in% chroms,]

rec_rate$Chromosome <- factor(rec_rate$Chromosome, levels=chroms)
effects$Chromosome <- factor(effects$Chromosome, levels=chroms)
excap$Chromosome <- factor(excap$Chromosome, levels=chroms)
bopa$Chromosome <- factor(bopa$Chromosome, levels=chroms)

#   Get get the deleterious SNPs
#   PROVEAN parameters determined by power analysis, LL
provean_cutoff <- -4.1528
sig <- 0.05

pph.del <- effects$PPH2 == "deleterious"
pro.del <- as.numeric(effects$PROVEAN) <= provean_cutoff
lrt.del <- effects$LogisticP_Masked <= sig
nonsense <- effects[(effects$Silent == "No" & effects$AA1 != "*" & effects$AA2 == "*"), c("Chromosome", "Position")]

deleterious <- effects[(lrt.del & pro.del & pph.del), c("Chromosome", "Position", "LogisticP_Masked")]

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
deleterious$LogisticP_Masked[deleterious$LogisticP_Masked==0] <- 1e-18
#   Set really high cM/Mb values to NA
rec_rate$Smoothed_cMMb[rec_rate$Smoothed_cMMb > 10] <- NA
#   Get rid of rows with NA for position (shouldn't happen, but ...)
effects$Position <- as.numeric(effects$Position)
effects <- effects[!is.na(effects$Position),]

pdf(file="RecRate_AllSNPs.pdf", 10, 6)
ggplot(effects) +
    geom_vline(aes(xintercept=Position/1000000), size=0.05, alpha=0.1, color="#a6cee3") +
    geom_line(aes(x=(Start+End)/2000000, y=NExCap/10), data=excap, size=0.75, color="#1f78b4", alpha=0.7) +
    geom_line(aes(x= (LeftBP+RightBP)/2000000, y=Smoothed_cMMb), data=rec_rate, color="#329f2a", size=1.1, alpha=0.7) +
    geom_point(aes(x=Position/1000000, y=0.25), size=2, color="#984ea3", shape=17, data=bopa) +
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

deleterious <- deleterious[deleterious$Chromosome %in% chroms,]
deleterious$Chromosome <- factor(deleterious$Chromosome, levels=chroms)
deleterious$LogisticP_Masked <- -log(deleterious$LogisticP_Masked, base=2)
deleterious$LogisticP_Masked[deleterious$LogisticP_Masked > 10] <- 10

pdf(file="RecRate_LRT_DelOnly2.pdf", 10, 6)
ggplot(deleterious, aes(x=Position/1000000, y=LogisticP_Masked)) +
    geom_line(aes(x=(Start+End)/2000000, y=NExCap/15), data=excap, size=1, color="#1f78b4", alpha=1) +
    geom_line(aes(x=(LeftBP+RightBP)/2000000, y=Smoothed_cMMb), data=rec_rate, color="#329f2a", size=1, alpha=1) +
    geom_point(pch=19, size=0.75, alpha=0.5, color="#ef8a62") +
    geom_point(aes(x=Position/1000000, y=10), data=nonsense, color="#cc0000", size=0.75, alpha=0.5) +
    geom_point(aes(x=Position/1000000, y=0.25), size=2, color="#984ea3", shape=17, data=bopa) +
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
# y <- -log(deleterious$LogisticP_Masked, base=100)
# x <- deleterious$Recomb
# plot(y~x, pch=19, cex=0.5, main="LRT P-value and Recombination Rate", xlab="Recombination Rate (cM/Mb)", ylab="-log100(P-value)")
# dev.off()
