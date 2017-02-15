#   Script to plot the overlap of GATK and FreeBayes SNPs on chromosome 7H in
#   the genomic prediction dataset. These are calculated with VCFTools with the
#   following parameters:
#       --vcf GATK_7H.vcf
#       --out GATK_FreeBayes_Comparison
#       --diff Genomic_Prediction_chr7H_2016-12-07.vcf
#       --diff-site
#   The discordance file is generated with the following parameters:
#       --vcf GATK_7H.vcf
#       --out GATK_FreeBayes_Comparison_Diploid
#       --diff Genomic_Prediction_chr7H_2016-12-19.vcf
#       --diff-site-discordance

#   use the scales library for transparency
library(scales)

#gt_presence <- read.table("/Volumes/LaCie/Genomic_Prediction/GATK_FB_Comparison.diff.sites_in_files", header=T)
#gt_discordance <- read.table("/Volumes/LaCie/Genomic_Prediction/GATK_FB_Comparison.diff.sites", header=T)

gt_presence <- read.table("/Volumes/Data_Disk/tmp/GP_NewCalls/G_F_Capture_7H.sites_in_files", header=T)
gt_discordance <- read.table("/Volumes/Data_Disk/tmp/GP_NewCalls/G_F_Capture_7H.sites", header=T)


#   We have to merge the positions columns to plot the data together
positions <- apply( 
    gt_presence[c("POS1", "POS2")],
    1,
    function(x) {
        if(x["POS1"] == x["POS2"]) {
            return(x["POS1"])
        }
        else if(x["POS2"] == ".") {
            return(x["POS1"])
        }
        else if(x["POS1"] == ".") {
            return(x["POS2"])
        }
    })

#   And we have to convert the IN_FILE variable to a value we can plot
infile <- sapply(
    gt_presence$IN_FILE,
    function(x) {
        if(x == 1) {
            return(0)
            }
        else if(x == 2) {
            return(2)
        }
        else if(x == "B") {
            return(1)
        }
    })

#pdf(file="GATK_FreeBayes_Comparison_Diploid_Parts.pdf", width=11, height=4)
png(file="GATK_FreeBayes_Comparison_Diploid_Parts.png", width=1600, height=400)
y <- jitter(infile)
x <- as.numeric(positions)/1000000
#   argument order: bottom, left, top, right
par(oma=c(0,2,0,0))
plot(y~x, pch=19, cex=0.1, axes=F, xlab="Physical Position (Mb)", ylab="", main="GATK FreeBayes Comparison on 7H", col=alpha("black", 0.2))
abline(v=((2**29)-1)/1000000, lwd=2, col="red")
axis(side=1, at=seq(0, 700, 50))
axis(side=2, at=c(0, 1, 2), labels=c("GATK\nOnly", "Both", "FreeBayes\nOnly"), las=2)
dev.off()


#   Plot discordance between the two across physical distance
#   We just want the sites that are in common between the two methods
common_discordance <- gt_discordance[gt_discordance$FILES == "B", ]
pdf(file="GATK_FreeBayes_Discordance_Diploid_Parts.pdf", width=11, height=4)
y <- jitter(common_discordance$DISCORDANCE)
x <- as.numeric(common_discordance$POS)/1000000
plot(y~x, pch=19, cex=0.01, axes=F, xlab="Physical Position (Mb)", ylab="Proportion Discordant Genotypes", main="Discordance in Genotypes Between GATK and FreeBayes on 7H")
lines(lowess(x, y), col="blue")
abline(v=((2**29)-1)/1000000, lwd=2, col="red")
axis(side=1, at=seq(0, 700, 50))
axis(side=2, at=seq(0, 1, 0.1))
dev.off()
