#   Make a plot of the number of variants per Mb by various SNP calling
#   algorithms, with the exome capture target density.

excap <- read.table("/Users/tomkono/Dropbox/Projects/DM_GenomicPrediction/Results/Annotations/ExomeCaptureTargets_per_Mb.txt")
gatk <- read.table("/Users/tomkono/Dropbox/Projects/DM_GenomicPrediction/Results/Annotations/GATK_Unfiltered_Positions.txt")
freebayes <- read.table("/Users/tomkono/Dropbox/Projects/DM_GenomicPrediction/Results/Annotations/FB_Unfiltered_Positions.txt")
#   Eventually have one for SAMTools mpileup

#   
