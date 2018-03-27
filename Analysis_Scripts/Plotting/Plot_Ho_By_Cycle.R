# Script to plot distributions of observed heterozygosity across cycles and
# across functional classes.

dat <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Summaries/Ho_ByCycle_ByFunc.txt", header=TRUE)
# Reorder the levels of the functional class
dat$FC <- factor(dat$FC, levels=c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"))

# Separate them by cycle
c1 <- dat[dat$Cycle=="C1",]
c2 <- dat[dat$Cycle=="C2",]
c3 <- dat[dat$Cycle=="C3",]

# THen make a plot by functional class
pdf(file="Ho_ByCycle_ByFunc.pdf", height=10.5, width=8)
par(mfrow=c(3, 1), mar=c(4, 4, 1, 2))
boxplot(c1$PHet ~ c1$FC, border=c("black", "blue", "green", "red"), xlab="Functional Class", ylab="Proportion Heterozygous", main="Cycle 1")
boxplot(c2$PHet ~ c2$FC, border=c("black", "blue", "green", "red"), xlab="Functional Class", ylab="Proportion Heterozygous", main="Cycle 2")
boxplot(c3$PHet ~ c3$FC, border=c("black", "blue", "green", "red"), xlab="Functional Class", ylab="Proportion Heterozygous", main="Cycle 3")
dev.off()
