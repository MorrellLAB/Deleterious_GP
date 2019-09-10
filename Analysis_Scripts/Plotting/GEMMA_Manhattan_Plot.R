# Make a Manhattan plot of -log(P_lrt) from the GEMMA results. These will be
# a little weird because the population is highly structured and the founders
# are pretty closely related.

# Read the GEMMA association results
yld <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/GEMMA/Yield_LMM.assoc.txt.gz", header=TRUE)
don <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/GEMMA/DON_LMM.assoc.txt.gz", header=TRUE)
hgt <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/GEMMA/Height_LMM.assoc.txt.gz", header=TRUE)

# Digest up the plotting a little bit. We want to plot the seven chromosomes
# side-by-side and give them distinct colors
# This is "Dark2" with seven levels from RColorBrewer
color_vector <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d")

# We will re-define the data to plot so that the chromosomes will be placed
# end-to-end. Divide by 1M to put it in megabase scale
chr1H_end <- as.numeric(tail(yld$ps[yld$chr == "chr1H"], n=1))
chr2H_end <- as.numeric(tail(yld$ps[yld$chr == "chr2H"], n=1)) + chr1H_end
chr3H_end <- as.numeric(tail(yld$ps[yld$chr == "chr3H"], n=1)) + chr2H_end
chr4H_end <- as.numeric(tail(yld$ps[yld$chr == "chr4H"], n=1)) + chr3H_end
chr5H_end <- as.numeric(tail(yld$ps[yld$chr == "chr5H"], n=1)) + chr4H_end
chr6H_end <- as.numeric(tail(yld$ps[yld$chr == "chr6H"], n=1)) + chr5H_end

yld_plt <- data.frame(
    X=c(
        yld$ps[yld$chr == "chr1H"],
        yld$ps[yld$chr == "chr2H"] + chr1H_end,
        yld$ps[yld$chr == "chr3H"] + chr2H_end,
        yld$ps[yld$chr == "chr4H"] + chr3H_end,
        yld$ps[yld$chr == "chr5H"] + chr4H_end,
        yld$ps[yld$chr == "chr6H"] + chr5H_end,
        yld$ps[yld$chr == "chr7H"] + chr6H_end) / 1000000,
    Y=-log10(p.adjust(yld$p_lrt, method="BH")),
    Chrom=c(
        rep("Chr1H", sum(yld$chr == "chr1H")),
        rep("Chr2H", sum(yld$chr == "chr2H")),
        rep("Chr3H", sum(yld$chr == "chr3H")),
        rep("Chr4H", sum(yld$chr == "chr4H")),
        rep("Chr5H", sum(yld$chr == "chr5H")),
        rep("Chr6H", sum(yld$chr == "chr6H")),
        rep("Chr7H", sum(yld$chr == "chr7H"))
        )
    )

don_plt <- data.frame(
    X=c(
        don$ps[don$chr == "chr1H"],
        don$ps[don$chr == "chr2H"] + chr1H_end,
        don$ps[don$chr == "chr3H"] + chr2H_end,
        don$ps[don$chr == "chr4H"] + chr3H_end,
        don$ps[don$chr == "chr5H"] + chr4H_end,
        don$ps[don$chr == "chr6H"] + chr5H_end,
        don$ps[don$chr == "chr7H"] + chr6H_end) / 1000000,
    Y=-log10(p.adjust(don$p_lrt, method="BH")),
    Chrom=c(
        rep("Chr1H", sum(don$chr == "chr1H")),
        rep("Chr2H", sum(don$chr == "chr2H")),
        rep("Chr3H", sum(don$chr == "chr3H")),
        rep("Chr4H", sum(don$chr == "chr4H")),
        rep("Chr5H", sum(don$chr == "chr5H")),
        rep("Chr6H", sum(don$chr == "chr6H")),
        rep("Chr7H", sum(don$chr == "chr7H"))
        )
    )

hgt_plt <- data.frame(
    X=c(
        hgt$ps[hgt$chr == "chr1H"],
        hgt$ps[hgt$chr == "chr2H"] + chr1H_end,
        hgt$ps[hgt$chr == "chr3H"] + chr2H_end,
        hgt$ps[hgt$chr == "chr4H"] + chr3H_end,
        hgt$ps[hgt$chr == "chr5H"] + chr4H_end,
        hgt$ps[hgt$chr == "chr6H"] + chr5H_end,
        hgt$ps[hgt$chr == "chr7H"] + chr6H_end) / 1000000,
    Y=-log10(p.adjust(hgt$p_lrt, method="BH")),
    Chrom=c(
        rep("Chr1H", sum(hgt$chr == "chr1H")),
        rep("Chr2H", sum(hgt$chr == "chr2H")),
        rep("Chr3H", sum(hgt$chr == "chr3H")),
        rep("Chr4H", sum(hgt$chr == "chr4H")),
        rep("Chr5H", sum(hgt$chr == "chr5H")),
        rep("Chr6H", sum(hgt$chr == "chr6H")),
        rep("Chr7H", sum(hgt$chr == "chr7H"))
        )
    )

# pdf(file="GEMMA_Manhattan.pdf", height=4, width=8)
png(file="GEMMA_Manhattan.png", res=150, height=600, width=1200)
par(mfrow=c(3, 1), mar=c(1, 3, 0.75, 0.1), mgp=c(2, 1, 0))
chroms <- c("Chr1H", "Chr2H", "Chr3H", "Chr4H", "Chr5H", "Chr6H", "Chr7H")
# Make a plot for Yield
plot(0, type="n", axes=FALSE, ylim=c(0, 4), xlim=c(0, 4600), xlab="", ylab="-log10(P)", main="Yield (kg/ha)")
sapply(
    seq_along(chroms),
    function(x) {
        cname <- chroms[x]
        color <- color_vector[x]
        points(
            yld_plt$Y[yld_plt$Chrom == cname] ~ yld_plt$X[yld_plt$Chrom == cname],
            pch=19,
            cex=0.25,
            col=color)
    })
abline(h=-log10(0.01), col="red", lwd=0.75, lty=3)
abline(h=-log10(0.1), col="blue", lwd=0.75, lty=3)
abline(
    v=c(chr1H_end, chr2H_end, chr3H_end, chr4H_end, chr5H_end, chr6H_end)/1000000,
    lwd=0.25,
    col="grey",
    lty=2)
axis(side=2)

# Make a plot for DON
plot(0, type="n", axes=FALSE, ylim=c(0, 4), xlim=c(0, 4600), xlab="", ylab="-log10(P)", main="DON (ppm)")
sapply(
    seq_along(chroms),
    function(x) {
        cname <- chroms[x]
        color <- color_vector[x]
        points(
            don_plt$Y[don_plt$Chrom == cname] ~ don_plt$X[don_plt$Chrom == cname],
            pch=19,
            cex=0.25,
            col=color)
    })
abline(h=-log10(0.01), col="red", lwd=0.75, lty=3)
abline(h=-log10(0.1), col="blue", lwd=0.75, lty=3)
abline(
    v=c(chr1H_end, chr2H_end, chr3H_end, chr4H_end, chr5H_end, chr6H_end)/1000000,
    lwd=0.25,
    col="grey",
    lty=2)
axis(side=2)

# Make a plot for Height
plot(0, type="n", axes=FALSE, ylim=c(0, 4), xlim=c(0, 4600), xlab="", ylab="-log10(P)", main="Height (cm)")
sapply(
    seq_along(chroms),
    function(x) {
        cname <- chroms[x]
        color <- color_vector[x]
        points(
            hgt_plt$Y[hgt_plt$Chrom == cname] ~ hgt_plt$X[hgt_plt$Chrom == cname],
            pch=19,
            cex=0.25,
            col=color)
    })
abline(h=-log10(0.01), col="red", lwd=0.75, lty=3)
abline(h=-log10(0.1), col="blue", lwd=0.75, lty=3)
abline(
    v=c(chr1H_end, chr2H_end, chr3H_end, chr4H_end, chr5H_end, chr6H_end)/1000000,
    lwd=0.25,
    col="grey",
    lty=2)
axis(side=2)

mtext(
    sub("Chr", "", chroms),
    side=1,
    at=c(
        mean(c(0, chr1H_end)),
        mean(c(chr1H_end, chr2H_end)),
        mean(c(chr2H_end, chr3H_end)),
        mean(c(chr3H_end, chr4H_end)),
        mean(c(chr4H_end, chr5H_end)),
        mean(c(chr5H_end, chr6H_end)),
        mean(c(chr6H_end, 4600000000))
        )/1000000,
    cex=0.75
    )

dev.off()
