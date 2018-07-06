# Make boxplots of the per-SNP additive and dominance components for the
# GVCBLUP output.

# Read the outputs for each trait and functional category
yld.nonc <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/GVCBLUP/Yld/Yld_Noncoding_Effects.txt", header=TRUE)
yld.syn <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/GVCBLUP/Yld/Yld_Synonymous_Effects.txt", header=TRUE)
yld.nons <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/GVCBLUP/Yld/Yld_Nonsynonymous_Effects.txt", header=TRUE)
yld.del <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/GVCBLUP/Yld/Yld_Deleterious_Effects.txt", header=TRUE)

don.nonc <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/GVCBLUP/DON/DON_Noncoding_Effects.txt", header=TRUE)
don.syn <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/GVCBLUP/DON/DON_Synonymous_Effects.txt", header=TRUE)
don.nons <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/GVCBLUP/DON/DON_Nonsynonymous_Effects.txt", header=TRUE)
don.del <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/GVCBLUP/DON/DON_Deleterious_Effects.txt", header=TRUE)

height.nonc <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/GVCBLUP/Height/Height_Noncoding_Effects.txt", header=TRUE)
height.syn <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/GVCBLUP/Height/Height_Synonymous_Effects.txt", header=TRUE)
height.nons <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/GVCBLUP/Height/Height_Nonsynonymous_Effects.txt", header=TRUE)
height.del <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/GVCBLUP/Height/Height_Deleterious_Effects.txt", header=TRUE)

# boxplot colors
box_cols <- c("black", "black", "blue", "blue", "green", "green", "red", "red")

# Put together data frames to plot the values
yld.toplot <- data.frame(
    Label=c(
        rep("NCA", nrow(yld.nonc)),
        rep("NCD", nrow(yld.nonc)),
        rep("SYA", nrow(yld.syn)),
        rep("SYD", nrow(yld.syn)),
        rep("NSA", nrow(yld.nons)),
        rep("NSD", nrow(yld.nons)),
        rep("DEA", nrow(yld.del)),
        rep("DED", nrow(yld.del))
        ),
    Value=c(
        yld.nonc$Effect_A2,
        yld.nonc$Effect_D2,
        yld.syn$Effect_A2,
        yld.syn$Effect_D2,
        yld.nons$Effect_A2,
        yld.nons$Effect_D2,
        yld.del$Effect_A2,
        yld.del$Effect_D2
        )
    )
yld.toplot$Label <- factor(yld.toplot$Label, levels=c("NCA", "NCD", "SYA",
                                                      "SYD", "NSA", "NSD",
                                                      "DEA", "DED"))

don.toplot <- data.frame(
    Label=c(
        rep("NCA", nrow(don.nonc)),
        rep("NCD", nrow(don.nonc)),
        rep("SYA", nrow(don.syn)),
        rep("SYD", nrow(don.syn)),
        rep("NSA", nrow(don.nons)),
        rep("NSD", nrow(don.nons)),
        rep("DEA", nrow(don.del)),
        rep("DED", nrow(don.del))
        ),
    Value=c(
        don.nonc$Effect_A2,
        don.nonc$Effect_D2,
        don.syn$Effect_A2,
        don.syn$Effect_D2,
        don.nons$Effect_A2,
        don.nons$Effect_D2,
        don.del$Effect_A2,
        don.del$Effect_D2
        )
    )
don.toplot$Label <- factor(don.toplot$Label, levels=c("NCA", "NCD", "SYA",
                                                      "SYD", "NSA", "NSD",
                                                      "DEA", "DED"))

height.toplot <- data.frame(
    Label=c(
        rep("NCA", nrow(height.nonc)),
        rep("NCD", nrow(height.nonc)),
        rep("SYA", nrow(height.syn)),
        rep("SYD", nrow(height.syn)),
        rep("NSA", nrow(height.nons)),
        rep("NSD", nrow(height.nons)),
        rep("DEA", nrow(height.del)),
        rep("DED", nrow(height.del))
        ),
    Value=c(
        height.nonc$Effect_A2,
        height.nonc$Effect_D2,
        height.syn$Effect_A2,
        height.syn$Effect_D2,
        height.nons$Effect_A2,
        height.nons$Effect_D2,
        height.del$Effect_A2,
        height.del$Effect_D2
        )
    )
height.toplot$Label <- factor(height.toplot$Label, levels=c("NCA", "NCD", "SYA",
                                                            "SYD", "NSA", "NSD",
                                                            "DEA", "DED"))


png(file="GVCBLUP_AD_Hists.png", height=3600, width=3600, res=300)
par(mfrow=c(3, 2), mar=c(4, 4, 2, 0), mgp=c(2, 1, 0))
h.nc <- hist(
    log(1+yld.nonc$Effect_A2),
    breaks=100,
    plot=FALSE
    )
h.nc$density <- h.nc$counts/sum(h.nc$counts)
h.s <- hist(
    log(1+yld.syn$Effect_A2),
    breaks=100,
    plot=FALSE
    )
h.s$density <- h.s$counts/sum(h.s$counts)
h.ns <- hist(
    log(1+yld.nons$Effect_A2),
    breaks=100,
    plot=FALSE
    )
h.ns$density <- h.ns$counts/sum(h.ns$counts)
h.del <- hist(
    log(1+yld.del$Effect_A2),
    breaks=100,
    plot=FALSE
    )
h.del$density <- h.del$counts/sum(h.del$counts)
plot(
    h.nc,
    freq=FALSE,
    border=NA,
    xlab="log(1+Per-SNP Additive Effect)",
    ylab="Proportion",
    ylim=c(0, 0.6),
    main="Yield",
    col=rgb(0, 0, 0, 0.25))
plot(
    h.s,
    freq=FALSE,
    border=NA,
    col=rgb(0, 0, 1, 0.25),
    add=TRUE)
plot(
    h.ns,
    freq=FALSE,
    border=NA,
    col=rgb(0, 1, 0, 0.25),
    add=TRUE)
plot(
    h.del,
    freq=FALSE,
    border=NA,
    col=rgb(1, 0, 0, 0.25),
    add=TRUE)
abline(v=0, lwd=2, col="grey", lty=3)
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"),
    fill=c(rgb(0, 0, 0, 0.25), rgb(0, 0, 1, 0.25), rgb(0, 1, 0, 0.25), rgb(1, 0, 0, 0.25)))

h.nc <- hist(
    log(1+yld.nonc$Effect_D2),
    breaks=100,
    plot=FALSE
    )
h.nc$density <- h.nc$counts/sum(h.nc$counts)
h.s <- hist(
    log(1+yld.syn$Effect_D2),
    breaks=100,
    plot=FALSE
    )
h.s$density <- h.s$counts/sum(h.s$counts)
h.ns <- hist(
    log(1+yld.nons$Effect_D2),
    breaks=100,
    plot=FALSE
    )
h.ns$density <- h.ns$counts/sum(h.ns$counts)
h.del <- hist(
    log(1+yld.del$Effect_D2),
    breaks=100,
    plot=FALSE
    )
h.del$density <- h.del$counts/sum(h.del$counts)
plot(
    h.nc,
    freq=FALSE,
    border=NA,
    xlab="log(1+Per-SNP Dominance Effect)",
    ylab="Proportion",
    ylim=c(0, 0.6),
    main="Yield",
    col=rgb(0, 0, 0, 0.25))
plot(
    h.s,
    freq=FALSE,
    border=NA,
    col=rgb(0, 0, 1, 0.25),
    add=TRUE)
plot(
    h.ns,
    freq=FALSE,
    border=NA,
    col=rgb(0, 1, 0, 0.25),
    add=TRUE)
plot(
    h.del,
    freq=FALSE,
    border=NA,
    col=rgb(1, 0, 0, 0.25),
    add=TRUE)
abline(v=0, lwd=2, col="grey", lty=3)
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"),
    fill=c(rgb(0, 0, 0, 0.25), rgb(0, 0, 1, 0.25), rgb(0, 1, 0, 0.25), rgb(1, 0, 0, 0.25)))

h.nc <- hist(
    log(1+don.nonc$Effect_A2),
    breaks=100,
    plot=FALSE
    )
h.nc$density <- h.nc$counts/sum(h.nc$counts)
h.s <- hist(
    log(1+don.syn$Effect_A2),
    breaks=100,
    plot=FALSE
    )
h.s$density <- h.s$counts/sum(h.s$counts)
h.ns <- hist(
    log(1+don.nons$Effect_A2),
    breaks=100,
    plot=FALSE
    )
h.ns$density <- h.ns$counts/sum(h.ns$counts)
h.del <- hist(
    log(1+don.del$Effect_A2),
    breaks=100,
    plot=FALSE
    )
h.del$density <- h.del$counts/sum(h.del$counts)
plot(
    h.nc,
    freq=FALSE,
    border=NA,
    xlab="log(1+Per-SNP Additive Effect)",
    ylab="Proportion",
    ylim=c(0, 0.6),
    main="DON Concentration",
    col=rgb(0, 0, 0, 0.25))
plot(
    h.s,
    freq=FALSE,
    border=NA,
    col=rgb(0, 0, 1, 0.25),
    add=TRUE)
plot(
    h.ns,
    freq=FALSE,
    border=NA,
    col=rgb(0, 1, 0, 0.25),
    add=TRUE)
plot(
    h.del,
    freq=FALSE,
    border=NA,
    col=rgb(1, 0, 0, 0.25),
    add=TRUE)
abline(v=0, lwd=2, col="grey", lty=3)
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"),
    fill=c(rgb(0, 0, 0, 0.25), rgb(0, 0, 1, 0.25), rgb(0, 1, 0, 0.25), rgb(1, 0, 0, 0.25)))

h.nc <- hist(
    log(1+don.nonc$Effect_D2),
    breaks=100,
    plot=FALSE
    )
h.nc$density <- h.nc$counts/sum(h.nc$counts)
h.s <- hist(
    log(1+don.syn$Effect_D2),
    breaks=100,
    plot=FALSE
    )
h.s$density <- h.s$counts/sum(h.s$counts)
h.ns <- hist(
    log(1+don.nons$Effect_D2),
    breaks=100,
    plot=FALSE
    )
h.ns$density <- h.ns$counts/sum(h.ns$counts)
h.del <- hist(
    log(1+don.del$Effect_D2),
    breaks=100,
    plot=FALSE
    )
h.del$density <- h.del$counts/sum(h.del$counts)
plot(
    h.nc,
    freq=FALSE,
    border=NA,
    xlab="log(1+Per-SNP Dominance Effect)",
    ylab="Proportion",
    ylim=c(0, 0.6),
    main="DON Concentration",
    col=rgb(0, 0, 0, 0.25))
plot(
    h.s,
    freq=FALSE,
    border=NA,
    col=rgb(0, 0, 1, 0.25),
    add=TRUE)
plot(
    h.ns,
    freq=FALSE,
    border=NA,
    col=rgb(0, 1, 0, 0.25),
    add=TRUE)
plot(
    h.del,
    freq=FALSE,
    border=NA,
    col=rgb(1, 0, 0, 0.25),
    add=TRUE)
abline(v=0, lwd=2, col="grey", lty=3)
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"),
    fill=c(rgb(0, 0, 0, 0.25), rgb(0, 0, 1, 0.25), rgb(0, 1, 0, 0.25), rgb(1, 0, 0, 0.25)))

h.nc <- hist(
    log(1+height.nonc$Effect_A2),
    breaks=100,
    plot=FALSE
    )
h.nc$density <- h.nc$counts/sum(h.nc$counts)
h.s <- hist(
    log(1+height.syn$Effect_A2),
    breaks=100,
    plot=FALSE
    )
h.s$density <- h.s$counts/sum(h.s$counts)
h.ns <- hist(
    log(1+height.nons$Effect_A2),
    breaks=100,
    plot=FALSE
    )
h.ns$density <- h.ns$counts/sum(h.ns$counts)
h.del <- hist(
    log(1+height.del$Effect_A2),
    breaks=50,
    plot=FALSE
    )
h.del$density <- h.del$counts/sum(h.del$counts)
plot(
    h.nc,
    freq=FALSE,
    border=NA,
    xlab="log(1+Per-SNP Additive Effect)",
    ylab="Proportion",
    ylim=c(0, 0.6),
    main="Height",
    col=rgb(0, 0, 0, 0.25))
plot(
    h.s,
    freq=FALSE,
    border=NA,
    col=rgb(0, 0, 1, 0.25),
    add=TRUE)
plot(
    h.ns,
    freq=FALSE,
    border=NA,
    col=rgb(0, 1, 0, 0.25),
    add=TRUE)
plot(
    h.del,
    freq=FALSE,
    border=NA,
    col=rgb(1, 0, 0, 0.25),
    add=TRUE)
abline(v=0, lwd=2, col="grey", lty=3)
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"),
    fill=c(rgb(0, 0, 0, 0.25), rgb(0, 0, 1, 0.25), rgb(0, 1, 0, 0.25), rgb(1, 0, 0, 0.25)))

h.nc <- hist(
    log(1+height.nonc$Effect_D2),
    breaks=100,
    plot=FALSE
    )
h.nc$density <- h.nc$counts/sum(h.nc$counts)
h.s <- hist(
    log(1+height.syn$Effect_D2),
    breaks=100,
    plot=FALSE
    )
h.s$density <- h.s$counts/sum(h.s$counts)
h.ns <- hist(
    log(1+height.nons$Effect_D2),
    breaks=100,
    plot=FALSE
    )
h.ns$density <- h.ns$counts/sum(h.ns$counts)
h.del <- hist(
    log(1+height.del$Effect_D2),
    breaks=100,
    plot=FALSE
    )
h.del$density <- h.del$counts/sum(h.del$counts)
plot(
    h.nc,
    freq=FALSE,
    border=NA,
    xlab="log(1+Per-SNP Dominance Effect)",
    ylab="Proportion",
    ylim=c(0, 0.6),
    main="Height",
    col=rgb(0, 0, 0, 0.25))
plot(
    h.s,
    freq=FALSE,
    border=NA,
    col=rgb(0, 0, 1, 0.25),
    add=TRUE)
plot(
    h.ns,
    freq=FALSE,
    border=NA,
    col=rgb(0, 1, 0, 0.25),
    add=TRUE)
plot(
    h.del,
    freq=FALSE,
    border=NA,
    col=rgb(1, 0, 0, 0.25),
    add=TRUE)
abline(v=0, lwd=2, col="grey", lty=3)
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"),
    fill=c(rgb(0, 0, 0, 0.25), rgb(0, 0, 1, 0.25), rgb(0, 1, 0, 0.25), rgb(1, 0, 0, 0.25)))
dev.off()

png(file="GVCBLUP_AD_Plots.png", height=3600, width=3600, res=300)
par(mfrow=c(2, 2), mar=c(4, 4, 2, 0), mgp=c(2, 1, 0))
boxplot(
    yld.toplot$Value ~ yld.toplot$Label,
    beside=TRUE,
    at=c(1, 2, 4, 5, 7, 8, 10, 11),
    border=box_cols,
    lwd=2,
    lty=c(1, 3),
    axes=FALSE,
    xlab="",
    ylab="Squared Effect Size",
    main="Yield")
axis(side=2)
axis(
    side=1,
    at=c(1.5, 4.5, 7.5, 10.5),
    labels=c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"))
legend("topright", c("Additive", "Dominance"), col="black", lwd=2, lty=c(1, 3))

boxplot(
    don.toplot$Value ~ don.toplot$Label,
    beside=TRUE,
    at=c(1, 2, 4, 5, 7, 8, 10, 11),
    border=box_cols,
    lwd=2,
    lty=c(1, 3),
    axes=FALSE,
    xlab="",
    ylab="Squared Effect Size",
    main="DON Concentration")
axis(side=2)
axis(
    side=1,
    at=c(1.5, 4.5, 7.5, 10.5),
    labels=c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"))
legend("topright", c("Additive", "Dominance"), col="black", lwd=2, lty=c(1, 3))

boxplot(
    height.toplot$Value ~ height.toplot$Label,
    beside=TRUE,
    at=c(1, 2, 4, 5, 7, 8, 10, 11),
    border=box_cols,
    lwd=2,
    lty=c(1, 3),
    axes=FALSE,
    xlab="",
    ylab="Squared Effect Size",
    main="Height")
axis(side=2)
axis(
    side=1,
    at=c(1.5, 4.5, 7.5, 10.5),
    labels=c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"))
legend("topright", c("Additive", "Dominance"), col="black", lwd=2, lty=c(1, 3))
dev.off()
