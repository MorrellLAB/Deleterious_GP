# Make a plot of the proportion of fixed sites by class and by C1 frequency.

nc_col <- '#2c7bb6'
syn_col <- '#abd9e9'
ns_col <- '#fdae61'
del_col <- '#d7191c'

# Read in the genotype frequency data
par_frq <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/Ancestral_State/GP_Ancestral.txt.gz", header=TRUE)
par_frq <- par_frq[,c("SNPID", "DAF")]
names(par_frq) <- c("SNP_ID", "C0_DAF")
freq <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Genotype_Freqs/DAF_By_Cycle.txt.gz", header=TRUE)
freq <- merge(freq, par_frq, by="SNP_ID")
# Remove rows with NA
keep <- apply(freq, 1, function(x) {
    if(any(is.na(x))) {
        return(FALSE)
    } else if(all(as.numeric(x[c("C0_DAF", "C1_DAF", "C2_DAF", "C3_DAF")]) == 0)) {
        return(FALSE)
    } else if(all(as.numeric(x[c("C0_DAF", "C1_DAF", "C2_DAF", "C3_DAF")]) == 1)) {
        return(FALSE)
    } else {
        return(TRUE)
    }
})
dim(freq)
freq <- freq[keep,]
dim(freq)
# Read in the functional class lists
nonc <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Noncoding.names.gz", header=FALSE)$V1)
syn <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Synonymous.names.gz", header=FALSE)$V1)
nonsyn <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Nonsynonymous.names.gz", header=FALSE)$V1)
del <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Deleterious.names.gz", header=FALSE)$V1)

# Slice down the matrix to the functional classes
nonc.dafs <- freq[freq$SNP_ID %in% nonc,]
syn.dafs <- freq[freq$SNP_ID %in% syn,]
nonsyn.dafs <- freq[freq$SNP_ID %in% nonsyn,]
del.dafs <- freq[freq$SNP_ID %in% del,]

# Make breakpoints for the frequency bins.
breaks <- seq(0, 1, by=0.1)

# Define a function to return 0/1 for the loss
loss <- function(x) {
    if(is.na(x)) {
        return(NA)
    } else if(x == 0.0) {
        return(1)
    } else {
        return(0)
    }
}

# Then cut the C1 frequencies by the bins
nonc.poly <- data.frame(
    Bin=cut(nonc.dafs$C0_DAF, breaks=breaks, include.lowest=TRUE),
    Lost=sapply(nonc.dafs$C3_DAF, loss)
    )

syn.poly <- data.frame(
    Bin=cut(syn.dafs$C0_DAF, breaks=breaks, include.lowest=TRUE),
    Lost=sapply(syn.dafs$C3_DAF, loss)
    )

nonsyn.poly <- data.frame(
    Bin=cut(nonsyn.dafs$C0_DAF, breaks=breaks, include.lowest=TRUE),
    Lost=sapply(nonsyn.dafs$C3_DAF, loss)
    )

del.poly <- data.frame(
    Bin=cut(del.dafs$C0_DAF, breaks=breaks, include.lowest=TRUE),
    Lost=sapply(del.dafs$C3_DAF, loss)
    )

# Count up the lost/not by frequency bin
nonc.tab <- table(nonc.poly)
syn.tab <- table(syn.poly)
nonsyn.tab <- table(nonsyn.poly)
del.tab <- table(del.poly)

# Convert to proportion
nonc.lprop <- apply(nonc.tab, 1, function(x) {
    return(x["1"]/(x["0"] + x["1"]))
    })
syn.lprop <- apply(syn.tab, 1, function(x) {
    return(x["1"]/(x["0"] + x["1"]))
    })
nonsyn.lprop <- apply(nonsyn.tab, 1, function(x) {
    return(x["1"]/(x["0"] + x["1"]))
    })
del.lprop <- apply(del.tab, 1, function(x) {
    return(x["1"]/(x["0"] + x["1"]))
    })

# Make a plot!
pdf(file="Prop_Loss_by_DAF.pdf", height=6, width=6)
par(mar=c(4, 4, 0, 0), mgp=c(2, 1, 0))
plot(
    nonc.lprop,
    type="b",
    col=nc_col,
    lwd=3,
    pch=19,
    xlab="DAF in Parents",
    ylab="Proportion of Variants Fixed for Ancestral Allele",
    axes=F,
    ylim=c(0, 0.5))
lines(syn.lprop, col=syn_col, lwd=3, pch=19, type="b")
lines(nonsyn.lprop, col=ns_col, lwd=3, pch=19, type="b")
lines(del.lprop, col=del_col, lwd=3, pch=19, type="b")
axis(side=2)
axis(side=1, at=1:(length(breaks)-1), labels=levels(nonc.poly$Bin))
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"), col=c(nc_col, syn_col, ns_col, del_col), pch=19)
dev.off()

# Let's plot the mean change in frequency by DAF in C1
delta_daf <- function(x) {
    return(as.numeric(x["C3_DAF"]) - as.numeric(x["C0_DAF"]))
}

nonc.ddaf <- data.frame(
    Bin=cut(nonc.dafs$C0_DAF, breaks=breaks, include.lowest=TRUE),
    DDAF=apply(nonc.dafs, 1, delta_daf)
    )

syn.ddaf <- data.frame(
    Bin=cut(syn.dafs$C0_DAF, breaks=breaks, include.lowest=TRUE),
    DDAF=apply(syn.dafs, 1, delta_daf)
    )

nonsyn.ddaf <- data.frame(
    Bin=cut(nonsyn.dafs$C0_DAF, breaks=breaks, include.lowest=TRUE),
    DDAF=apply(nonsyn.dafs, 1, delta_daf)
    )

del.ddaf <- data.frame(
    Bin=cut(del.dafs$C0_DAF, breaks=breaks, include.lowest=TRUE),
    DDAF=apply(del.dafs, 1, delta_daf)
    )

nonc.mddaf <- sapply(levels(nonc.ddaf$Bin), function(x) {
    k <- nonc.ddaf[nonc.ddaf$Bin == x, ]
    return(mean(k$DDAF))})
syn.mddaf <- sapply(levels(syn.ddaf$Bin), function(x) {
    k <- syn.ddaf[syn.ddaf$Bin == x, ]
    return(mean(k$DDAF))})
nonsyn.mddaf <- sapply(levels(nonsyn.ddaf$Bin), function(x) {
    k <- nonsyn.ddaf[nonsyn.ddaf$Bin == x, ]
    return(mean(k$DDAF))})
del.mddaf <- sapply(levels(del.ddaf$Bin), function(x) {
    k <- del.ddaf[del.ddaf$Bin == x, ]
    return(mean(k$DDAF))})

pdf(file="DeltaDAF_by_DAF.pdf", height=6, width=6)
par(mar=c(4, 4, 0, 0), mgp=c(2, 1, 0))
plot(
    nonc.mddaf,
    type="b",
    col=nc_col,
    lwd=3,
    pch=19,
    xlab="DAF in Parents",
    ylab="Median Change in DAF From C1 to C3",
    axes=F,
    ylim=c(-0.25, 0.25))
lines(syn.mddaf, col=syn_col, lwd=3, pch=19, type="b")
lines(nonsyn.mddaf, col=ns_col, lwd=3, pch=19, type="b")
lines(del.mddaf, col=del_col, lwd=3, pch=19, type="b")
abline(h=0, lwd=1, lty=2)
axis(side=2)
axis(side=1, at=1:(length(breaks)-1), labels=levels(nonc.poly$Bin))
legend("top", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"), col=c(nc_col, syn_col, ns_col, del_col), pch=19)
dev.off()


# How about the proportion of variants decreasing DAF?
decreasing <- function(x) {
    if(as.numeric(x["C3_DAF"]) < as.numeric(x["C0_DAF"])) {
        return(1)
    } else {
        return(0)
    }
}

nonc.dec <- data.frame(
    Bin=cut(nonc.dafs$C0_DAF, breaks=breaks, include.lowest=TRUE),
    DDAF=apply(nonc.dafs, 1, decreasing)
    )

syn.dec <- data.frame(
    Bin=cut(syn.dafs$C0_DAF, breaks=breaks, include.lowest=TRUE),
    DDAF=apply(syn.dafs, 1, decreasing)
    )

nonsyn.dec <- data.frame(
    Bin=cut(nonsyn.dafs$C0_DAF, breaks=breaks, include.lowest=TRUE),
    DDAF=apply(nonsyn.dafs, 1, decreasing)
    )

del.dec <- data.frame(
    Bin=cut(del.dafs$C0_DAF, breaks=breaks, include.lowest=TRUE),
    DDAF=apply(del.dafs, 1, decreasing)
    )

nonc.dectab <- table(nonc.dec)
syn.dectab <- table(syn.dec)
nonsyn.dectab <- table(nonsyn.dec)
del.dectab <- table(del.dec)

# Convert to proportion
nonc.decprop <- apply(nonc.dectab, 1, function(x) {
    return(x["1"]/(x["0"] + x["1"]))
    })
syn.decprop <- apply(syn.dectab, 1, function(x) {
    return(x["1"]/(x["0"] + x["1"]))
    })
nonsyn.decprop <- apply(nonsyn.dectab, 1, function(x) {
    return(x["1"]/(x["0"] + x["1"]))
    })
del.decprop <- apply(del.dectab, 1, function(x) {
    return(x["1"]/(x["0"] + x["1"]))
    })

# Make a plot!
pdf(file="Prop_Decreasing_by_DAF.pdf", height=6, width=6)
par(mar=c(4, 4, 0, 0), mgp=c(2, 1, 0))
plot(
    nonc.decprop,
    type="b",
    col=nc_col,
    lwd=3,
    pch=19,
    xlab="DAF in Parents",
    ylab="Proportion of Variants Decreasing in DAF",
    axes=F,
    ylim=c(0, 1))
lines(syn.decprop, col=syn_col, lwd=3, pch=19, type="b")
lines(nonsyn.decprop, col=ns_col, lwd=3, pch=19, type="b")
lines(del.decprop, col=del_col, lwd=3, pch=19, type="b")
axis(side=2)
axis(side=1, at=1:(length(breaks)-1), labels=levels(nonc.poly$Bin))
legend("topright", c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"), col=c(nc_col, syn_col, ns_col, del_col), pch=19)
dev.off()
