# Make a boxplot of the differences in DAF in random and selected panels for
# C1, C2, C3, and all, separated by functional class

nc_col <- '#2c7bb6'
sy_col <- '#abd9e9'
ns_col <- '#fdae61'
de_col <- '#d7191c'

# Read the frequencies in the panels
ran_sel <- read.table("/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/Ancestral_State/Ran_Sel_DAFs.txt", header=TRUE)
# Read the functional annltations
nc <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Noncoding.names.gz", header=FALSE)$V1)
sy <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Synonymous.names.gz", header=FALSE)$V1)
ns <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Nonsynonymous.names.gz", header=FALSE)$V1)
de <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Deleterious.names.gz", header=FALSE)$V1)

# Slice up the frequencies
c0.nc <- ran_sel$C0_par[ran_sel$SNP_ID %in% nc]
c0.sy <- ran_sel$C0_par[ran_sel$SNP_ID %in% sy]
c0.ns <- ran_sel$C0_par[ran_sel$SNP_ID %in% ns]
c0.de <- ran_sel$C0_par[ran_sel$SNP_ID %in% de]

c1.nc <- ran_sel$C1_sel[ran_sel$SNP_ID %in% nc] - ran_sel$C1_ran[ran_sel$SNP_ID %in% nc]
c1.sy <- ran_sel$C1_sel[ran_sel$SNP_ID %in% sy] - ran_sel$C1_ran[ran_sel$SNP_ID %in% sy]
c1.ns <- ran_sel$C1_sel[ran_sel$SNP_ID %in% ns] - ran_sel$C1_ran[ran_sel$SNP_ID %in% ns]
c1.de <- ran_sel$C1_sel[ran_sel$SNP_ID %in% de] - ran_sel$C1_ran[ran_sel$SNP_ID %in% de]

c2.nc <- ran_sel$C2_sel[ran_sel$SNP_ID %in% nc] - ran_sel$C2_ran[ran_sel$SNP_ID %in% nc]
c2.sy <- ran_sel$C2_sel[ran_sel$SNP_ID %in% sy] - ran_sel$C2_ran[ran_sel$SNP_ID %in% sy]
c2.ns <- ran_sel$C2_sel[ran_sel$SNP_ID %in% ns] - ran_sel$C2_ran[ran_sel$SNP_ID %in% ns]
c2.de <- ran_sel$C2_sel[ran_sel$SNP_ID %in% de] - ran_sel$C2_ran[ran_sel$SNP_ID %in% de]

c3.nc <- ran_sel$C3_sel[ran_sel$SNP_ID %in% nc] - ran_sel$C3_ran[ran_sel$SNP_ID %in% nc]
c3.sy <- ran_sel$C3_sel[ran_sel$SNP_ID %in% sy] - ran_sel$C3_ran[ran_sel$SNP_ID %in% sy]
c3.ns <- ran_sel$C3_sel[ran_sel$SNP_ID %in% ns] - ran_sel$C3_ran[ran_sel$SNP_ID %in% ns]
c3.de <- ran_sel$C3_sel[ran_sel$SNP_ID %in% de] - ran_sel$C3_ran[ran_sel$SNP_ID %in% de]

all.nc <- ran_sel$all_sel[ran_sel$SNP_ID %in% nc] - ran_sel$all_ran[ran_sel$SNP_ID %in% nc]
all.sy <- ran_sel$all_sel[ran_sel$SNP_ID %in% sy] - ran_sel$all_ran[ran_sel$SNP_ID %in% sy]
all.ns <- ran_sel$all_sel[ran_sel$SNP_ID %in% ns] - ran_sel$all_ran[ran_sel$SNP_ID %in% ns]
all.de <- ran_sel$all_sel[ran_sel$SNP_ID %in% de] - ran_sel$all_ran[ran_sel$SNP_ID %in% de]

# Make a data frame out of it
toplot <- data.frame(
    DAF=c(
        c0.nc, c0.sy, c0.ns, c0.de,
        c1.nc, c1.sy, c1.ns, c1.de,
        c2.nc, c2.sy, c2.ns, c2.de,
        c3.nc, c3.sy, c3.ns, c3.de,
        all.nc, all.sy, all.ns, all.de
        ),
    Val=c(
        rep("c0.nc", length(c0.nc)),
        rep("c0.sy", length(c0.sy)),
        rep("c0.ns", length(c0.ns)),
        rep("c0.de", length(c0.de)),
        rep("c1.nc", length(c1.nc)),
        rep("c1.sy", length(c1.sy)),
        rep("c1.ns", length(c1.ns)),
        rep("c1.de", length(c1.de)),
        rep("c2.nc", length(c2.nc)),
        rep("c2.sy", length(c2.sy)),
        rep("c2.ns", length(c2.ns)),
        rep("c2.de", length(c2.de)),
        rep("c3.nc", length(c3.nc)),
        rep("c3.sy", length(c3.sy)),
        rep("c3.ns", length(c3.ns)),
        rep("c3.de", length(c3.de)),
        rep("all.nc", length(all.nc)),
        rep("all.sy", length(all.sy)),
        rep("all.ns", length(all.ns)),
        rep("all.de", length(all.de))
    )
)

# Set the order
toplot$Val <- factor(toplot$Val, levels=c(
        "c0.nc",
        "c0.sy",
        "c0.ns",
        "c0.de",
        "c1.nc",
        "c1.sy",
        "c1.ns",
        "c1.de",
        "c2.nc",
        "c2.sy",
        "c2.ns",
        "c2.de",
        "c3.nc",
        "c3.sy",
        "c3.ns",
        "c3.de",
        "all.nc",
        "all.sy",
        "all.ns",
        "all.de"))

# Make the plot!
pdf(file="Ran_Sel_DeltaDAF.pdf", height=4, width=8)
par(mar=c(4, 4, 2, 0.1))
boxplot(
    toplot$DAF ~ toplot$Val,
    at=c(
        1, 2, 3, 4,
        6, 7, 8, 9,
        11, 12, 13, 14,
        16, 17, 18, 19,
        21, 22, 23, 24),
    lwd=1,
    border=rep(c(nc_col, sy_col, ns_col, de_col), 5),
    pch=19,
    cex=0.6,
    axes=F,
    ylab="Sel DAF - Ran DAF")
abline(h=0, lwd=1, lty=3, col="grey")
axis(side=2)
axis(side=1,
    at=c(
        mean(c(1, 2, 3, 4)),
        mean(c(6, 7, 8, 9)),
        mean(c(11, 12, 13, 14)),
        mean(c(16, 17, 18, 19)),
        mean(c(21, 22, 23, 24))),
    labels=c("C0", "C1", "C2", "C3", "All"))
legend("topright",
    c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"),
    lwd=2,
    col=c(nc_col, sy_col, ns_col, de_col),
    ncol=2)
dev.off()
