# Make boxplots of the per-line burden of derived alleles for the four
# functional classes that we are analysing

nonc <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Burden/Full_T3_Noncoding_Burden.profile.gz", header=TRUE)
syn <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Burden/Full_T3_Synonymous_Burden.profile.gz", header=TRUE)
nonsyn <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Burden/Full_T3_Nonsynonymous_Burden.profile.gz", header=TRUE)
deleterious <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Burden/Full_T3_Deleterious_Burden.profile.gz", header=TRUE)

# Define a function to assign C0, C1, C2, C3. The value that will be passed
# to this function is the family ID
assign_cycle <- function(x) {
    if(x == 0) {
        return("C0")
    } else if(grepl("MS10S3", as.character(x))) {
        return("C1")
    } else if(grepl("MS11S2", as.character(x))) {
        return("C2")
    } else if(grepl("MS12_", as.character(x))) {
        return("C3")
    } else {
        return(NA)
    }
}

nonc$Cycle <- sapply(nonc$FID, assign_cycle)
syn$Cycle <- sapply(syn$FID, assign_cycle)
nonsyn$Cycle <- sapply(nonsyn$FID, assign_cycle)
deleterious$Cycle <- sapply(deleterious$FID, assign_cycle)

# Put the scores together into a single data frame
toplot <- data.frame(
    Value=c(nonc$SCORE, syn$SCORE, nonsyn$SCORE, deleterious$SCORE),
    Cycle=c(nonc$Cycle, syn$Cycle, nonsyn$Cycle, deleterious$Cycle),
    Func=c(
        rep("Noncoding", nrow(nonc)),
        rep("Synonymous", nrow(syn)),
        rep("Nonsynonymous", nrow(nonsyn)),
        rep("Deleterious", nrow(deleterious))
        )
    )
toplot$Cycle <- factor(toplot$Cycle, levels=c("C0", "C1", "C2", "C3"))
toplot$Func <- factor(toplot$Func, levels=c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"))

pdf(file="Burden_by_FuncClass.pdf", height=8, width=8)
boxplot(
    toplot$Value ~ toplot$Func+toplot$Cycle,
    border=c("black", "blue", "green", "red"),
    lwd=2,
    at=c(
        1, 2, 3, 4,
        6, 7, 8, 9,
        11, 12, 13, 14,
        16, 17, 18, 19
        ),
    axes=FALSE,
    ylab="Average Burden of Derived Alleles",
    ylim=c(0.1, 0.5),
    )
axis(side=2)
axis(side=1,
    at=c(2.5, 7.5, 12.5, 17.5),
    labels=c("Parents", "C1", "C2", "C3"))
legend(
    "topright",
    c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious"),
    col=c("black", "blue", "green", "red"),
    lwd=3,
    ncol=2,
    cex=0.8)
dev.off()
