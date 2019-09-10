# Make a plot of dSNPs/sSNPs for each window of the Munoz consensus map. We
# will plot the ratio against the smoothed recombination rate estimate.

setwd("/Users/tomkono/Dropbox/GP_Revisions")

syn_count <- read.table("BOPA_Synonymous_Count.txt", header=FALSE)
del_count <- read.table("BOPA_Deleterious_Count.txt", header=FALSE)

# Remove rows that are on the unmapped chromosome
chr_un <- grepl("chrUn", as.character(syn_count$V1))
syn_count <- syn_count[!chr_un,]
del_count <- del_count[!chr_un,]

# Extract the counts of sSNPs and dSNPs for each window
sSNPs <- as.numeric(
    unlist(
        lapply(
            strsplit(
                as.character(syn_count$V1),
                "[|]"),
            function(x) {
                as.numeric(x[1])
            }
            )
        )
    )

dSNPs <- as.numeric(
    unlist(
        lapply(
            strsplit(
                as.character(del_count$V1),
                "[|]"),
            function(x) {
                as.numeric(x[1])
            }
            )
        )
    )

# Divide the two counts. This will make a lot of NaN and some Inf values.
d_s <- dSNPs/sSNPs
keep <- is.finite(d_s)

d_s <- d_s[keep]
rec_rate <- syn_count$V4[keep]
# Cap the recombination rate at 20 cM/Mb
rec_rate[rec_rate > 20] <- 20
l_d_s <- log2(1+d_s)
l_rec_rate <- log2(1+rec_rate)

cor.test(l_d_s, l_rec_rate, use="pairwise.complete.obs", method="spearman")

# Make a plot
pdf(file="log_dSNPs_per_sSNP_cMMb.pdf", 6, 6)
par(mar=c(3, 3, 0.25, 0.25), mgp=c(2, 1, 0))
plot(
    l_d_s ~ l_rec_rate,
    xlab="log2(cM/Mb Estimate + 1)",
    ylab="log2(dSNPs per sSNP + 1)",
    main="")
legend("topright", c("rho=-0.096", "p=1.25E-5"))
dev.off()

rr_trim <- as.numeric(del_count$V4)
cor.test(log2(1+dSNPs), log2(1+rr_trim), use="pairwise.complete.obs", method="spearman")
# Make a plot of just dSNPs across recombination rate
pdf(file="log_dSNPs_per_cMMb_Window.pdf", 6, 6)
par(mar=c(3, 3, 0.25, 0.25), mgp=c(2, 1, 0))
plot(
    log2(1+dSNPs) ~ log2(1+rr_trim),
    xlab="log2(cM/Mb Estimate + 1)",
    ylab="log2(dSNP Count + 1)",
    main="")
dev.off()
