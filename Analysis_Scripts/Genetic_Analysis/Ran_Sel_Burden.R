# Make a table of the mean/sd homozygous genotypes by class
dat <- read.csv("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Burden/Ran_Sel_Homozygous_Derived_Counts.csv", header=TRUE)

c0.nc.mean <- mean(dat$NC[dat$Cycle == "C0"])
c0.nc.sd <- sd(dat$NC[dat$Cycle == "C0"])
c0.sy.mean <- mean(dat$SY[dat$Cycle == "C0"])
c0.sy.sd <- sd(dat$SY[dat$Cycle == "C0"])
c0.ns.mean <- mean(dat$NS[dat$Cycle == "C0"])
c0.ns.sd <- sd(dat$NS[dat$Cycle == "C0"])
c0.de.mean <- mean(dat$DE[dat$Cycle == "C0"])
c0.de.sd <- sd(dat$DE[dat$Cycle == "C0"])

c1.ran.nc.mean <- mean(dat$NC[dat$Cycle == "C1" & dat$Type == "Ran"])
c1.ran.nc.sd <- sd(dat$NC[dat$Cycle == "C1" & dat$Type == "Ran"])
c1.ran.sy.mean <- mean(dat$SY[dat$Cycle == "C1" & dat$Type == "Ran"])
c1.ran.sy.sd <- sd(dat$SY[dat$Cycle == "C1" & dat$Type == "Ran"])
c1.ran.ns.mean <- mean(dat$NS[dat$Cycle == "C1" & dat$Type == "Ran"])
c1.ran.ns.sd <- sd(dat$NS[dat$Cycle == "C1" & dat$Type == "Ran"])
c1.ran.de.mean <- mean(dat$DE[dat$Cycle == "C1" & dat$Type == "Ran"])
c1.ran.de.sd <- sd(dat$DE[dat$Cycle == "C1" & dat$Type == "Ran"])

c1.sel.nc.mean <- mean(dat$NC[dat$Cycle == "C1" & dat$Type == "Sel"])
c1.sel.nc.sd <- sd(dat$NC[dat$Cycle == "C1" & dat$Type == "Sel"])
c1.sel.sy.mean <- mean(dat$SY[dat$Cycle == "C1" & dat$Type == "Sel"])
c1.sel.sy.sd <- sd(dat$SY[dat$Cycle == "C1" & dat$Type == "Sel"])
c1.sel.ns.mean <- mean(dat$NS[dat$Cycle == "C1" & dat$Type == "Sel"])
c1.sel.ns.sd <- sd(dat$NS[dat$Cycle == "C1" & dat$Type == "Sel"])
c1.sel.de.mean <- mean(dat$DE[dat$Cycle == "C1" & dat$Type == "Sel"])
c1.sel.de.sd <- sd(dat$DE[dat$Cycle == "C1" & dat$Type == "Sel"])

c2.ran.nc.mean <- mean(dat$NC[dat$Cycle == "C2" & dat$Type == "Ran"])
c2.ran.nc.sd <- sd(dat$NC[dat$Cycle == "C2" & dat$Type == "Ran"])
c2.ran.sy.mean <- mean(dat$SY[dat$Cycle == "C2" & dat$Type == "Ran"])
c2.ran.sy.sd <- sd(dat$SY[dat$Cycle == "C2" & dat$Type == "Ran"])
c2.ran.ns.mean <- mean(dat$NS[dat$Cycle == "C2" & dat$Type == "Ran"])
c2.ran.ns.sd <- sd(dat$NS[dat$Cycle == "C2" & dat$Type == "Ran"])
c2.ran.de.mean <- mean(dat$DE[dat$Cycle == "C2" & dat$Type == "Ran"])
c2.ran.de.sd <- sd(dat$DE[dat$Cycle == "C2" & dat$Type == "Ran"])

c2.sel.nc.mean <- mean(dat$NC[dat$Cycle == "C2" & dat$Type == "Sel"])
c2.sel.nc.sd <- sd(dat$NC[dat$Cycle == "C2" & dat$Type == "Sel"])
c2.sel.sy.mean <- mean(dat$SY[dat$Cycle == "C2" & dat$Type == "Sel"])
c2.sel.sy.sd <- sd(dat$SY[dat$Cycle == "C2" & dat$Type == "Sel"])
c2.sel.ns.mean <- mean(dat$NS[dat$Cycle == "C2" & dat$Type == "Sel"])
c2.sel.ns.sd <- sd(dat$NS[dat$Cycle == "C2" & dat$Type == "Sel"])
c2.sel.de.mean <- mean(dat$DE[dat$Cycle == "C2" & dat$Type == "Sel"])
c2.sel.de.sd <- sd(dat$DE[dat$Cycle == "C2" & dat$Type == "Sel"])

c3.ran.nc.mean <- mean(dat$NC[dat$Cycle == "C3" & dat$Type == "Ran"])
c3.ran.nc.sd <- sd(dat$NC[dat$Cycle == "C3" & dat$Type == "Ran"])
c3.ran.sy.mean <- mean(dat$SY[dat$Cycle == "C3" & dat$Type == "Ran"])
c3.ran.sy.sd <- sd(dat$SY[dat$Cycle == "C3" & dat$Type == "Ran"])
c3.ran.ns.mean <- mean(dat$NS[dat$Cycle == "C3" & dat$Type == "Ran"])
c3.ran.ns.sd <- sd(dat$NS[dat$Cycle == "C3" & dat$Type == "Ran"])
c3.ran.de.mean <- mean(dat$DE[dat$Cycle == "C3" & dat$Type == "Ran"])
c3.ran.de.sd <- sd(dat$DE[dat$Cycle == "C3" & dat$Type == "Ran"])

c3.sel.nc.mean <- mean(dat$NC[dat$Cycle == "C3" & dat$Type == "Sel"])
c3.sel.nc.sd <- sd(dat$NC[dat$Cycle == "C3" & dat$Type == "Sel"])
c3.sel.sy.mean <- mean(dat$SY[dat$Cycle == "C3" & dat$Type == "Sel"])
c3.sel.sy.sd <- sd(dat$SY[dat$Cycle == "C3" & dat$Type == "Sel"])
c3.sel.ns.mean <- mean(dat$NS[dat$Cycle == "C3" & dat$Type == "Sel"])
c3.sel.ns.sd <- sd(dat$NS[dat$Cycle == "C3" & dat$Type == "Sel"])
c3.sel.de.mean <- mean(dat$DE[dat$Cycle == "C3" & dat$Type == "Sel"])
c3.sel.de.sd <- sd(dat$DE[dat$Cycle == "C3" & dat$Type == "Sel"])

# Make a table to print it out
toprint <- data.frame(
    C0.Mean=c(c0.nc.mean, c0.sy.mean, c0.ns.mean, c0.de.mean),
    C0.Sd=c(c0.nc.sd, c0.sy.sd, c0.ns.sd, c0.de.sd),
    C1.Ran.Mean=c(c1.ran.nc.mean, c1.ran.sy.mean, c1.ran.ns.mean, c1.ran.de.mean),
    C1.Ran.Sd=c(c1.ran.nc.sd, c1.ran.sy.sd, c1.ran.ns.sd, c1.ran.de.sd),
    C1.Sel.Mean=c(c1.sel.nc.mean, c1.sel.sy.mean, c1.sel.ns.mean, c1.sel.de.mean),
    C1.Sel.Sd=c(c1.sel.nc.sd, c1.sel.sy.sd, c1.sel.ns.sd, c1.sel.de.sd),
    C2.Ran.Mean=c(c2.ran.nc.mean, c2.ran.sy.mean, c2.ran.ns.mean, c2.ran.de.mean),
    C2.Ran.Sd=c(c2.ran.nc.sd, c2.ran.sy.sd, c2.ran.ns.sd, c2.ran.de.sd),
    C2.Sel.Mean=c(c2.sel.nc.mean, c2.sel.sy.mean, c2.sel.ns.mean, c2.sel.de.mean),
    C2.Sel.Sd=c(c2.sel.nc.sd, c2.sel.sy.sd, c2.sel.ns.sd, c2.sel.de.sd),
    C3.Ran.Mean=c(c3.ran.nc.mean, c3.ran.sy.mean, c3.ran.ns.mean, c3.ran.de.mean),
    C3.Ran.Sd=c(c3.ran.nc.sd, c3.ran.sy.sd, c3.ran.ns.sd, c3.ran.de.sd),
    C3.Sel.Mean=c(c3.sel.nc.mean, c3.sel.sy.mean, c3.sel.ns.mean, c3.sel.de.mean),
    C3.Sel.Sd=c(c3.sel.nc.sd, c3.sel.sy.sd, c3.sel.ns.sd, c3.sel.de.sd)
    )
rownames(toprint) <- c("Noncoding", "Synonymous", "Nonsynonymous", "Deleterious")
write.csv(toprint, file="Ran_Sel_Homozygous_Derived_Counts.csv", quote=FALSE, row.names=TRUE)
