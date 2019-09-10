# This is a simple script to randomly mask a given proportion of a phenotype
# vector. This was built for the purpose of cross-validation of GEMMA results.
# Takes the PLINK fam file (with phenotype) as an argument

# Mask 20% of the vector
mask_prop <- 0.2

# Take arguments
args <- commandArgs(TRUE)
phen_dat <- args[1]

# Read the table
phen <- read.table(phen_dat, header=FALSE, sep=" ")

# Find out the number to mask, rounding UP to the next integer
nmask <- ceiling(nrow(phen) * mask_prop)

# Get sample the IDs to mask via sampling
mask <- sample(1:nrow(phen), size=nmask, replace=FALSE)

# Mask them by setting to NA (missing)
phen$V6[mask] <- NA

# Cast to matrix
phen <- as.matrix(phen)

# And write it to stdout
write.table(phen, file="", quote=F, sep=" ", col.names=FALSE, row.names=FALSE)
