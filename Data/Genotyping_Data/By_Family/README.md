# Deleterious_GP/Data/Genotyping_Data/By_Family
Genotyping data for Deleterious Mutations 2. 

Contains PLINK PED/MAP files for each family in each cycle. The PED files here
are different than the ones present in the unified PED file in that each family
has been independently cleaned of strand flip errors. The following procedure
was taken (e.g., for family MS10S3001 in Cycle 1):

```
grep "$(cut -f 3 MS10S3001.ped | uniq)\b" C0.ped >> Parents.ped
grep "$(cut -f 4 MS10S3001.ped | uniq)\b" C0.ped >> Parents.ped
plink2 --file Parents --make-bed --allow-extra-chr --out Parents_Bin
plink2 --file MS10S3001 --make-bed --allow-extra-chr --out Prog_Bin
plink2 --bfile Parents_Bin --bmerge Prog_Bin --allow-extra-chr --out Merged
plink2 --bfile Prog_Bin --flip Merged.missnp --make-bed --out Prog_Bin_Flipped --allow-extra-chr
plink2 --bfile Parents_Bin --bmerge Prog_Bin_Flipped --out Merged_Flipped --allow-extra-chr
plink2 --bfile Merged_Flipped --recode --out MS10S3001_Merged --allow-extra-chr
python Fix_Map_Order.py MS10S3001.map MS10S3001_Merged.map MS10S3001_Merged.ped > MS10S3001_Reordered.ped
```

Use PLINK2 to merge the parents with the progeny, and check for potential flip
errors. This is done by identifying triallelic sites. Then, flip the problematic
SNPs in the progeny, and re-merge. PLINK2 re-orders the markers based on
physical position, so we have to re-order based on genetic position.

PLINK2 merge doesn't capture all flip errors - G/C and A/T SNPs are still
potentially flipped. We can't really solve this problem, though, without any
external information. What we can do, though, is check for obviously wrong
SNPs. If both parents are homozygous for G, and the progeny all genotype as C,
we can flip those. For the sake of doing actual anaylsis, the monomorphic SNPs
are noninformative, but we keep them for convenience, since we can then
use the same MAP file for each family.
