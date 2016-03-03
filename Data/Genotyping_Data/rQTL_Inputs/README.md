# Deleterious_GP/Data/Genotyping_Data/rQTL_Inputs
Contains [r/QTL](http://www.rqtl.org/) input files in the "csv" format, for
each cycle. Phenotypic data is not present in these data files yet, and the
phenotype column is filled with missing values.

## To Read into r/QTL
The individuals in each cycle are F3s derived from crosses among individuals
in the previous cycle. r/QTL cannot handle F3 in itself (the default is to
treat the data as F2 data), but you can read the data in as a backcross with
0 generations of backcross and 3 generations of selfing:

```r
cross <- read.cross(
    format="csv",
    file="Par1_x_Par2.csv",
    BC.gen=0,
    F.gen=3)
```

Note that if you read the data in this way, then certain functions will not
work. `checkAlleles()` and similar linkage map construction functions are only
defined for F2 cross objects.
