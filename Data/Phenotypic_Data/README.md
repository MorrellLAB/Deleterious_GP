# Data/Phenotypic_Data
This directory contains the phenotypic data used in the genomic prediction
experiment. The data in here has been spatially adjusted by Tyler Tiede, in
October 2016. He used the following parameters to adjust the data:

- Outliers removed based on studentized residuals with alpha = 0.01
- Spatial variation for yield and DON accounted for by a moving average covariate
- Spatial variation for height accounted for with with random effects for row, column, row_blk, and column_blk (all nested within trials)
- Only significant effects (tested by p.value for fixed effects and LRT for random effects) within trials fit in the multi-environment analysis
- Significant line_name:trial effects for yield and DON, not height.... variation accounted for by random term in model

Scripts/package to perform spatial correction are available by request from
Tyler.

A note on line names:
The line names have been modified from their original sheet so that they are
consistent between phenotypic data files. The following changes have been made:

- G10W was replaced with MS10S3
- MS11S3 was replaced with MS11S2
- MS10S4 was replaced with MS10S3
- Line names ending in '-dd' were replaced with '-0dd', where d is any integer.
- M135(FEG141-20) was replaced with FEG141-20
- M138(FEG153-25) was replaced with M138
- M139(FEG153-58) was replaced with FEG153-58
- M141(FEG175-57) was replaced with FEG175-57
- M136(FEG147-37) was replaced with M136
- M140(FEG154-47) was replaced with FEG154-47
- Check variety names were converted to Title Case from ALLCAPS.
- A dummy row of NA was added for Quest, which does not appear in the adjusted phenotypic data file.
- Dummy row added for FILLER in adjusted file
- Dummy row added for MS12_2158-023 in adjusted file
- Dummy row added for MS12_2108-006 in adjusted file

Raw phenotypic data files without the line name modifications above are
available upon request.
