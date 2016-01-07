# Deleterious_GP/Data
Small data files for Deleterious Mutations 2. Contains the following files:
- Capture_Sample_Name_Translation.txt
    - Translation between the 'Peter_XY' name format and the actual sample name.
- Covered_Regions.bed.bz2
    - BED file describing regions that have moderate to high sequence depth in
    Deleterious Mutations 1. Created by parsing
    [bedtools](http://bedtools.readthedocs.org/en/latest/) 'genomecov' output
    and identifying regions with at least 20 reads, and at least 200bp long.
    This should be used for estimating sequence-based population genetic
    descriptive statistics.
- Phenotypic_Data/
    - Contains yield and DON (FHB proxy) best linear unbiased predictions (BLUPs),
    and best linear unbiased estimations (BLUEs) as well as 5-location raw data.
- Population_Schematic.pdf
    - A sketch of the history of the population being studied. Hand-drawn by
    Tyler Tiede. May be re-drawn in OmniGraffle eventually.
