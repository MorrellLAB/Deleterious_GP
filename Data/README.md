# Deleterious_GP/Data
Small data files for Deleterious Mutations 2. Contains the following files:
- Capture_Sample_Name_Translation.txt
	- Translation between the 'Peter_XY' name format and the actual sample name.
- Covered_Regions.bed.bz2
	- BED file describing regions that have moderate to high sequence depth in
	Deleterious Mutations 1. Created by parsing bedtools output and identifying
	regions with at least 20 reads, and at least 200bp long.
