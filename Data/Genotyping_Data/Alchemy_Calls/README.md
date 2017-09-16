# ALCHEMY Genotype Calling
There were issues with the genotyping data downloaded from T3: parental A/B
allele calls were incompatible with the observed genotypes segregating in their
progeny families. We downloaded raw intensity data for the BOPA experiments
from T3, and from S. Chao at a USDA lab in Fargo, ND.

## Workflow
1. Download the data for the experiments. The parental lines are included on the
   genotyping runs listed in this table. Use both BOPA1 and BOPA2 chips:

   | Chip    | Parent (Pedigree Name)    | Synonym (Genotyping Name) |
   |---------|---------------------------|---------------------------|
   | 2008 MN | FEG183-52                 | 08MN-32                   |
   | 2008 BA | 6B06-1132                 | 08BA-45                   |
   | 2008 ND | ND20448                   | 08N6-01                   |
   | 2008 ND | ND24906                   | 08N6-04                   |
   | 2008 ND | ND25160                   | 08N6-13                   |
   | 2008 ND | ND25986                   | 08N6-47                   |
   | 2008 ND | ND26036                   | 08N6-62                   |
   | 2008 ND | ND26104                   | 08N6-73                   |
   | 2007 MN | FEG153-58                 | 07-MN-04                  |
   | 2007 MN | FEG154-47                 | 07-MN-06                  |
   | 2007 MN | FEG175-57                 | 07-MN-52                  |
   | 2007 BA | 6B04-0290                 | 07-BA-11                  |
   | 2007 BA | 6B05-0922                 | 07-BA-39                  |
   | 2007 ND | ND25652                   | 07-N6-12                  |
   | 2007 ND | ND25728                   | 07-N6-88                  |
   | 2006 MN | FEG141-20                 | FEG141-20                 |
   | 2006 BA | 6B01-2218                 | 6B01-2218                 |
   | 2006 BA | 6B03-4304                 | 6B03-4304                 |
   | 2006 BA | 6B03-4478                 | 6B03-4478                 |

2. Import the data into GenomeStudio version 3. You must place the manifest,
   sample sheet, and intensity files into the same directory, and import them
   as a new experiment in GenomeStudio.
3. Install the ALCHEMY GenomeStudio plugin.
4. Export the ALCHEMY input files as a "custom report" export from the
   GenomeStudio menu. Export a custom report for each experiment, and place them
   into the same directory.
5. Edit `Create_Alchemy_Inputs.sh` and enter the path to the custom export
   report directory, A-B allele translation CSV, and the Python script that
   generates ALCHEMY SNP map and sample map files. Run the script.
6. Edit `Alchemy_Calls.sh` and enter the paths to the custom GenomeStudio
   reports directory, ALCHEMY SNP map and sample map directory, and the desired
   output directory. Run the script. This should be run on MSI, where a module
   for ALCHEMY is available.
7. Repeat steps 2-7 for the progeny genotyping data. Be sure to remove the
   problematic samples from Cycle 2 genotyping data: the low signal values for
   both channels stalls the ALCHEMY program.
8. Edit the paths in `Batch_Alchemy_to_PLINK.sh` and run the script. It will
   combine BOPA1 and BOPA2 data into a single file for each program. Then, it
   will sort the maps, interpolate genetic and physical positions, then
   reorder the PED to match the MAPs. Finally, it will trim the combined data
   to just the parents of the population, and maek the names consistent with
   the rest of the datasets.

When trimming down genotype data for the parental lines, use the line names
listed in the `Synonym` column in the above table.
