[![DOI](https://zenodo.org/badge/1091757605.svg)](https://doi.org/10.5281/zenodo.17634564)

# ptMAGs
Code associated with the manuscript "Identification of a deep-branching lineage of algae using environmental plastid genomes" (Jamy et al 2025).

Preprint here: https://www.biorxiv.org/content/10.1101/2025.01.16.633336v2

## Directory structure
The respository contains a bunch of Jupyter notebooks in different folders. Where relevant, the input files for scripts used to generate plots are provided (along with the paths of where to find them).

1. [00_plastid_MAGS](https://github.com/burki-lab/ptMAGs/tree/main/00_plastid_MAGS)
Here, we annotated the plastid MAGs and calculated stats about their completeness, size, etc. 

2. [01_phylogenies](https://github.com/burki-lab/ptMAGs/tree/main/01_phylogenies)
Here, we ran all phylogenies in the manuscript except for mitochondrial and psbO phylogeneies (those can be found in `04_search_lepto_host`)

3. [02_biogeography](https://github.com/burki-lab/ptMAGs/tree/main/02_biogeography)
No Jupyter notebook here, but we generated plots for leptophyte distribution. 

4. [03_comp_genomics](https://github.com/burki-lab/ptMAGs/tree/main/03_comp_genomics)
Here, we compared gene content across red plastid genomes and looked at synteny across the leptophyte plastid genomes.

5. [src](https://github.com/burki-lab/ptMAGs/tree/main/src)
All the little python scripts and R scripts used at various points. 

6. [uppmax_scripts](https://github.com/burki-lab/ptMAGs/tree/main/uppmax_scripts/script_bin)
All the Slurm scripts used to submit jobs to our cluster (uppmax). Have a look here to see the exact commands used for running phylogenies, and more.

7. [PATHS.json](https://github.com/burki-lab/ptMAGs/blob/main/PATHS.json)
A json file specifying where all files were located.

8. [SLURMSCRIPTS.json](https://github.com/burki-lab/ptMAGs/blob/main/SLURMSCRIPTS.json)
A json file pointing to some Slurm scripts that were used regularly. 

The respoistory is slightly messy (might be an understatement), but hopefully contains (almost) everything needed to reproduce the analyses. 

