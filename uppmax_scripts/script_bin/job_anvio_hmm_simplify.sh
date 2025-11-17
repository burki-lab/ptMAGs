#!/bin/bash
#SBATCH -A naiss2024-5-197
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00 
#SBATCH --mail-user mahwash.jamy@slu.se
#SBATCH --mail-type=FAIL
#SBATCH -J simplify
#SBATCH -o /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.out
#SBATCH -e /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.err

# modules
module load conda
source conda_init.sh
conda activate /crex/proj/naiss2023-6-81/conda/envs/anvio-8

## Generate contigs database 
anvi-script-reformat-fasta $1 -o $2 -l 0 --simplify-names
