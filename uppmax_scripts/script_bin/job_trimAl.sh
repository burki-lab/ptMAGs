#!/bin/bash
#SBATCH -A uppmax2025-2-76
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00-02:00:00
#SBATCH --mail-user mahwash.jamy@slu.se
#SBATCH --mail-type=NONE
#SBATCH -J job_trimAl
#SBATCH -o /crex/proj/naiss2023-6-81/Mahwash/beta-Cyclocitral/uppmax_scripts/out_bin/slurm-%A.out
#SBATCH -e /crex/proj/naiss2023-6-81/Mahwash/beta-Cyclocitral/uppmax_scripts/out_bin/slurm-%A.err

# modules
module load bioinfo-tools
module load trimAl

# commands
## command used for inferring reference only trees
##/sw/bioinfo/trimAl/1.4.1/rackham/trimal -in $1 -out $2 -gt 0.8 -fasta

## command for mags plus references - SGTs
/sw/bioinfo/trimAl/1.4.1/rackham/trimal -in $1 -out $2 -gt 0.1 -fasta

## command for mags plus references - concatenated phylogeny
#/sw/bioinfo/trimAl/1.4.1/rackham/trimal -in $1 -out $2 -gt 0.2 -fasta
