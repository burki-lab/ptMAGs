#!/bin/bash
#SBATCH -A naiss2024-5-197
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00-01:00:00
#SBATCH --mail-user mahwash.jamy@slu.se
#SBATCH --mail-type=NONE
#SBATCH -J job_trimAl_16S
#SBATCH -o /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.out
#SBATCH -e /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.err

# modules
module load bioinfo-tools
module load trimAl

# commands
## old 16S 
#/sw/bioinfo/trimAl/1.4.1/rackham/trimal -in $1 -out $2 -gt 0.3 -st 0.001

/sw/bioinfo/trimAl/1.4.1/rackham/trimal -in $1 -out $2 -gt 0.1