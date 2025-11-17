#!/bin/bash
#SBATCH -A naiss2024-5-197
#SBATCH -M snowy
#SBATCH -p node
#SBATCH -n 16
#SBATCH -C mem512GB
#SBATCH -t 7-10:00:00 
#SBATCH --mail-user mahwash.jamy@slu.se
#SBATCH --mail-type=FAIL
#SBATCH -J megahit
#SBATCH -o /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.out
#SBATCH -e /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.err

## Load modules
module load bioinfo-tools
module load megahit/1.2.9

## Assemble!
megahit -1 $1 -2 $2 --min-contig-len $3 -o $4 -t 16
