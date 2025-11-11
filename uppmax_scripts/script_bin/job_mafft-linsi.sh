#!/bin/bash
#SBATCH -A naiss2024-5-197
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 00-06:00:00
#SBATCH --mail-user mahwash.jamy@slu.se
#SBATCH --mail-type=NONE
#SBATCH -J job_MAFFT_linsi
#SBATCH -o /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.out
#SBATCH -e /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.err

# modules
module load bioinfo-tools
module load MAFFT

# commands
/sw/apps/bioinfo/MAFFT/7.407/rackham/bin/mafft-linsi --thread 6 --reorder --adjustdirection $1 > $2
