#!/bin/bash
#SBATCH -A naiss2024-5-197
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 06:00:00
#SBATCH --mail-user mahwash.jamy@slu.se
#SBATCH --mail-type=FAIL
#SBATCH -J job_2023_10_17_RAxML
#SBATCH -o /crex/proj/naiss2023-6-81/Mahwash/beta-Cyclocitral/uppmax_scripts/out_bin/slurm-%A.out
#SBATCH -e /crex/proj/naiss2023-6-81/Mahwash/beta-Cyclocitral/uppmax_scripts/out_bin/slurm-%A.err

# modules
module load bioinfo-tools

# commands (using raxml-ng v1.2.0)
~/bin/raxml-ng --all --msa $1 --model GTR+G --prefix $2 --seed $RANDOM --threads 3 --bs-trees 100 --bs-metric fbp

