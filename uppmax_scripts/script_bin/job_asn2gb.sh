#!/bin/bash
#SBATCH -A uppmax2025-2-76
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH --mail-user mahwash.jamy@slu.se
#SBATCH --mail-type=NONE
#SBATCH -J job_2023_10_12_asn2gb
#SBATCH -o /crex/proj/naiss2023-6-81/Mahwash/beta-Cyclocitral/uppmax_scripts/out_bin/slurm-%A.out
#SBATCH -e /crex/proj/naiss2023-6-81/Mahwash/beta-Cyclocitral/uppmax_scripts/out_bin/slurm-%A.err

# commands
## asn2gb downloaded from https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/asn2gb/
/home/mahja/beta-Cyclocitral/src/asn2gb -i $1 -o $2
