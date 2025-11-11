#!/bin/bash
#SBATCH -A uppmax2025-2-76
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00-04:00:00
#SBATCH --mail-user mahwash.jamy@slu.se
#SBATCH --mail-type=NONE
#SBATCH -J job_prequal
#SBATCH -o /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.out
#SBATCH -e /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.err

# commands
gene=$(basename $1 | cut -f1 -d '.')
/crex/proj/naiss2023-6-81/bin/prequal/prequal -filterthresh 0.95 $1
mv "$1".filtered "$2"/"$gene".fasta
mv "$1".filtered.PP "$2"/"$gene".fasta.filtered.PP
mv "$1".warning "$2"/"$gene".fasta.warning

