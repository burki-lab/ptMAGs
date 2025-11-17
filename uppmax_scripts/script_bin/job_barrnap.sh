#!/bin/bash
#SBATCH -A naiss2024-5-197
#SBATCH -M rackham
#SBATCH -p node
#SBATCH -n 8
#SBATCH -t 2:00:00 
#SBATCH --mail-user mahwash.jamy@slu.se
#SBATCH --mail-type=FAIL
#SBATCH -J barrnap
#SBATCH -o /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.out
#SBATCH -e /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.err

## Load modules
module load bioinfo-tools
module load barrnap/0.9

## Run command
barrnap --kingdom euk --threads 8 --lencutoff 0.1 --outseq $2 $1
