#!/bin/bash
#SBATCH -A uppmax2025-2-76
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5-00:00:00
#SBATCH --mail-user mahwash.jamy@slu.se
#SBATCH --mail-type=NONE
#SBATCH -J job_bmge
#SBATCH -o /crex/proj/naiss2023-6-81/Mahwash/beta-Cyclocitral/uppmax_scripts/out_bin/slurm-%A.out
#SBATCH -e /crex/proj/naiss2023-6-81/Mahwash/beta-Cyclocitral/uppmax_scripts/out_bin/slurm-%A.err

# modules
module load bioinfo-tools

# commands
## Thomas's command
##/crex/proj/naiss2023-6-81/conda/envs/ptMAGs/bin/bmge -g 0.2 -b 10 -t AA -m BLOSUM35 -i $1 -of $2

## standard command used throughout
##/crex/proj/naiss2023-6-81/conda/envs/ptMAGs/bin/bmge -g 0.8 -t AA -m BLOSUM35 -i $1 -of $2

## standard command adapted for mitochondria
/crex/proj/naiss2023-6-81/conda/envs/ptMAGs/bin/bmge -g 0.8 -t AA -m BLOSUM30 -i $1 -of $2

## for testing. Gap threshold = 0.6
##/crex/proj/naiss2023-6-81/conda/envs/ptMAGs/bin/bmge -g 0.6 -t AA -m BLOSUM35 -i $1 -of $2

## for testing. Gap threshold = 0.4
##/crex/proj/naiss2023-6-81/conda/envs/ptMAGs/bin/bmge -g 0.4 -t AA -m BLOSUM35 -i $1 -of $2

## for testing. Gap threshold = 0.2
##/crex/proj/naiss2023-6-81/conda/envs/ptMAGs/bin/bmge -g 0.2 -t AA -m BLOSUM35 -i $1 -of $2

## stationary-based trimming to remove compositionally heterogenous sites
##/crex/proj/naiss2023-6-81/conda/envs/ptMAGs/bin/bmge -t AA -s FAST -h 0:1 -g 1 -i $1 -of $2