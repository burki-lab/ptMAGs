#!/bin/bash
#SBATCH -A uppmax2025-2-76
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 5-00:00:00 
#SBATCH --mail-user mahwash.jamy@slu.se
#SBATCH --mail-type=FAIL
#SBATCH -J job_AU
#SBATCH -o /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.out
#SBATCH -e /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.err

# modules
module load bioinfo-tools
module load iqtree/2.2.2.6-omp-mpi

## cpREV+C60+G
##iqtree2 -s $1 --prefix $2 -T 16 -m cpREV+C60+G -z $3 -n 0 -zb 10000 -au

## lg-meow80-g model 
iqtree2 -s $1 --prefix $2 -T 16 -m LG+MEOW60+G -mdef $3 -mwopt -te $4 -z $5 -zb 10000 -au   

## cat-pmsf model
##/home/mahja/bin/iqtree-2.3.6-Linux-intel/bin/iqtree2 -s $1 --prefix $2 -T 16 -mdef $3 -m cat+G4 -fs $4 -te $5 -z $6 -n 0 -zb 10000 -au

## ghost 8 categories model
##iqtree2 -s $1 --prefix $2 -T 16 -mdef $3 -m LG+MEOW60+H8 -mwopt -te $4 -z $5 -zb 10000 -au -optlen BFGS