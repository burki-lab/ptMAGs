#!/bin/bash
#SBATCH -A naiss2024-5-197
#SBATCH -M snowy
#SBATCH -p node
#SBATCH -n 16
#SBATCH -C mem256GB
#SBATCH -t 2-0:00:00 
#SBATCH --mail-user mahwash.jamy@slu.se
#SBATCH --mail-type=FAIL
#SBATCH -J job_ghost
#SBATCH -o /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.out
#SBATCH -e /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.err

# modules
module load bioinfo-tools
##module load iqtree/2.2.2.6-omp-mpi
module load iqtree/1.6.12-omp-mpi


## LG+MEOW(60,20)+H4 model
##/home/mahja/bin/iqtree-2.3.6-Linux-intel/bin/iqtree2 -s $1 --prefix $2 -T 16 -m LG+MEOW60+H4 -mdef $3 -mwopt -te $4 -optlen BFGS

## LG+MEOW(60,20)+H6 model
##/home/mahja/bin/iqtree-2.3.6-Linux-intel/bin/iqtree2 -s $1 --prefix $2 -T 16 -m LG+MEOW60+H6 -mdef $3 -te $4 -optlen BFGS 

## LG+MEOW(60,20)+H8 model
##/home/mahja/bin/iqtree-2.3.6-Linux-intel/bin/iqtree2 -s $1 --prefix $2 -T 16 -m LG+MEOW60+H8 -mdef $3 -mwopt -te $4 -optlen BFGS 

## model test for number of rate classes in GHOST model
iqtree-omp -s $1 -pre $2 -nt 16 -m TESTONLY -mdef $3 -mset LG -madd LG+MEOW60+H4,LG+MEOW60+H6,LG+MEOW60+H8 -mwopt -optlen BFGS

