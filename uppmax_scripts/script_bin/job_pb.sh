#!/bin/bash
#SBATCH -A uppmax2025-2-76
#SBATCH -M rackham
#SBATCH -p node
#SBATCH -n 16
#SBATCH -t 7-00:00:00 
#SBATCH --mail-user mahwash.jamy@slu.se
#SBATCH --mail-type=FAIL
#SBATCH -J job_pb
#SBATCH -o /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.out
#SBATCH -e /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.err

# modules
module load bioinfo-tools
module load phylobayesmpi/1.8

## For constrained tree for CAT-PMSF. Step 1 of CAT-PMSF
##mpirun -np 14 pb_mpi -cat -gtr -d $1 -T $2 $3

## Continue run
mpirun -np 14 pb_mpi $1

## Regular run
##mpirun -np 14 pb_mpi -cat -gtr -dgam 4 -dc -d $1 $2
