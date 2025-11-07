#!/bin/bash
#SBATCH -A uppmax2025-2-76
#SBATCH -M rackham
#SBATCH -p node
#SBATCH -n 16
#SBATCH -t 3-00:00:00 
#SBATCH --mail-user mahwash.jamy@slu.se
#SBATCH --mail-type=FAIL
#SBATCH -J job_IQTREE
#SBATCH -o /crex/proj/naiss2023-6-81/Mahwash/beta-Cyclocitral/uppmax_scripts/out_bin/slurm-%A.out
#SBATCH -e /crex/proj/naiss2023-6-81/Mahwash/beta-Cyclocitral/uppmax_scripts/out_bin/slurm-%A.err

# modules
module load bioinfo-tools
module load iqtree/2.2.2.6-omp-mpi
##module load iqtree/1.6.12-omp-mpi


# commands
## For SGTs - single core
##iqtree2 -s $1 --prefix $2 -B 1000 -T 1 -bnni -m MFP

## LG4X model
##iqtree2 -s $1 --prefix $2 -B 1000 -T 16 -bnni -m LG4X

## For concatenated tree with full dataset
##iqtree2 -s $1 --prefix $2 -B 1000 -T 16 -bnni -m MFP -mset LG,cpREV -mrate E,I,G,I+G

## For concatenated tree with subset of the data
##iqtree2 -s $1 --prefix $2 -B 1000 -T 16 -bnni -mset LG,cpREV -madd LG+C60,LG+C60+G,cpREV+C60,cpREV+C60+G
### The command above did not test the C60 models, so we tried with with iqtree v1
##iqtree-omp -s $1 -pre $2 -bb 1000 -nt 16 -bnni -m TEST -mset LG,cpREV -madd LG+C60,LG+C60+G,cpREV+C60,cpREV+C60+G

## cpREV+C60+G
##iqtree2 -s $1 --prefix $2 -B 1000 -T 16 -bnni -m cpREV+C60+G -mwopt

## LG+C60+G
iqtree2 -s $1 --prefix $2 -B 1000 -T 16 -bnni -m LG+C60+G -mwopt

## Model test only
##iqtree-omp -s $1 -pre $2 -bb 1000 -nt 8 -bnni -m TESTONLY -mset LG,cpREV -madd LG+C60,LG+C60+G,cpREV+C60,cpREV+C60+G

## EX_EHO + G4
##iqtree2 -s $1 --prefix $2 -B 1000 -T 16 -bnni -m EX_EHO+G4

## cpREV+C60+G+F
##iqtree2 -s $1 --prefix $2 -B 1000 -T 16 -bnni -m cpREV+C60+G+F

## Poisson+C60+G
##iqtree2 -s $1 --prefix $2 -B 1000 -T 16 -bnni -m Poisson+C60+G

## cpREV+C60+R9
##iqtree2 -s $1 --prefix $2 -B 1000 -T 16 -bnni -m cpREV+C60+R9

## model test cpREV+MAM60+G vs cpREV+C60+G
##iqtree-omp -s $1 -pre $2 -nt 16 -m TESTONLY -mdef $3 -mset LG,WAG,cpREV -madd cpREV+C60+G,cpREV+ESmodel+G,LG+C60+G,LG+ESmodel+G,WAG+C60+G,WAG+ESmodel+G -mwopt

## cpREV+MAM60+G model
##iqtree2 -s $1 --prefix $2 -B 1000 -T 16 -bnni -m cpREV+ESmodel+G -mdef $3 

## LG+MAM60+G model
##iqtree2 -s $1 --prefix $2 -B 1000 -T 16 -bnni -m LG+ESmodel+G -mdef $3 

## LG+MEOW(60,20)+G model
##iqtree2 -s $1 --prefix $2 -B 1000 -T 16 -bnni -m LG+MEOW60+G -mdef $3 -mwopt 

## LG+MEOW(60,20)+G model - constrained tree search
##iqtree2 --runs 2 -s $1 --prefix $2 -T 16 -m LG+MEOW60+G -mdef $3 -g $4 -mwopt 

## LG+C60 model for sr4 recoded data - constrained
##iqtree2 --runs 5 -s $1 --prefix $2 -T 16 -m xmC60SR4 -mdef $3 -g $4 -mwopt 

## cpREV-PMSF
##iqtree2 -s $1 --prefix $2 -B 1000 -T 16 -bnni -m cpREV+G4 -fs $3

## CATPMSF - step 1
## Run tree with cpREV+G model or best fit model
##iqtree2 -s $1 --prefix $2 -T 16 -m cpREV+G
##iqtree2 -s $1 --prefix $2 -T 16 -m LG+G

## CATPMSF - step 3 - 5000 ufb
##/home/mahja/bin/iqtree-2.3.6-Linux-intel/bin/iqtree2 -s $1 --prefix $2 -B 5000 -T 16 -bnni -mdef $3 -m cat+G4 -fs $4
##/home/mahja/bin/iqtree-2.3.6-Linux-intel/bin/iqtree2 -s $1 --prefix $2 -B 1000 -T 16 -bnni -m LG+G4 -fs $3 -mwopt

## CATPMSF - step 3 - 100 fbp
##/home/mahja/bin/iqtree-2.3.6-Linux-intel/bin/iqtree2 -s $1 --prefix $2 -b 100 -T 16 -mdef $3 -m cat+G4 -fs $4 -mwopt

## CATPMSF model - constrained tree search
##/home/mahja/bin/iqtree-2.3.6-Linux-intel/bin/iqtree2 --runs 5 -s $1 --prefix $2 -T 16 -mdef $3 -m cat+G4 -fs $4 -g $5 -mwopt

## CHOOSING THE BEST MEOW MODEL 
## BIC CRITERIA
##iqtree-omp -s $1 -pre $2 -nt 16 -m TESTONLY -mdef $3 -mset LG,cpREV -madd cpREV+MEOW40+G,cpREV+MEOW60+G,cpREV+MEOW80+G,LG+MEOW40+G,LG+MEOW60+G,LG+MEOW80+G -mwopt


## 2-FOLD CORRECTED CROSS VALIDATION
### Calculate maximized log-likelihood for the whole alignment, X and tree, T and adjustable parameters, theta.
##iqtree2 -s $1 --prefix $2 -T 16 -m cpREV+ESmodel+G -mdef $3 -mwopt


