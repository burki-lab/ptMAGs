#!/bin/bash
#SBATCH -A naiss2024-5-197
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 2:00:00 
#SBATCH --mail-user mahwash.jamy@slu.se
#SBATCH --mail-type=FAIL
#SBATCH -J job_gfmix
#SBATCH -o /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.out
#SBATCH -e /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.err

# Load modules
module load R/4.2.1
module load R_packages/4.2.1

# set options
phylip=$1
tree=$2
iqtree=$3
freq=$4
root=$5

# Echo input tree file
echo "$tree"

## Run gfmix (provided in the `src` folder in the repository)
/home/mahja/bin/garp/gfmix_custombins -s "$phylip" -t "$tree" -d -i "$iqtree" -f "$freq" -r "$root" -gclass INFKTDSE -fclass RCVWHMAQ
  
