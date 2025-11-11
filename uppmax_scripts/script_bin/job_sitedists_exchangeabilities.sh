#!/bin/bash
#SBATCH -A naiss2024-5-197
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00 
#SBATCH --mail-user mahwash.jamy@slu.se
#SBATCH --mail-type=FAIL
#SBATCH -J pb_site_exchangeabilities
#SBATCH -o /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.out
#SBATCH -e /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.err

# conda env
conda activate ptMAGs

# modules
module load bioinfo-tools
module load phylobayesmpi/1.8

readpb_mpi -ss -x 1800 1 $1
python /home/mahja/beta-Cyclocitral/src/convert-site-dists.py "$1".siteprofiles

readpb_mpi -rr -x 1800 1 $1
python /home/mahja/beta-Cyclocitral/src/convert-exchangeabilities.py "$1".meanrr