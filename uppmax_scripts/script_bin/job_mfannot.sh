#!/bin/bash
#SBATCH -A uppmax2025-2-76
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 2-00:00:00
#SBATCH --mail-user mahwash.jamy@slu.se
#SBATCH --mail-type=NONE
#SBATCH -J job_mfannot
#SBATCH -o /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.out
#SBATCH -e /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.err
# modules
module load bioinfo-tools

# commands
MFANNOT_LIB_PATH=$3
## for plastids
##singularity run /crex/proj/naiss2023-6-81/Thomas/data/analysis_plastomes/mfannot/mfannot_latest.sif /mfannot/mfannot --partial $1 -o $2 --sqn -g 11

## for mitochondria (code 1 or 4 as appropriate)
singularity run /crex/proj/naiss2023-6-81/Thomas/data/analysis_plastomes/mfannot/mfannot_latest.sif /mfannot/mfannot --partial $1 -o $2 --sqn -g 4
