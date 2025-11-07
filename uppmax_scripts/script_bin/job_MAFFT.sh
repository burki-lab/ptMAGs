#!/bin/bash
#SBATCH -A uppmax2025-2-76
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 00-03:00:00
#SBATCH --mail-user mahwash.jamy@slu.se
#SBATCH --mail-type=NONE
#SBATCH -J job_MAFFT
#SBATCH -o /crex/proj/naiss2023-6-81/Mahwash/beta-Cyclocitral/uppmax_scripts/out_bin/slurm-%A.out
#SBATCH -e /crex/proj/naiss2023-6-81/Mahwash/beta-Cyclocitral/uppmax_scripts/out_bin/slurm-%A.err

# modules
module load bioinfo-tools
module load MAFFT

# commands
## for references only, we used mafft auto
##/sw/apps/bioinfo/MAFFT/7.407/rackham/bin/mafft --auto $3 --thread 10 --reorder --adjustdirection $1 > $2

## mafft-ginsi
/sw/apps/bioinfo/MAFFT/7.407/rackham/bin/mafft --globalpair --maxiterate 1000 $3 --thread 10 --unalignlevel 0.6 --reorder --adjustdirection $1 > $2

## mafft-linsi
##/sw/apps/bioinfo/MAFFT/7.407/rackham/bin/mafft-linsi --reorder --adjustdirection --thread 2 $1 > $2
