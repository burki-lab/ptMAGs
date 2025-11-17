#!/bin/bash
#SBATCH -A naiss2024-5-197
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 10:00:00
#SBATCH --mail-user mahwash.jamy@slu.se
#SBATCH --mail-type=FAIL
#SBATCH -J job_2023_10_22_blastn
#SBATCH -o /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.out
#SBATCH -e /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.err

# modules
module load bioinfo-tools
module load blast

## takes in three variables
## $1 = the fasta file (with the absolute path)
## $2 = the blastdb file (with the absolute path) 
## $3 = output directory

# commands
taxon=$(basename $1 | cut -f 1 -d '.')
db=$(basename $2 | cut -f 1 -d '.')

blastn -query $1 -num_threads 6 -evalue 1e-10 -db $2 -out "$3"/"$taxon"__"$db".blastout -outfmt '6 std salltitles'

