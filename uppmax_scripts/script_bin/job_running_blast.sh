#!/bin/bash
#SBATCH -A uppmax2025-2-76
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 10:00:00
#SBATCH --mail-user mahwash.jamy@slu.se
#SBATCH --mail-type=FAIL
#SBATCH -J job_2023_10_17_blast
#SBATCH -o /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.out
#SBATCH -e /crex/proj/naiss2023-6-81/Mahwash/ptMAGs/uppmax_scripts/out_bin/slurm-%A.err

# modules
module load bioinfo-tools
module load blast

## takes in two variables
## $1 = the fasta file (with the absolute path)
## $2 = the blastdb_name.txt file containing the names of the blastdbs (with absolute paths)

# commands
cat $2 | while read line; do
    taxon=$(basename $line | cut -f 1 -d '.')
    blastp -query $1 -num_threads 6 -evalue 1e-02 -db $line -out "$1"__"$taxon".blastout -outfmt '6 std salltitles'
done

#cat $2 | while read line; do
#    taxon=$(basename $line | cut -f 1 -d '.')
#    blastp -query $1 -num_threads 6 -evalue 1e-01 -db $line -out "$1"__"$taxon".blastout -outfmt '6 std salltitles'
#done

