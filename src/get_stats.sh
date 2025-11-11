#!/bin/bash

## A quick and dirty script to get gene, taxa, and occupancy stats for a phylogenomic matrix
## Usage `get_stats.sh {path to directory containing fasta files for every gene}`.

## Go to user provided directory
cd $1

## Check if the "stats" directory exists. If not, create it.
if [ ! -d "stats" ]; then
    mkdir "stats"
fi

## Get taxa stats
printf "ID\tTaxonomy\tGene_count\n" > stats/taxa_stats.tsv

grep "mag=" *fasta | \
sed -E 's/.*>(.*)/\1/' | \
sort | \
uniq -c | \
sed -E 's/\s+([0-9]+) (mag=.*taxonomy=)(.*)/\2\3\t\3\t\1/' \
>> stats/taxa_stats.tsv

grep "taxo=" *fasta | \
sed -E 's/.*>(.*)/\1/' | \
sort | \
uniq -c | \
sed -E 's/\s+([0-9]+) (taxo=.*)/\2\tref\t\1/' \
>> stats/taxa_stats.tsv

## Get gene stats
printf "Gene\tNo_of_taxa\tPerc_of_taxa\n" > stats/gene_stats.tsv

Total=$(grep ">" *fasta | sed -E 's/.*>(.*)/\1/' | sort -u | wc -l)

for i in *fasta
do 
    gene=$(basename $i | cut -f1 -d '.')
    num=$(grep -c ">" $i)
    perc=$(echo $num $Total | awk '{percentage = ($1 / $2) * 100; printf "%.2f%\n", percentage}')
    printf "$gene\t$num\t$perc%\n" >> stats/gene_stats.tsv
done


## Get occupancy table
gene_list=($(ls *.fasta | sed -E 's/\.fasta//'))

# Generate the header line with gene names
header="Taxa\t${gene_list[@]}"

# Initialize the content with taxa and zeros
content=""

# Collect all unique taxa across all fasta files
all_taxa=$(awk '/^>/{print substr($1, 2)}' "${gene_list[@]/%/.fasta}" | sort -u)

# Iterate through taxa and genes, check presence, and append 0 or 1
for taxa in ${all_taxa}; do
    line="${taxa}"
    for gene in "${gene_list[@]}"; do
        presence=$(grep -q "${taxa}" "${gene}.fasta" && echo 1 || echo 0)
        line+="\t${presence}"
    done
    content+="${line}\n"
done

# Combine header and content and write to occupancy.tsv
echo -e "$header\n$content" > stats/occupancy.tsv