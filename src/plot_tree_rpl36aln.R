################################################################################
#                         Plot alignment next to phylogeny                     #
################################################################################

## Here we want to plot the rpl36 alignment next to the subset phylogeny based 
## on 93 genes

# Load packages
library("tidyverse")
library("ggmsa")
library("ggtree")


# Set working directory
setwd("ptMAGs/01_phylogenies/results/")

## read tree
tree <- read.tree("concat/subset/trees/v14/catgtr/concat_107t_93g_prequal_ginsi_bmge_chain2-3_2000-9000_0.047.rooted.tre")


## Plot with alignment
msaplot(p=ggtree(tree), fasta="rpl36/alignments/trimmed/v3/rpl36.mafft.trimmed.fasta", width = 0.5) +
  geom_treescale()

fasta <- "rpl36/alignments/trimmed/v3/rpl36.mafft.trimmed.fasta"


ggmsa(fasta, color = "Chemistry_AA", font = NULL)
ggmsa(fasta, color = "Shapely_AA", font = NULL)
ggmsa(fasta, color = "Taylor_AA", font = NULL)
ggmsa(fasta, color = "Zappo_AA", font = NULL)
ggmsa(fasta, color = "LETTER", font = NULL)
