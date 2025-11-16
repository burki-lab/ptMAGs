################################################################################
#             Heatmap of gene presence absence in plastid genomes              #
################################################################################

## load required package
library("phytools")
library("tidyverse")
library("pheatmap")
library("dplyr")

# set working directory
setwd("ptMAGs/03_comp_genomics/results/gene_content/")

# 1. Cleaning the data frame
## Load the manually curated HOGs file
gene_presence_absence <- readxl::read_excel("orthofinder/N0_orthofinder_HOGs_manually_checked.xlsx")

## Remove HOGs with less than 3 taxa
gene_presence_absence <- gene_presence_absence[gene_presence_absence$No_of_taxa > 2, ]
### We are left with 237 HOGs with sequences from at least 3 taxa. 

## Remove columns that are not needed
gene_presence_absence <- gene_presence_absence %>% select(-c(HOG, OG, `Gene Tree Parent Clade`, No_of_taxa, Synonym, Note))

gene_presence_absence <- gene_presence_absence %>% select(c(Gene, `Lepto_01_TARA_B110000977_METAG-scaffold_147`))

## Replace NAs with 0, and other text with 1 to indicate presence/absence
### Set gene names as row names so that they do not get replaced by 1s.
gene_presence_absence <- data.frame(gene_presence_absence, row.names = 1)
gene_presence_absence[is.na(gene_presence_absence)] <- 0
gene_presence_absence[gene_presence_absence != '0'] = '1'

## Transpose dataframe so that row names are taxa, and column names are gene names
df <- t(gene_presence_absence)
df <- as.data.frame(df)

## Set rownames as first column now
df <- rownames_to_column(df, "Taxa")

## Replace row names so that they are shorter
### Read in correspondence file
correspondence <- readxl::read_excel("orthofinder/replace_taxa_names.xlsx")

### Merge the two dataframes
df <- inner_join(df, correspondence)

## Subset to remove original taxon names
df <- df[-1]
  
## Set column with tip labels as row names
df <- data.frame(df, row.names = 239)


# 2. Heatmap with UPGMA clustering

## Subset to remove group
df_subset <- df[-238]

## Convert all 0s and 1s to numeric
df_subset <- df_subset %>%
  mutate(across(everything(), as.numeric))

## Plot heatmap
pheatmap(df_subset,
         cluster_rows = TRUE,   
         cluster_cols = TRUE,   
         scale = "none",        # No scaling as it's binary data (0s and 1s)
         color = c("white", "blue"),  
         show_rownames = TRUE, 
         show_colnames = TRUE,
         clustering_method = "average",
         legend = FALSE,
         fontsize = 5) 

## Annotate taxonomic groups
### Keep only the new taxonomic names and taxonomic group
correspondence <- correspondence[-1]

### Set the names as row names
correspondence <- data.frame(correspondence, row.names = 2)

### Set colours
cols<- list(Group = c("Cryptophyta" = "#d960a3", 
        "Haptophyta" = "#ff6e00", 
        "Leptophyte" = "#e22f34", 
        "Ochrophyta" = "#ffb800",
        "Rhodophyta" = "#ff7c70"))


## without gene names
pheatmap(df_subset,
         cluster_rows = TRUE,   
         cluster_cols = TRUE,   
         scale = "none",        # No scaling as it's binary data (0s and 1s)
         color = c("#ece7f2", "#2b8cbe"),  
         show_rownames = TRUE, 
         show_colnames = FALSE,
         annotation_row = correspondence,
         clustering_method = "average",
         legend = FALSE,
         fontsize_row = 8,
         annotation_colors = cols,
         annotation_legend = FALSE)


## with gene names
pheatmap(df_subset,
         cluster_rows = TRUE,   
         cluster_cols = TRUE,   
         scale = "none",        # No scaling as it's binary data (0s and 1s)
         color = c("#ece7f2", "#2b8cbe"),  
         show_rownames = TRUE, 
         show_colnames = TRUE,
         annotation_row = correspondence,
         clustering_method = "average",
         legend = FALSE,
         fontsize_row = 8,
         fontsize_col = 5,
         annotation_colors = cols,
         annotation_legend = FALSE)

