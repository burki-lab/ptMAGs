################################################################################################
#                       Rscript Comparison between MAGs and references                         #
################################################################################################

## load required package
library("ggplot2")
library("dplyr")
library("hrbrthemes")
library("wesanderson")
library("ggpubr")

## set working directory
setwd("ptMAGs/00_plastid_MAGS/results/stats/")

# 1. Completeness level

## We use the presence/absence of 44 core to estimate completeness of the MAGs and the references. 

### Load file 
completeness <- read.csv("completeness_taxa.txt",
                         sep = "\t",
                         header = TRUE)

completeness$Type <- factor(completeness$Type , levels=c("ref", "mag"))

### Make boxplot
p1 <- ggplot(completeness, aes(x=Type, y=Completeness.percentage, fill=Type)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = wes_palette("Royal1")) + 
  geom_jitter(color="black", size=0.6, alpha=0.9) +
  theme_minimal() +
  xlab("") + ylab("Completeness (%)") +
  scale_x_discrete(labels=c("ref" = "References", "mag" = "MAGs")) +
  theme(legend.position = "none") +
  ggtitle("Plastid genome completeness based on 44 core genes")

p1


## Let us looks at completeness as a function of taxonomy

### Load file
taxonomy <- read.csv("taxonomy_taxa.txt",
                         sep = "\t",
                         header = TRUE)

### Merge dataframes
df_merge <- merge(completeness, taxonomy, by="Taxon") 

### Exclude Paulinella
df_merge <- filter(df_merge, Group!="Paulinella")


### Plot barplot
p2 <- ggplot(df_merge, aes(x=Group, y=Completeness.percentage, fill=Type)) + 
  geom_boxplot() +
  scale_fill_manual(values = wes_palette("Royal1"), labels = c("Reference", "MAG")) + 
  theme_minimal() +
  ylab("Completeness (%)") +
  xlab("") +
  facet_wrap(~Group, scale="free") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank()) +
  ggtitle("Completeness by group")

p2


## Merge the completeness plot
p <- ggarrange(p1, p2, ncol = 1, labels = "auto")
p

ggsave(p, filename = "completeness.pdf", 
       device = cairo_pdf)

# 2. Size

## How big are the MAGs and references on average?

### Load file 
size <- read.csv("size_taxa.txt",
                         sep = "\t",
                         header = TRUE)

size$Type <- factor(completeness$Type , levels=c("ref", "mag"))

size <- size %>% 
  mutate(Size_kbp = Size/1e3)

### Make boxplot
#### Force non-scientific notation
options(scipen = 999)

p1 <- ggplot(size, aes(x=Type, y=Size_kbp, fill=Type)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = wes_palette("Royal1")) + 
  geom_jitter(color="black", size=0.7, alpha=0.8) +
  theme_minimal() +
  xlab("") + ylab("Size (kbp)") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_discrete(labels=c("ref" = "References", "mag" = "MAGs")) +
  theme(legend.position = "none")

p1

### Redo boxplot without Paulinella!
p1 <- size %>% 
  filter(Taxon!="NC_039737") %>%
  ggplot(aes(x=Type, y=Size_kbp, fill=Type)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = wes_palette("Royal1")) + 
  geom_jitter(color="black", size=0.7, alpha=0.8) +
  theme_minimal() +
  xlab("") + ylab("Size (kbp)") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_discrete(labels=c("ref" = "References", "mag" = "MAGs")) +
  theme(legend.position = "none") +
  ggtitle("Plastid genome size")

p1

## Let's look at size by taxonomic group

### Merge dataframes
df_merge <- merge(size, taxonomy, by="Taxon") 

### Exclude Paulinella
df_merge <- filter(df_merge, Group!="Paulinella")

p2 <- ggplot(df_merge, aes(x=Group, y=Size_kbp, fill=Type)) + 
  geom_boxplot() +
  scale_fill_manual(values = wes_palette("Royal1"), labels = c("Reference", "MAG")) + 
  theme_minimal() +
  xlab("") + ylab("Size (kbp)") +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~Group, scale="free") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank()) +
  ggtitle("Genome size across groups")

p2

## Merge the plots
p <- ggarrange(p1, p2, ncol = 1, labels = "auto")
p

ggsave(p, filename = "size.pdf", 
       path = "ptMAGs/00_plastid_MAGS/results/stats/",
       device = cairo_pdf)



# 3. GC content

## What is the GC content of the MAGs and references on average?

### Load file 
gc <- read.csv("gc_taxa.txt",
                 sep = "\t",
                 header = TRUE)

gc$Type <- factor(completeness$Type , levels=c("ref", "mag"))

### Make boxplot
p <- ggplot(gc, aes(x=Type, y=GCcontent, fill=Type)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = wes_palette("Royal1")) + 
  geom_jitter(color="black", size=0.7, alpha=0.8) +
  theme_minimal() +
  xlab("") + ylab("GC content (%)") +
  scale_x_discrete(labels=c("ref" = "References", "mag" = "MAGs")) +
  theme(legend.position = "none")

p

## Probably a lot more interesting to do this by taxonomic group!

### Merge dataframes
df_merge <- merge(gc, taxonomy, by="Taxon") 

### Exclude Paulinella
df_merge <- filter(df_merge, Group!="Paulinella")

p <- ggplot(df_merge, aes(x=Group, y=GCcontent, fill=Type)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = wes_palette("Royal1"), labels = c("Reference", "MAG")) +
  geom_jitter(color="black", size=0.7, alpha=0.6) +
  theme_minimal() +
  ylab("GC Content (%)") +
  xlab("") +
  facet_wrap(~Group, scale="free") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank()) +
  ggtitle("GC content of plastid genomes")

p

ggsave(p, filename = "gc.pdf", 
       path = "ptMAGs/00_plastid_MAGS/results/stats/",
       device = cairo_pdf)



### A bit of a difference but not too much in my opinion. 
### Let's just look at taxonomic groups without considering type. 

p <- ggplot(df_merge, aes(x=reorder(Group, GCcontent), y=GCcontent, fill=Group)) + 
  geom_boxplot() +
  theme_minimal() +
  xlab("") + ylab("GC content (%)") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#8ba976",
                             "#8ba976",
                             "#ffb900",
                             "#ee72b2",
                             "#ffb900",
                             "#8ba976",
                             "#54a2a9",
                             "#ff8e49",
                             "#f93647",
                             "#ffb900",
                             "#006376",
                             "#ff7f78", 
                             "#434e75"))
                             
                             
                        
p


# 4. No. of genes

## How many genes are encoded on the references and MAGs on average?

### Load file
genes <- read.csv("no-of-genes_taxa.txt",
               sep = "\t",
               header = TRUE)

genes$Type <- factor(genes$Type , levels=c("ref", "mag"))

### Make boxplot
p <- ggplot(genes, aes(x=Type, y=Number.of.genes, fill=Type)) + 
  geom_boxplot() +
  scale_fill_manual(values = wes_palette("Royal1")) + 
  theme_minimal() +
  xlab("") + ylab("Number of genes") +
  scale_x_discrete(labels=c("ref" = "References", "mag" = "MAGs")) +
  theme(legend.position = "none")

p

### Redo boxplot without Paulinella!
p1 <- genes %>% 
  filter(Taxon!="NC_039737") %>%
  ggplot(aes(x=Type, y=Number.of.genes, fill=Type)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = wes_palette("Royal1")) + 
  geom_jitter(color="black", size=0.7, alpha=0.8) +
  theme_minimal() +
  xlab("") + ylab("Number of genes") +
  scale_x_discrete(labels=c("ref" = "References", "mag" = "MAGs")) +
  theme(legend.position = "none") +
  ggtitle("Number of genes in plastid genomes")

p1


## Probably a lot more interesting to do this by taxonomic group!

### Merge dataframes
df_merge <- merge(genes, taxonomy, by="Taxon") 

### Exclude Paulinella
df_merge <- filter(df_merge, Group!="Paulinella")

p2 <- ggplot(df_merge, aes(x=Group, y=Number.of.genes, fill=Type)) + 
  geom_boxplot() +
  scale_fill_manual(values = wes_palette("Royal1"), labels = c("Reference", "MAG")) +
  theme_minimal() +
  ylab("Number of genes") +
  xlab("") +
  facet_wrap(~Group, scale="free") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank()) +
  ggtitle("Number of genes in different plastid groups")

p2


## Merge the gene numbers plot
p <- ggarrange(p1, p2, ncol = 1, labels = "auto")
p

ggsave(p, filename = "number-of-genes.pdf", 
       path = "ptMAGs/00_plastid_MAGS/results/stats/",
       device = cairo_pdf)

# 5. Abundance

## How abundant are the references and MAGs on average in the global oceans?

### Load file 
abundance <- read.csv("abundance_taxa.txt", sep = "\t", header = TRUE)


abundance$MAG <- factor(abundance$MAG , levels=c("ref", "mag"))

### Make boxplot
p <- ggplot(abundance, aes(x=MAG, y=Total_abundance, fill=MAG)) + 
  geom_boxplot() +
  scale_fill_manual(values = wes_palette("Royal1")) + 
  theme_minimal() +
  xlab("") + ylab("Abundance (coverage)") +
  scale_x_discrete(labels=c("ref" = "References", "mag" = "MAGs")) +
  theme(legend.position = "none")

p

### Let's remove non-marine alga 
abundance <- abundance[abundance$Total_abundance > 0, ]

## That leaves 720 taxa

### Regular boxplot
p <- ggplot(abundance, aes(x=MAG, y=Total_abundance, fill=MAG)) + 
  geom_boxplot() +
  scale_fill_manual(values = wes_palette("Royal1")) + 
  theme_minimal() +
  xlab("") + ylab("Abundance (coverage)") +
  scale_x_discrete(labels=c("ref" = "References", "mag" = "MAGs")) +
  theme(legend.position = "none")

p

### Log scale
p1 <- ggplot(abundance, aes(x=MAG, y=Total_abundance, fill=MAG)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  scale_fill_manual(values = wes_palette("Royal1")) + 
  theme_minimal() +
  scale_y_log10() +
  xlab("") + ylab("Abundance (mean coverage, log scale)") +
  scale_x_discrete(labels=c("ref" = "References", "mag" = "MAGs")) +
  theme(legend.position = "none") +
  ggtitle("Abundance of plastid genomes in the sunlit ocean")

p1


### Let's look at it by group
p2 <- ggplot(abundance, aes(x=Clade, y=Total_abundance, fill=MAG)) + 
  geom_boxplot() +
  scale_fill_manual(values = wes_palette("Royal1"), labels = c("Reference", "MAG")) + 
  theme_minimal() +
  ylab("Abundance (mean coverage)") +
  xlab("") +
  facet_wrap(~Clade, scale="free") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank(),
        legend.position = "none") +
  ggtitle("Abundance of different plastid groups")

p2

## Merge the coverage plot
p <- ggarrange(p1, p2, ncol = 1, labels = "auto")
p

ggsave(p, filename = "coverage.pdf", 
       path = "/ptMAGs/00_plastid_MAGS/results/stats/",
       device = cairo_pdf)
