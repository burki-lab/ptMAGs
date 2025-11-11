################################################################################
#                           Amino acid usage across taxa                       #
################################################################################

## load required package
library("factoextra")
library("dplyr")
library("tidyr")
library("ggdendro")
library("dendextend")


## setwd
setwd("ptMAGs/01_phylogenies/results/concat/subset/trees/v14/")

# 93 gene dataset with 107 taxa
## load data
freq <- read.csv(file = "comp_het/concat_107_93g_aa_freq_group.txt", header = TRUE, sep = "\t")
## set first column as row names
freq <- data.frame(freq, row.names = 1)

## 1. Bar plot
### Convert to long format
freq_long <- freq %>% 
  pivot_longer(cols = `Freq.A.`:`Freq.H.`, 
    names_to = "amino_acid",
    values_to = "freq")

### Plot
ggplot(freq_long, aes(x = amino_acid, y = freq, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Amino acid frequencies per group", x = "Amino Acid", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

## Highest variability in amino acid I, A, N, and V.

## 2. PCA
### extract quantitative variables only
freq.subset <- freq[1:20]

### compute PCA
res.pca <- prcomp(freq.subset, scale = FALSE)

### generate scree plot
fviz_eig(res.pca, main="Scree plot - 93 genes")

### graph of variables
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             title = "93 genes")

## Points out amino acids I, K, V, A, and N

groups <- as.factor(freq$Group)

fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             addEllipses = TRUE, # Concentration ellipses
             palette = c("#d960a3",  "#ff7c70", "#51939a", "#ff6e00", "#e22f34", "#ffb800", "#ff7c70", "#7f9a66"),
             ellipse.type = "confidence",
             #ellipse.level = 0.68, 
             legend.title = "Groups",
             repel = FALSE,
             geom = "point",
             title = "93 genes"
)


# Create a PCA biplot with both variables and individuals
fviz_pca_biplot(res.pca,
                label = "var", # Display variable labels
                col.ind = groups, # Color individuals by groups
                addEllipses = TRUE, # Add confidence ellipses for groups
                palette = c("#d960a3", "#ff7c70", "#51939a", "#ff6e00", "#e22f34", "#ffb800", "#ff7c70", "#7f9a66"),
                ellipse.type = "confidence",
                legend.title = "Groups",
                repel = TRUE, # Avoid text overlapping
                title = "PCA - Biplot for 93 genes")


### The PCA results show cryptophytes and haptophytes clustering closer together. 
### Both cryptophytes and haptophytes have higher V frequencies, and lower N frequencies.


## 3. Recoding strategies

### Dayhoff 15
## Change columns so they reflect Dayhoff 15 recoding
## 0 = D, E, and Q
## 1 = M and L
## 2 = I and V
## 3 = F and Y
## The rest of the amino acids are treated individually so we do not need to touch those columns

freq.d15 <- freq.subset %>% mutate(Freq.0 = Freq.D. + Freq.E. + Freq.Q.)
freq.d15 <- freq.d15 %>% mutate(Freq.1 = Freq.M. + Freq.L.)
freq.d15 <- freq.d15 %>% mutate(Freq.2 = Freq.I. + Freq.V.)
freq.d15 <- freq.d15 %>% mutate(Freq.3 = Freq.F. + Freq.Y.)

## drop the columns we dont need
freq.d15 <- freq.d15 %>% select(-(Freq.D.:Freq.E.))
freq.d15 <- freq.d15 %>% select(-(Freq.Q.))
freq.d15 <- freq.d15 %>% select(-(Freq.M.))
freq.d15 <- freq.d15 %>% select(-(Freq.L.))
freq.d15 <- freq.d15 %>% select(-(Freq.I.))
freq.d15 <- freq.d15 %>% select(-(Freq.V.))
freq.d15 <- freq.d15 %>% select(-(Freq.F.))
freq.d15 <- freq.d15 %>% select(-(Freq.Y.))


## compute PCA
res.pca <- prcomp(freq.d15, scale = FALSE)

### generate scree plot
fviz_eig(res.pca, main="Scree plot - D15")

### graph of variables
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             title = "Dayhoff 15 recoding")

## Points out amino acids A, K, 0, 3, and N

groups <- as.factor(freq$Group)

fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             addEllipses = TRUE, # Concentration ellipses
             palette = c("#d960a3",  "#ff7c70", "#51939a", "#ff6e00", "#e22f34", "#ffb800", "#ff7c70", "#7f9a66"),
             ellipse.type = "confidence",
             #ellipse.level = 0.68,
             legend.title = "Groups",
             repel = TRUE,
             geom = "point",
             title = "Dayhoff 15"
)


# Create a PCA biplot with both variables and individuals
fviz_pca_biplot(res.pca,
                label = "var", # Display variable labels
                col.ind = groups, # Color individuals by groups
                addEllipses = FALSE, # Add confidence ellipses for groups
                palette = c("#d960a3", "#ff7c70", "#51939a", "#ff6e00", "#e22f34", "#ffb800", "#ff7c70", "#7f9a66"),
                #ellipse.type = "confidence",
                ellipse.level = 0.68,
                legend.title = "Groups",
                repel = TRUE, # Avoid text overlapping
                title = "PCA - Biplot for Dayhoff15 recoding")

## Dayhoff 15 recoding results in clusters that overlap with each other more, with the
## exception of Viridiplantae which cluster very far away from the rest

## Dayhoff 12 recoding

## Change columns so they reflect Dayhoff 12 recoding
## 0 = D, E, and Q, same as in D15
## 1 = M, L, I, and V, add two columns of D15
## 2 = F and Y, same as in D15 (Freq.3)
## 3 = K, H, and R

## Columns
## 0 = D, E, and Q
## 1 = M, L, I, and V
## 2 = K, H, and R
## 3 = F and Y
## The rest of the amino acids are treated individually so we do not need to touch those columns

freq.d12 <- freq.d15 %>% mutate(Freq.1 = Freq.1 + Freq.2)
freq.d12 <- freq.d12 %>% mutate(Freq.2 = Freq.K. + Freq.H. + Freq.R.)

## drop the columns we dont need
freq.d12 <- freq.d12 %>% select(-(Freq.K.))
freq.d12 <- freq.d12 %>% select(-(Freq.H.))
freq.d12 <- freq.d12 %>% select(-(Freq.R.))

## compute PCA
res.pca <- prcomp(freq.d12, scale = FALSE)

### generate scree plot
fviz_eig(res.pca, main="Scree plot - D12")

### graph of variables
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             title = "Dayhoff 12 recoding")

## Points out amino acids A, 0, 3, and N

groups <- as.factor(freq$Group)

fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             addEllipses = TRUE, # Concentration ellipses
             palette = c("#d960a3",  "#ff7c70", "#51939a", "#ff6e00", "#e22f34", "#ffb800", "#ff7c70", "#7f9a66"),
             ellipse.type = "confidence",
             #ellipse.level = 0.68,
             legend.title = "Groups",
             repel = TRUE,
             geom = "point",
             title = "Dayhoff 12"
)


# Create a PCA biplot with both variables and individuals
fviz_pca_biplot(res.pca,
                label = "var", # Display variable labels
                col.ind = groups, # Color individuals by groups
                addEllipses = TRUE, # Add confidence ellipses for groups
                palette = c("#d960a3", "#ff7c70", "#51939a", "#ff6e00", "#e22f34", "#ffb800", "#ff7c70", "#7f9a66"),
                ellipse.type = "confidence",
                legend.title = "Groups",
                repel = TRUE, # Avoid text overlapping
                title = "PCA - Biplot for Dayhoff12 recoding")

## Dayhoff 12 somehow looks worse than Dayhoff 15. Viridiplantae continues to be a an oddball.


## Dayhoff 9 recoding

## Change columns so they reflect Dayhoff 9 recoding
## 0 = D, E, Q, H, and N
## 1 = I, L, M, and V, add two columns of D15
## 2 = F and Y, same as in D15 (Freq.3)
## 3 = A, S, and T
## 4 = K and R

## Columns
## 0 = D, E, Q, H, and N
## 1 = M, L, I, and V
## 2 = A, S, and T
## 3 = F and Y
## 4 = K and R
## The rest of the amino acids are treated individually so we do not need to touch those columns

freq.d9 <- freq.d15 %>% mutate(Freq.0 = Freq.0 + Freq.H. + Freq.N.)
freq.d9 <- freq.d9 %>% mutate(Freq.1 = Freq.1 + Freq.2)
freq.d9 <- freq.d9 %>% mutate(Freq.2 = Freq.A. + Freq.S. + Freq.T.)
freq.d9 <- freq.d9 %>% mutate(Freq.4 = Freq.K. + Freq.R.)

## drop the columns we dont need
freq.d9 <- freq.d9 %>% select(-(Freq.H.))
freq.d9 <- freq.d9 %>% select(-(Freq.N.))
freq.d9 <- freq.d9 %>% select(-(Freq.A.))
freq.d9 <- freq.d9 %>% select(-(Freq.S.))
freq.d9 <- freq.d9 %>% select(-(Freq.T.))
freq.d9 <- freq.d9 %>% select(-(Freq.K.))
freq.d9 <- freq.d9 %>% select(-(Freq.R.))

## compute PCA
res.pca <- prcomp(freq.d9, scale = FALSE)

### generate scree plot
fviz_eig(res.pca, main="Scree plot - D9")

### graph of variables
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             title = "Dayhoff 9 recoding")

## Points out amino acids 2, and 3

groups <- as.factor(freq$Group)

fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             addEllipses = TRUE, # Concentration ellipses
             palette = c("#d960a3",  "#ff7c70", "#51939a", "#ff6e00", "#e22f34", "#ffb800", "#ff7c70", "#7f9a66"),
             ellipse.type = "confidence",
             #ellipse.level = 0.68,
             legend.title = "Groups",
             repel = TRUE,
             geom = "point",
             title = "Dayhoff 9"
)


# Create a PCA biplot with both variables and individuals
fviz_pca_biplot(res.pca,
                label = "var", # Display variable labels
                col.ind = groups, # Color individuals by groups
                addEllipses = TRUE, # Add confidence ellipses for groups
                palette = c("#d960a3", "#ff7c70", "#51939a", "#ff6e00", "#e22f34", "#ffb800", "#ff7c70", "#7f9a66"),
                ellipse.type = "confidence",
                legend.title = "Groups",
                repel = TRUE, # Avoid text overlapping
                title = "PCA - Biplot for Dayhoff9 recoding")

## Dayhoff 6 recoding

## Change columns so they reflect Dayhoff 6 recoding
## 0 = A, G, P, S, T
## 1 = D, E, N, Q
## 2 = H, K, R
## 3 = I, L, M, V
## 4 = F, W, Y
## 5 = C


freq.d6 <- freq.subset %>% mutate(Freq.0 = Freq.A. + Freq.G. + Freq.P. + Freq.S. + Freq.T.)
freq.d6 <- freq.d6 %>% mutate(Freq.1 = Freq.D. + Freq.E. + Freq.N. + Freq.Q.)
freq.d6 <- freq.d6 %>% mutate(Freq.2 = Freq.H. + Freq.K. + Freq.R.)
freq.d6 <- freq.d6 %>% mutate(Freq.3 = Freq.I. + Freq.L. + Freq.M. + Freq.V.)
freq.d6 <- freq.d6 %>% mutate(Freq.4 = Freq.F. + Freq.W. + Freq.Y.)

## drop the columns we dont need
freq.d6 <- freq.d6 %>% select(-(Freq.A.:Freq.T.))
freq.d6 <- freq.d6 %>% select(-(Freq.N.:Freq.H.))


## compute PCA
res.pca <- prcomp(freq.d6, scale = FALSE)

### generate scree plot
fviz_eig(res.pca, main="Scree plot - D6")

### graph of variables
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             title = "Dayhoff 6 recoding")

## Points out amino acids 0

groups <- as.factor(freq$Group)

fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             addEllipses = TRUE, # Concentration ellipses
             palette = c("#d960a3",  "#ff7c70", "#51939a", "#ff6e00", "#e22f34", "#ffb800", "#ff7c70", "#7f9a66"),
             ellipse.type = "confidence",
             #ellipse.level = 0.68,
             legend.title = "Groups",
             repel = TRUE,
             geom = "point",
             title = "Dayhoff 6"
)


# Create a PCA biplot with both variables and individuals
fviz_pca_biplot(res.pca,
                label = "var", # Display variable labels
                col.ind = groups, # Color individuals by groups
                addEllipses = TRUE, # Add confidence ellipses for groups
                palette = c("#d960a3", "#ff7c70", "#51939a", "#ff6e00", "#e22f34", "#ffb800", "#ff7c70", "#7f9a66"),
                #ellipse.type = "confidence",
                ellipse.level = 0.68,
                legend.title = "Groups",
                repel = TRUE, # Avoid text overlapping
                title = "PCA - Biplot for Dayhoff6 recoding")


## Not so great

## SR6 recoding
## 0 = APST
## 1 = DENG
## 2 = QKR
## 3 = MIVL
## 4 = WC
## 5 = FYH

freq.sr6 <- freq.subset %>% mutate(Freq.0 = Freq.A. + Freq.P. + Freq.S. + Freq.T.)
freq.sr6 <- freq.sr6 %>% mutate(Freq.1 = Freq.D. + Freq.E. + Freq.N. + Freq.G.)
freq.sr6 <- freq.sr6 %>% mutate(Freq.2 = Freq.Q. + Freq.K. + Freq.R.)
freq.sr6 <- freq.sr6 %>% mutate(Freq.3 = Freq.I. + Freq.L. + Freq.M. + Freq.V.)
freq.sr6 <- freq.sr6 %>% mutate(Freq.4 = Freq.C. + Freq.W.)
freq.sr6 <- freq.sr6 %>% mutate(Freq.5 = Freq.F. + Freq.Y. + Freq.H.)

## drop the columns we dont need
freq.sr6 <- freq.sr6 %>% select(-(Freq.A.:Freq.H.))

## compute PCA
res.pca <- prcomp(freq.sr6, scale = FALSE)

### generate scree plot
fviz_eig(res.pca, main="Scree plot - SR6")

### graph of variables
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             title = "SR 6 recoding")

## Points out amino acids 0, 1, 5

groups <- as.factor(freq$Group)

fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             addEllipses = TRUE, # Concentration ellipses
             palette = c("#d960a3",  "#ff7c70", "#51939a", "#ff6e00", "#e22f34", "#ffb800", "#ff7c70", "#7f9a66"),
             ellipse.type = "confidence",
             #ellipse.level = 0.68,
             legend.title = "Groups",
             repel = TRUE,
             geom = "point",
             title = "SR 6"
)


# Create a PCA biplot with both variables and individuals
fviz_pca_biplot(res.pca,
                label = "var", # Display variable labels
                col.ind = groups, # Color individuals by groups
                addEllipses = TRUE, # Add confidence ellipses for groups
                palette = c("#d960a3", "#ff7c70", "#51939a", "#ff6e00", "#e22f34", "#ffb800", "#ff7c70", "#7f9a66"),
                ellipse.type = "confidence",
                legend.title = "Groups",
                repel = TRUE, # Avoid text overlapping
                title = "PCA - Biplot for SR6 recoding")


## Not so great


## SR4 recoding
## A = AGNPST
## C = CHWY
## G = DEKQR
## T = FILMV

freq.sr4 <- freq.subset %>% mutate(Freq.1. = Freq.A. + Freq.G. + Freq.N. + Freq.P. + Freq.S. + Freq.T.)
freq.sr4 <- freq.sr4 %>% mutate(Freq.2. = Freq.C. + Freq.H. + Freq.W. + Freq.Y.)
freq.sr4 <- freq.sr4 %>% mutate(Freq.3. = Freq.D. + Freq.E. + Freq.K. + Freq.Q. + Freq.R.)
freq.sr4 <- freq.sr4 %>% mutate(Freq.4. = Freq.F. + Freq.I. + Freq.L. + Freq.M. + Freq.V.)

## drop the columns we dont need
freq.sr4 <- freq.sr4 %>% select(-(Freq.A.:Freq.H.))

## compute PCA
res.pca <- prcomp(freq.sr4, scale = FALSE)

### generate scree plot
fviz_eig(res.pca, main="Scree plot - SR4")

### graph of variables
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             title = "SR 4 recoding")

## Points out amino acids 3

groups <- as.factor(freq$Group)

fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             addEllipses = TRUE, # Concentration ellipses
             palette = c("#d960a3",  "#ff7c70", "#51939a", "#ff6e00", "#e22f34", "#ffb800", "#ff7c70", "#7f9a66"),
             ellipse.type = "confidence",
             #ellipse.level = 0.68,
             legend.title = "Groups",
             repel = TRUE,
             geom = "point",
             title = "SR 4"
)


# Create a PCA biplot with both variables and individuals
fviz_pca_biplot(res.pca,
                label = "var", # Display variable labels
                col.ind = groups, # Color individuals by groups
                addEllipses = TRUE, # Add confidence ellipses for groups
                palette = c("#d960a3", "#ff7c70", "#51939a", "#ff6e00", "#e22f34", "#ffb800", "#ff7c70", "#7f9a66"),
                ellipse.type = "confidence",
                legend.title = "Groups",
                repel = TRUE, # Avoid text overlapping
                title = "PCA - Biplot for SR4 recoding")


## Not so great



## 3. UPGMA clustering
### Calculate the distance matrix using Euclidean distance
dist_matrix <- dist(freq.subset, method = "euclidean")

### Perform hierarchical clustering
hclust_res <- hclust(dist_matrix, method = "average")

### Convert the hclust object into a dendrogram for ggplot2
dendro_data <- dendro_data(as.dendrogram(hclust_res), type = "rectangle")

dendro_ggplot <- ggdendrogram(dendro_data, rotate = FALSE, theme_dendro = TRUE)

### Plot the dendrogram
dendro_ggplot + 
  theme(axis.text.x = element_text(size = 5))

### Get group info
freq <- read.csv(file = "comp_het/concat_107_93g_aa_freq_group.txt", header = TRUE, sep = "\t")

groups <- as.factor(freq$Group)

group_data <- data.frame(label = freq$NAME,
                         group = groups)


### Merge group data with the labels data from `ggdendro`
dendro_data$labels <- dendro_data$labels %>%
  left_join(group_data, by = "label")


# Create the dendrogram plot
dendro_ggplot <- ggplot() + 
  geom_segment(data = dendro_data$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = dendro_data$labels, 
             aes(x = x, y = y, color = group), 
             size = 3) + 
  theme_minimal() +
  scale_color_manual(values=c("#d960a3", "#ff7c70", "#51939a", "#ff6e00", "#e22f34", "#ffb800", "#ff7c70", "#7f9a66")) +
  theme(axis.text.y = element_text(size = 10)) +
  labs(color = "Group", title = "Hierarchical Clustering with Group Annotations")

# Print the dendrogram with group annotations
print(dendro_ggplot)
