################################################################################
#                     Consensus trees of alternate topologies                  #
################################################################################


# load required package
library("ape")
library("phytools")


# set working directory
setwd("ptMAGs/01_phylogenies/results/concat/subset/trees/v14/meow/")

# read in trees
new_sister_h_c_complex_mono <- read.newick(file = "concat_107t_93g_LG-MEOW6020-G_new-sister-h-c_complex-mono.treefile")
new_sister_h_c_complex_nonmono <- read.newick(file = "concat_107t_93g_LG-MEOW6020-G_new-sister-h-c_complex-non-mono.treefile")
new_sister_c_complex_mono <- read.newick(file = "concat_107t_93g_LG-MEOW6020-G_new-sister-c_complex-mono.treefile")
new_sister_c_complex_nonmono <- read.newick(file = "concat_107t_93g_LG-MEOW6020-G_new-sister-c_complex-non-mono_2.treefile")
new_sister_h_complex_mono <- read.newick(file = "concat_107t_93g_LG-MEOW6020-G_new-sister-h_complex-mono_2.treefile")
new_sister_h_complex_nonmono <- read.newick(file = "concat_107t_93g_LG-MEOW6020-G_new-sister-h_complex-non-mono.treefile")

# root the trees
## outgroup
glauco <- c("taxo_Eukaryota-Glaucocystophyceae-Cyanophora-paradoxa-accession_NC-001675",
            "taxo_Eukaryota-Glaucocystophyceae-Cyanophora-biloba-accession_NC-038216",
            "taxo_Eukaryota-Glaucocystophyceae-Glaucocystis-sp-BBH-accession_MF167424",
            "taxo_Eukaryota-Glaucocystophyceae-Cyanoptyche-gloeocystis-accession_MF167427")

new_sister_h_c_complex_mono <- root(new_sister_h_c_complex_mono, outgroup = glauco, resolve.root = TRUE)
new_sister_h_c_complex_nonmono <- root(new_sister_h_c_complex_nonmono, outgroup = glauco, resolve.root = TRUE)
new_sister_c_complex_mono <- root(new_sister_c_complex_mono, outgroup = glauco, resolve.root = TRUE)
new_sister_c_complex_nonmono <- root(new_sister_c_complex_nonmono, outgroup = glauco, resolve.root = TRUE)
new_sister_h_complex_mono <- root(new_sister_h_complex_mono, outgroup = glauco, resolve.root = TRUE)
new_sister_h_complex_nonmono <- root(new_sister_h_complex_nonmono, outgroup = glauco, resolve.root = TRUE)

# Likelihood ratio test under GTR+CAT-PMSF model
## Here the best likelihood score is of the new_sister_h_c_complex_mono tree, so we treat it as the alternative tpology in all cases
## We compare this topology with all the other topologies

## new_sister_h_c_complex_mono vs. new_sister_h_c_complex_nonmono
## Null hypothesis topology = new_sister_h_c_complex_nonmono
## Alternate hypothesis topology = new_sister_h_c_complex_mono

con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_h_c_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 3
A_t0 <- 3^d

p_value <- 1 - pchisq(12.212, df = d)
p_value ## 0.006691122

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 0.1657878

## new_sister_h_c_complex_mono vs. new_sister_h_complex_mono
con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_h_complex_mono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 2
A_t0 <- 3^d

p_value <- 1 - pchisq(21.676, df = d)
p_value ## 1.963887e-05

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 0.0001767359

## new_sister_h_c_complex_mono vs. new_sister_h_complex_nonmono
con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_h_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 4
A_t0 <- 3^d

p_value <- 1 - pchisq(40.74, df = d)
p_value ## 3.042467e-08

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value 

## new_sister_h_c_complex_mono vs. new_sister_c_complex_mono
con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_c_complex_mono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 1
A_t0 <- 3^d

p_value <- 1 - pchisq(31.18, df = d)
p_value ## 2.351766e-08

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 7.033525e-08

## new_sister_h_c_complex_mono vs. new_sister_c_complex_nonmono
con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_c_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 4
A_t0 <- 3^d

p_value <- 1 - pchisq(51.73, df = d)
p_value ## 1.570941e-10

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 1.272462e-08



# Likelihood ratio test under LG+MEOW(60,20)+G model
## Here the best likelihood score is of the new_sister_h_c_complex_mono tree, so we treat it as the alternative topology in all cases
## We compare this topology with all the other topologies

con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_h_c_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 3
A_t0 <- 3^d

p_value <- 1 - pchisq(2.622, df = d)
p_value ## 0.4536457

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 0.9999999

## new_sister_h_c_complex_mono vs. new_sister_h_complex_mono
con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_h_complex_mono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 2
A_t0 <- 3^d

p_value <- 1 - pchisq(7.358, df = d)
p_value ## 0.02524821

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 0.2055869

## new_sister_h_c_complex_mono vs. new_sister_h_complex_nonmono
con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_h_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 4
A_t0 <- 3^d

p_value <- 1 - pchisq(21.46, df = d)
p_value ## 0.0002566364

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value 

## new_sister_h_c_complex_mono vs. new_sister_c_complex_mono
con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_c_complex_mono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 1
A_t0 <- 3^d

p_value <- 1 - pchisq(3.874, df = d)
p_value ## 0.04903952

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 0.1400219

## new_sister_h_c_complex_mono vs. new_sister_c_complex_nonmono
con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_c_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 4
A_t0 <- 3^d

p_value <- 1 - pchisq(40.184, df = d)
p_value ## 3.965272e-08

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 3.211866e-06



# Likelihood ratio test under LG+MEOW(60,20)+H4 model
## Here the best likelihood score is of the new_sister_h_c_complex_mono tree, so we treat it as the alternative topology in all cases
## We compare this topology with all the other topologies

con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_h_c_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 3
A_t0 <- 3^d

p_value <- 1 - pchisq(0.232, df = d)
p_value ## 0.9722652

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 1

## new_sister_h_c_complex_mono vs. new_sister_h_complex_mono
con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_h_complex_mono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 2
A_t0 <- 3^d

p_value <- 1 - pchisq(15.808, df = d)
p_value ## 0.0003692635

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 0.003318467

## new_sister_h_c_complex_mono vs. new_sister_h_complex_nonmono
con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_h_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 4
A_t0 <- 3^d

p_value <- 1 - pchisq(25.95, df = d)
p_value ## 3.238776e-05

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 0.002620013

## new_sister_h_c_complex_mono vs. new_sister_c_complex_mono
con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_c_complex_mono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 1
A_t0 <- 3^d

p_value <- 1 - pchisq(16.14, df = d)
p_value ## 5.882835e-05

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 0.0001764747

## new_sister_h_c_complex_mono vs. new_sister_c_complex_nonmono
con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_c_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 4
A_t0 <- 3^d

p_value <- 1 - pchisq(49.852, df = d)
p_value ## 3.877137e-10

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 3.140481e-08



# Likelihood ratio test under LG+MEOW(60,20)+H6 model
## Here the best likelihood score is of the new_sister_h_complex_mono tree, so we treat it as the alternative topology in all cases
## We compare this topology with all the other topologies

con_tree <- consensus(new_sister_h_complex_mono, new_sister_h_c_complex_mono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 2
A_t0 <- 3^d

p_value <- 1 - pchisq(36.634, df = d)
p_value ## 1.109246e-08

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 9.983216e-08

## new_sister_h_complex_mono vs. new_sister_c_complex_mono
con_tree <- consensus(new_sister_h_complex_mono, new_sister_c_complex_mono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 2
A_t0 <- 3^d

p_value <- 1 - pchisq(16.91, df = d)
p_value ## 0.0002128336

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 0.001913872

## new_sister_h_complex_mono vs. new_sister_h_complex_nonmono
con_tree <- consensus(new_sister_h_complex_mono, new_sister_h_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 2
A_t0 <- 3^d

p_value <- 1 - pchisq(29.398, df = d)
p_value ## 4.133381e-07

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 3.720037e-06

## new_sister_h_complex_mono vs. new_sister_c_complex_nonmono
con_tree <- consensus(new_sister_h_complex_mono, new_sister_c_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 3
A_t0 <- 3^d

p_value <- 1 - pchisq(30.318, df = d)
p_value ## 1.183029e-06

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 3.194129e-05

## new_sister_h_c_complex_mono vs. new_sister_h_c_complex_nonmono
con_tree <- consensus(new_sister_h_complex_mono, new_sister_h_c_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 3
A_t0 <- 3^d

p_value <- 1 - pchisq(22.052, df = d)
p_value ## 6.3626e-05

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 0.001716482



# Likelihood ratio test under LG+MEOW(60,20)+H8 model
## Here the best likelihood score is of the new_sister_h_c_complex_mono tree, so we treat it as the alternative topology in all cases
## We compare this topology with all the other topologies

con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_h_c_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 3
A_t0 <- 3^d

p_value <- 1 - pchisq(15.702, df = d)
p_value ## 0.001305184

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 0.03464849

## new_sister_h_c_complex_mono vs. new_sister_h_complex_mono
con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_h_complex_mono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 2
A_t0 <- 3^d

p_value <- 1 - pchisq(10.764, df = d)
p_value ## 0.004598615

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 0.04063435

## new_sister_h_c_complex_mono vs. new_sister_h_complex_nonmono
con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_h_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 4
A_t0 <- 3^d

p_value <- 1 - pchisq(39.808, df = d)
p_value ## 4.74277e-08

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 3.841636e-06

## new_sister_h_c_complex_mono vs. new_sister_c_complex_mono
con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_c_complex_mono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 1
A_t0 <- 3^d

p_value <- 1 - pchisq(8.912, df = d)
p_value ## 0.002833028

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 0.008475029

## new_sister_h_c_complex_mono vs. new_sister_c_complex_nonmono
con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_c_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 4
A_t0 <- 3^d

p_value <- 1 - pchisq(47.366, df = d)
p_value ## 1.279386e-09

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 1.036303e-07



# Likelihood ratio test under LG+MEOW(60,20)+G + GFmix model
## Here the best likelihood score is of the new_sister_h_complex_mono tree, so we treat it as the alternative topology in all cases
## We compare this topology with all the other topologies

con_tree <- consensus(new_sister_h_complex_mono, new_sister_h_c_complex_mono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 2
A_t0 <- 3^d

p_value <- 1 - pchisq(68.751599, df = d)
p_value ## 1.221245e-15

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 1.099121e-14

## new_sister_h_complex_mono vs. new_sister_c_complex_mono
con_tree <- consensus(new_sister_h_complex_mono, new_sister_c_complex_mono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 2
A_t0 <- 3^d

p_value <- 1 - pchisq(31.751912, df = d)
p_value ## 1.273972e-07

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 1.146574e-06

## new_sister_h_complex_mono vs. new_sister_h_complex_nonmono
con_tree <- consensus(new_sister_h_complex_mono, new_sister_h_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 2
A_t0 <- 3^d

p_value <- 1 - pchisq(63.60649, df = d)
p_value ## 1.54321e-14

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 1.388889e-13

## new_sister_h_complex_mono vs. new_sister_c_complex_nonmono
con_tree <- consensus(new_sister_h_complex_mono, new_sister_c_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 3
A_t0 <- 3^d

p_value <- 1 - pchisq(129.98509, df = d)
p_value ## 0

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 0

## new_sister_h_c_complex_mono vs. new_sister_h_c_complex_nonmono
con_tree <- consensus(new_sister_h_complex_mono, new_sister_h_c_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 3
A_t0 <- 3^d

p_value <- 1 - pchisq(39.152206, df = d)
p_value ## 1.611463e-08

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 4.350949e-07



# Likelihood ratio test for the statoinary-trimmed alignment under LG+MEOW(60,20)+G model
## Here the best likelihood score is of the new_sister_hc_complex_mono tree, so we treat it as the alternative topology in all cases
## We compare this topology with all the other topologies

con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_h_c_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 3
A_t0 <- 3^d

p_value <- 1 - pchisq(14.496, df = d)
p_value ## 0.002302169

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 0.06033347

## new_sister_h_c_complex_mono vs. new_sister_h_complex_mono
con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_h_complex_mono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 2
A_t0 <- 3^d

p_value <- 1 - pchisq(12.148, df = d)
p_value ## 0.002301947

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 0.02052778

## new_sister_h_c_complex_mono vs. new_sister_h_complex_nonmono
con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_h_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 4
A_t0 <- 3^d

p_value <- 1 - pchisq(29.492, df = d)
p_value ## 6.209603e-06

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 0.0005028529

## new_sister_h_c_complex_mono vs. new_sister_c_complex_mono
con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_c_complex_mono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 1
A_t0 <- 3^d

p_value <- 1 - pchisq(15.852, df = d)
p_value ## 6.849411e-05

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 0.0002054682

## new_sister_h_c_complex_mono vs. new_sister_c_complex_nonmono
con_tree <- consensus(new_sister_h_c_complex_mono, new_sister_c_complex_nonmono, check.labels = TRUE, rooted = TRUE)
d <- Nnode(new_sister_h_c_complex_mono) - Nnode(con_tree)
d ## degree of freedom = 4
A_t0 <- 3^d

p_value <- 1 - pchisq(38.6, df = d)
p_value ## 1.279386e-09

### Bonferroni correction
bon_p_value <- 1 - (1 - p_value)^A_t0
bon_p_value ## 6.824897e-06



