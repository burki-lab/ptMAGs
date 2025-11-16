################################################################################
#               Biogeography analyses based on metagenomic data                #
################################################################################

# Load packages
library("ggplot2")
library("dplyr")
library("tidyr")
library("qs")
library("Biostrings")
library("sf")
library("rnaturalearth")
library("hrbrthemes")
library("ggpubr")

# set working directory
setwd("ptMAGs/02_biogeography/data/")

# Question 1. Are leptophytes found predominantly in the surface or DCM layer? Are they found in the mesopelagic?

## Load file
occurrence <- read.csv("Lepto_occurence.txt", header=TRUE, sep="\t")
length(unique(occurrence$Station)) ## 147 stations

## Filter out filters from uncommon depth layers (IZZ, ZZZ, and MXL)
occurrence <- occurrence %>% filter(Depth %in% c("SUR", "DCM", "MES"))

length(unique(occurrence$Station)) ## 146 stations


## Order filter depth by depth
occurrence$Depth <- factor(occurrence$Depth , levels=c("SUR", "DCM", "MES"))
length(unique(occurrence$Station)) ## 146 stations

## Plot mean_coverage in ocean layers - with data points
## Followed https://stackoverflow.com/questions/75767454/hide-0-values-zero-counts-in-plots-utilizing-geom-count
## and removed points that are 0 in geom_point

## Replace the names of the ptMAGs
occurrence$ID[occurrence$ID == 'REFM_CHLORO_00001'] <- 'Lepto-01'
occurrence$ID[occurrence$ID == 'TARA_CHLORO_00332'] <- 'Lepto-02'
occurrence$ID[occurrence$ID == 'TARA_CHLORO_00478'] <- 'Lepto-03'
occurrence$ID[occurrence$ID == 'TARA_CHLORO_00158'] <- 'Lepto-04'


p1 <- ggplot(occurrence, aes(x=ID, y=Mean_coverage, fill=Depth)) +
  geom_boxplot(position=position_dodge(preserve = "single"), outlier.shape = NA) +
  scale_y_continuous(name="Mean Coverage (log scale)", trans='log10', labels = scales::comma) +
  scale_fill_manual(values=c("#9ecae1", "#3182bd", "#084594")) +
  geom_point(aes(y= ifelse(Mean_coverage == 0, NA, Mean_coverage)), position=position_jitterdodge(), size=0.6, na.rm = TRUE) +
  labs(x = "") +
  theme_ipsum() +
  ggtitle("Leptophytes mean coverage across depths") +
  theme(plot.title = element_text(size = 12, face = "bold"))
p1


## Zoom into most abundant MAG
p2 <- occurrence %>% 
        filter(ID == "Lepto-01") %>%
        filter(Mean_coverage > 0) %>%
        ggplot(aes(x=Depth, y=Mean_coverage, fill=Depth)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(position=position_jitterdodge(), size=0.6) +
        scale_y_continuous(name="Mean Coverage", labels = scales::comma) +
        scale_fill_manual(values=c("#9ecae1", "#3182bd", "#084594")) +
        labs(x = "") +
        theme_ipsum() +
        ggtitle("Lepto-01 mean coverage across depths") +
        theme(plot.title = element_text(size = 12, face = "bold"),
              legend.position = "none")

## Plot together
p <- ggarrange(p1, p2, labels = "auto", ncol = 1)
p

ggsave(p, filename = "depth.pdf", 
       path = "ptMAGs/02_biogeography/results/metagenomes",
       device = cairo_pdf)

### Based on these plots, we can say that leptophytes are present in both the DCM and surface layers (and veryy little in the mesopelagic layers), 
### with a greater tendency to be found in surface layers.


################################################################################

# Question 2. Which size fraction are leptophytes found in?

## We keep filters of interest, removing bacterial size filters, and removing 
## filters with a very broad range

occurrence_size <- 
  occurrence %>% filter(Filter_Size %in% c("0.22-3", "0.8-5", "3-20", "5-20", "20-180"))

## This leaves 2508 filters

## Order filter size in ascending order
occurrence_size$Filter_Size <- factor(occurrence_size$Filter_Size , levels=c("0.22-3", "0.8-5", "3-20", "5-20", "20-180"))

## Plot mean_coverage in size fractions 
p <- occurrence_size %>%
  ggplot(aes(x=Filter_Size, y=Mean_coverage, fill=Filter_Size)) +
  geom_boxplot(data = ~subset(., Mean_coverage > 0), position=position_dodge(preserve = "single"), outlier.shape = NA) +
  scale_x_discrete("Filter_Size", drop = FALSE) +
  scale_y_continuous(name="Mean Coverage", labels = scales::comma) +
  scale_fill_manual(values=c("#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8", "#253494")) +
  geom_point(aes(y= ifelse(Mean_coverage == 0, NA, Mean_coverage)), position=position_jitterdodge(), size=1, na.rm = TRUE) +
  labs(x = "") +
  theme_ipsum() +
  ggtitle("Leptophytes coverage across depths") +
  theme(plot.title = element_text(size = 12, face = "bold"))+
  facet_wrap(~ID, scale="free") + 
  ggtitle("Leptophytes abundance across size fractions")

p


ggsave(p, filename = "size.pdf", 
       path = "ptMAGs/02_biogeography/results/metagenomes",
       device = cairo_pdf)

### Based on this plot, we can say that leptophytes are very small, larger than 1.6 microns and smaller than 5 microns
### They can potentially leak into the 0.2-16 micron size fraction (but barely).
### They are also NOT present in the 20-180 micron size fraction or the 5-20 micron size fraction
### They are however, present in the 0.22-3 micron, 0.8-5 micron, and 3-20 micron size fraction.
### I would therefore place their size to be around 2-4 microns roughly. 

## To check if these results are biased by sampling bias, I checked the number of filters in each size fraction
## 552 0-0.2
## 80 0.1-0.2
## 72 0.2-0.45
## 152 0.2-1.6
## 580 0.22-3
## 84 0.45-0.8
## 32 0.8-20
## 4 0.8-200
## 632 0.8-2000
## 92 0.8-3
## 676 0.8-5
## 4 1.6-20
## 804 180-2000
## 776 20-180
## 16 20-200
## 180 3-20
## 112 3-2000
## 408 5-20


################################################################################

# Question 3. Where are leptophytes found?

## Load station coordinates
coordinates <- read.csv("stations.postions", header=TRUE, sep="\t")

## Change the station field so it is consistent with the occurrence data file.
coordinates$Station <- gsub("TARA_", "", coordinates$Station)
coordinates$Station <- as.numeric(coordinates$Station)

## Merge the dataframes
occurrence_coord <- merge(occurrence, coordinates, by="Station")

length(unique(occurrence_coord$Station)) ## 146 stations

## Filter points where mean coverage = 0
occurrence_coord <- occurrence_coord %>% filter(Mean_coverage > 0)

length(unique(occurrence_coord$Station)) ## 103 stations

## Plot the stations where leptophytes are found
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

theme_set(theme_bw())

## All Leptophytes
p1 <- ggplot(data = world) +
  geom_sf() +
  coord_sf(expand = FALSE) +
  labs(x = "", y= "") +
  ggtitle("All leptophytes") + 
  theme(plot.title = element_text(size=12)) +
  geom_point(data = occurrence_coord, 
             aes(x=Longitude, y=Latitude, size=Mean_coverage),
             color="#e22e34") +
  labs(size = "Mean coverage")

# Lepto-01
p2 <- ggplot(data = world) +
  geom_sf() +
  coord_sf(expand = FALSE) +
  labs(x = "", y= "") +
  ggtitle("Lepto-01") + 
  theme(plot.title = element_text(size=12)) +
  geom_point(data = occurrence_coord %>%
               filter(ID == "Lepto-01"), 
             aes(x=Longitude, y=Latitude, size=Mean_coverage),
             color="#e22e34") +
  labs(size = "Mean coverage")

# Lepto-02
p3 <- ggplot(data = world) +
  geom_sf() +
  coord_sf(expand = FALSE) +
  labs(x = "", y= "") +
  ggtitle("Lepto-02") + 
  theme(plot.title = element_text(size=12)) +
  geom_point(data = occurrence_coord %>%
               filter(ID == "Lepto-02"), 
             aes(x=Longitude, y=Latitude, size=Mean_coverage),
             color="#e22e34") +
  labs(size = "Mean coverage")

## Lepto-03
p4 <- ggplot(data = world) +
  geom_sf() +
  coord_sf(expand = FALSE) +
  labs(x = "", y= "") +
  ggtitle("Lepto-03") + 
  theme(plot.title = element_text(size=12)) +
  geom_point(data = occurrence_coord %>%
               filter(ID == "Lepto-03"), 
             aes(x=Longitude, y=Latitude, size=Mean_coverage),
             color="#e22e34") +
  labs(size = "Mean coverage")

# Lepto-04
p5 <- ggplot(data = world) +
  geom_sf() +
  coord_sf(expand = FALSE) +
  labs(x = "", y= "") +
  ggtitle("Lepto-04") + 
  theme(plot.title = element_text(size=12)) +
  geom_point(data = occurrence_coord %>%
               filter(ID == "Lepto-04"), 
             aes(x=Longitude, y=Latitude, size=Mean_coverage),
             color="#e22e34") +
  labs(size = "Mean coverage")



p <- ggarrange(p2, p3, p4, p5, labels = "auto")
p

ggsave(p, filename = "lepto_occurence.pdf", 
       path = "ptMAGs/02_biogeography/results/metagenomes",
       device = cairo_pdf)

### We can say that leptophytes are present globally but is rare.
### Seem to be relatively more abundant in the arctic.
### Lepto-01 is more abundant than all other leptophytes.



################################################################################

## Summarise mean coverage across stations
### All leptophytes
occurence_stations <- occurrence %>% 
  group_by(Station) %>% 
  summarise(Mean_coverage = sum(Mean_coverage))

length(unique(occurence_stations$Station)) ## 146 stations in total

#### Drop stations with no coverage
occurence_stations <- occurence_stations %>% filter(Mean_coverage > 0)

length(unique(occurence_stations$Station)) ## 103 stations

mean(occurence_stations$Mean_coverage) ## 7.212
median(occurence_stations$Mean_coverage) ## 2.26

### Lepto-01
occurrence_lepto1 <- occurrence %>% 
  filter(ID == "REFM_CHLORO_00001") %>% 
  filter(Mean_coverage > 0)

length(unique(occurrence_lepto1$Station)) ## 20 stations

mean(occurrence_lepto1$Mean_coverage) ## 6.5
median(occurrence_lepto1$Mean_coverage) ## 1.7



