## Global epiphytes paper
# Creating a dataset with count of  native tdwg occurrences for all angiosperms

## load required packages ----
library(tidyverse)
library(readr)

## Read in data ----

# Download WCVP at https://powo.science.kew.org/

# WCVP names 
wcvp_names <- read.table("wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8") 
#WCVP distributions
wcvp_dist <- read.table("wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8") 

# Epilist
# download from supp. info in https://doi.org/10.1002/ecy.3326
epilist <- read.csv('EpiList Final revised.csv',  sep=";") 

## lifeform crosswalk from Humphreys et al. 2019
wcvp_lifeforms <- read_csv(file = 'lifeform_mapping.csv')
#add row for 'Parasitic', 'Unknown'
wcvp_lifeforms <- wcvp_lifeforms %>% add_row(lifeform_description = "Parasitic", humphreys_lifeform = "herbaceous perennial")
wcvp_lifeforms <- wcvp_lifeforms %>% add_row(lifeform_description = "", humphreys_lifeform = "Unknown")

######################################################################### 

## Filter species and assign lifeform ----

# add lifeform matching 
hum <- left_join(wcvp_names,wcvp_lifeforms,by="lifeform_description")

# filter out non-accepted sp and hybrids
hum <- hum %>% mutate(epiphyte = ifelse(humphreys_lifeform=='epiphyte', 'epiphyte', 'non-epiphyte')) %>%
  filter(taxon_status == 'Accepted' & species_hybrid=='' & genus_hybrid=='' & taxon_rank=='Species')

# label hemiepiphytes as non-epiphytes
hemi <- epilist %>% filter(Hemi == 'H' )
hum <- hum %>% mutate(epiphyte = ifelse(taxon_name %in% hemi$Species, 'non-epiphyte', epiphyte ))

# add epiphytes missing from WCVP using epilist
true_epi <- epilist %>% filter( Hemi != 'H')
hum <- hum %>% mutate(epiphyte = ifelse(taxon_name %in% true_epi$Species, 'epiphyte', epiphyte))

# remove non-angiosperms
hum <- subset(hum,
              !(family %in% c("Lycopodiaceae", "Ophioglossaceae", 
                              "Selaginellaceae","Aspleniaceae","Equisetaceae",
                              "Isoetaceae","Psilotaceae","Marattiaceae","Osmundaceae",
                              "Hymenophyllaceae", "Gleicheniaceae","Dipteridaceae",
                              "Matoniaceae","Schizaeaceae","Marsileaceae",
                              "Salviniaceae","Cyatheaceae","Cystodiaceae",
                              "Lonchitidaceae", "Lindsaeaceae","Saccolomataceae",
                              "Dennstaedtiaceae","Pteridaceae","Aspleniaceae",
                              "Polypodiaceae", 'Zamiaceae', 'Cupressaceae', 'Taxaceae',
                              'Cycadaceae','Ginkgoaceae', 'Ephedraceae', 'Gnetaceae', 
                              'Welwitschiaceae', 'Pinaceae', 'Araucariaceae', 'Podocarpaceae',
                              'Sciadopityaceae')))

table(hum$epiphyte)
######################################################################### 

##  merge with distribution data

# clean distribution data
dist <- wcvp_dist %>% filter(area_code_l3 != '' & introduced==0 & location_doubtful==0) 

# merge with species list
dist <- dist %>%
  merge(hum, by='plant_name_id') %>% 
  mutate(name_author = paste(taxon_name, taxon_authors))

# count number of occupied tdwg level 3 regions (botanical countries) per sp
count_wcvp <- dist %>%
  group_by(taxon_name, name_author, family , epiphyte) %>% 
  dplyr::summarise(tdwg_regions = n(), .groups = 'drop') %>%
  filter(taxon_name != '') #remove pesky erroneous record

