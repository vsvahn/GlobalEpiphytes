# Regressions for Global epiphytes

# load in libraries
library(ape)
library(phylolm)
library(tidyverse)
library(readr)

#load in range size data
ranges <- read_csv('range_sizes.csv') %>%
  mutate(lifeform = lifeform %>% fct_relevel('non-epiphyte'))

# download phylogenetic trees from https://doi.org/10.5281/zenodo.7600341
phylo <- read.tree('Angiosperms_100trees.tre')

# one phylogenetic tree for pruning datasets
phylo <- phylo_list[[1]]

#### Botanical countries analysis ----

#find families with 10+ epiphytes
epi.fams <- ranges %>%
  filter(lifeform=='epiphyte') %>%
  group_by(family) %>%
  dplyr::summarise(n = n()) %>%
  filter(n>9)  

###split data into groups to be analyzed
ranges.list <- ranges %>%
  filter(family %in% epi.fams$family) %>%
  group_by(family) %>%
  group_split 
ranges.list[[25]] <- ranges #manually add all-angiosperm dataframe

# function to run GLM regression for botanical country datasets
ranges.lm <- function(range_data){
  (glm.count <- glm(tdwg_regions ~ lifeform, data=range_data, family='quasipoisson'))
  sum <- summary(glm.count)
  sum$data <- range_data$family[1]
  return(sum)
}

#apply function
ranges.lm.results <- lapply(ranges.list, ranges.lm)
ranges.lm.results[[25]]$data <- 'Angiosperms'


#### EOO analysis  ----

## make EOO dataframe
eoo <- ranges %>% drop_na(EOOkm2) %>% filter(EOOkm2 > 0) %>% # remove species with no EOO
  mutate(logEOO= log(EOOkm2)) # get log of EOO
eoo <- as.data.frame(eoo)
# prune tree for phylogenetic regression
drop_1 <- phylo$tip.label[!phylo$tip.label %in% eoo$tip.label]
tree_1 <- drop.tip(phylo, as.character(drop_1))
# assign unique names to node labels - solve problem of duplicated nodes (null)
tree_1 <- makeLabel(tree_1)
# match row names with the tip labels in the tree
# sort trait data and phylogeny in the same species order, prune dataset to tree
eoo <- as.data.frame(eoo)
rownames(eoo) <- eoo$tip.label
eoo <- eoo[match(tree_1$tip.label,rownames(eoo)),]

# split EOO data by groups to be analyzed
eoo.list <- eoo %>%
  filter(family %in% epi.fams$family) %>%
  group_by(family) %>%
  group_split 
eoo.list[[25]] <- eoo #manually add all-angiosperm dataframe


## ordinary regression

# function to run OLS for EOO datasets
eoo.lm <- function(range_data){
  (lm.eoo <- glm(logEOO ~ lifeform, data=range_data))
  sum <- summary(lm.eoo)
  sum$family <- range_data$family[1]
  return(sum)
}

#apply function
eoo.lm.results <- lapply(eoo.list, eoo.lm)
eoo.lm.results[[25]]$family <- 'Angiosperms'


## PGLS regression 

# test phylogenetic regression with 1 tree
#(pgls.eoo <- phylolm(logEOO ~ lifeform, data=eoo,phy=tree_1,model='lambda'))


# function that calculates mean and sd of regression values for 100 trees
eoo.loop <- function(x) {
  #create empty dataframe
  eoo.trees <- data.frame(matrix(ncol = 9, nrow = 100))
  colnames(eoo.trees) <- c('lambda','intercept_estimate', 'intercept_std_error','intercept_t','intercept_p',
                           'non_epiphyte_estimate', 'non_epiphyte_std_error', 'non_epiphyte_t', 'non_epiphyte_p')
  #prepare eoo df for phylo regression
  x$tip.label <- sub(" ","_", x$accepted_name)
  x <- as.data.frame(x)
  # calculate regression info for each tree, save to dataframe
  for( i in 1:length(phylo_list)){
    phylo <- phylo_list[[i]]
    phylo.spp <- sub("_"," ", phylo$tip.label)
    drop_1 <- phylo$tip.label[!phylo$tip.label %in% x$tip.label]
    tree_1 <- drop.tip(phylo, as.character(drop_1))
    tree_1 <- makeLabel(tree_1)
    rownames(x) <- x$tip.label
    x <- x[match(tree_1$tip.label,rownames(x)),]
    pgls <- phylolm(logEOO ~ lifeform, data=x,phy=tree_1,model='lambda')
    eoo.trees$lambda[i] <- pgls[[3]]
    eoo.trees$intercept_estimate[i] <- summary(pgls)$coef[[1]]
    eoo.trees$intercept_std_error[i] <- summary(pgls)$coef[[3]]
    eoo.trees$intercept_t[i] <- summary(pgls)$coef[[5]]
    eoo.trees$intercept_p[i] <- summary(pgls)$coef[[7]]
    eoo.trees$epiphyte_estimate[i] <- summary(pgls)$coef[[2]]
    eoo.trees$epiphyte_std_error[i] <- summary(pgls)$coef[[4]]
    eoo.trees$epiphyte_t[i] <- summary(pgls)$coef[[6]]
    eoo.trees$epiphyte_p[i] <- summary(pgls)$coef[[8]]
    eoo.trees <- as.data.frame(eoo.trees)
    
  }
  
  # calculate mean, sd of regression values for each tree
  eoo.sum <- eoo.trees %>% summarise_if(is.numeric, list(mean,sd))
  # add family
  eoo.sum$family <- x$family[1]
  return(eoo.sum)
}

#apply function 
eoo.results <- lapply(eoo.list,eoo.loop)
#bind list elements together
eoo.results.df <- data.table::rbindlist(eoo.results, fill=TRUE, use.names=TRUE)
names(eoo.results.df) <- gsub(x = names(eoo.results.df), pattern = "fn1", replacement = "mean")  
names(eoo.results.df) <- gsub(x = names(eoo.results.df), pattern = "fn2", replacement = "sd")  
#make family correct 
eoo.results.df$family[25] <- 'Angiosperms'

#### Specimen count analysis ----

spec <- ranges %>% drop_na(recordCount) %>%  # remove species with no records
  mutate(logRecord= log(recordCount)) #%>% #get log of EOO
# prune tree for phylogenetic regression
drop_3 <- phylo$tip.label[!phylo$tip.label %in% spec$tip.label]
tree_3 <- drop.tip(phylo, as.character(drop_3))
# assign unique names to node labels - solve problem of duplicated nodes (null)
tree_3 <- makeLabel(tree_3)
# match row names with the tip labels in the tree
# sort trait data and phylogeny in the same species order
spec <- as.data.frame(spec)
rownames(spec) <- spec$tip.label
spec <- spec[match(tree_3$tip.label,rownames(spec)),]

# split specimen count data by groups to be analyzed
spec.list <- spec %>%
  filter(family %in% epi.fams$family) %>%
  group_by(family) %>%
  group_split 
spec.list[[25]] <- spec  #manually add all-angiosperm dataframe


## ordinary regression

# function to run OLS for specimen count datasets
spec.lm <- function(range_data){
  (lm.spec <- glm(logRecord ~ lifeform, data=range_data))
  sum <- summary(lm.spec)
  sum$data <- range_data$family[1]
  return(sum)
}

#apply function
spec.lm.results <- lapply(spec.list, spec.lm)
spec.lm.results[[25]]$data <- 'Angiosperms'


## PGLS regression

# test phylogenetic regression
#(pgls.spec <- phylolm(logRecord ~ lifeform, data=spec,phy=tree_3,model='lambda'))


#function that calculates mean and sd for 100 trees
spec.loop <- function(x) {
  #create empty dataframe
  spec.trees <- data.frame(matrix(ncol = 9, nrow = 100))
  colnames(spec.trees) <- c('lambda','intercept_estimate', 'intercept_std_error','intercept_t','intercept_p',
                            'non_epiphyte_estimate', 'non_epiphyte_std_error', 'non_epiphyte_t', 'non_epiphyte_p')
  #prepare eoo df for phylo regression
  x$tip.label <- sub(" ","_", x$accepted_name)
  x <- as.data.frame(x)
  # calculate regression info for each tree, save to dataframe
  for( i in 1:length(phylo_list)){
    phylo <- phylo_list[[i]]
    phylo.spp <- sub("_"," ", phylo$tip.label)
    drop_1 <- phylo$tip.label[!phylo$tip.label %in% x$tip.label]
    tree_1 <- drop.tip(phylo, as.character(drop_1))
    tree_1 <- makeLabel(tree_1)
    rownames(x) <- x$tip.label
    x <- x[match(tree_1$tip.label,rownames(x)),]
    pgls <- phylolm(logRecord ~ lifeform, data=x,phy=tree_1,model='lambda')
    spec.trees$lambda[i] <- pgls[[3]]
    spec.trees$intercept_estimate[i] <- summary(pgls)$coef[[1]]
    spec.trees$intercept_std_error[i] <- summary(pgls)$coef[[3]]
    spec.trees$intercept_t[i] <- summary(pgls)$coef[[5]]
    spec.trees$intercept_p[i] <- summary(pgls)$coef[[7]]
    spec.trees$epiphyte_estimate[i] <- summary(pgls)$coef[[2]]
    spec.trees$epiphyte_std_error[i] <- summary(pgls)$coef[[4]]
    spec.trees$epiphyte_t[i] <- summary(pgls)$coef[[6]]
    spec.trees$epiphyte_p[i] <- summary(pgls)$coef[[8]]
    spec.trees <- as.data.frame(spec.trees)
    
  }
  
  # calculate mean, sd of regression values for each tree
  spec.sum <- spec.trees %>% summarise_if(is.numeric, list(mean,sd))
  # add family 
  spec.sum$family <- x$family[1]
  return(spec.sum)
}

#apply function 
spec.results <- lapply(spec.list,spec.loop)
#bind list elements together
spec.results.df <- data.table::rbindlist(spec.results, fill=TRUE, use.names=TRUE)
names(spec.results.df) <- gsub(x = names(spec.results.df), pattern = "fn1", replacement = "mean")  
names(spec.results.df) <- gsub(x = names(spec.results.df), pattern = "fn2", replacement = "sd")  
#make family correct 
spec.results.df$family[25] <- 'Angiosperms'

#### prep dataframes with just tropical species for sensitivity analysis ----

#read in lat data for all sp
lat <- read.csv('01_raw_data/lat.zones.csv')
lat <- lat %>% select(c(taxon_name, abs.lat.mean, latitude.zone))

#match lat to datasets
ranges.lat <- merge(ranges, lat, by.x='accepted_name', by.y='taxon_name')
eoo.lat <- merge(eoo, lat, by.x='accepted_name', by.y='taxon_name')
spec.lat <- merge(spec, lat, by.x='accepted_name', by.y='taxon_name')

#Create tropical versions of datasets 
ranges.trop <- subset(ranges.lat, latitude.zone=='Tropical')
eoo.trop <- subset(eoo.lat, latitude.zone=='Tropical')
spec.trop <- subset(spec.lat, latitude.zone=='Tropical')

#run same analyses as above


#### prep dataframes with just families containing epiphytes for sensitivity analysis ----

# restrict datasets to families with epiphytes
epi <- ranges %>% filter(lifeform == 'epiphyte')
ranges.epi <- ranges %>% filter(family %in% epi$family)
eoo.sens <- eoo %>% filter(family %in% epi$family)
spec.sens <- spec %>% filter(family %in% epi$family)

# run regressions

#specimen count
spec.loop(spec.sens) #phyloLM
lm.eoo <- lm(logEOO ~ lifeform, data=spec.sens) #OLS

eoo.loop(eoo.sens) #phyloLM
lm.eoo <- lm(logEOO ~ lifeform, data=eoo.sens) #OLS




