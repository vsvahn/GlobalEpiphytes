## Global epiphytes 
# Downloading, cleaning GBIF occurrence data, calculating EOO and specimen count metrics

#### load packages ----

library(rgbif)
library(tidyverse)
library(rWCVP)
library(sf)
library(readr)
library(CoordinateCleaner)
library(rCAT)


# fill in your gbif credentials 
user <- '' # gbif username 
pwd <- '' # gbif password
email <- '' # your email 


###############################################################################


#### get GBIF data  ----

# get WCVP data from wcvp_summary script
count_wcvp <- count_wcvp 
# filter tp species in < 5 tdwg regions
count_wcvp_4 <- count_wcvp %>% filter(count_wcvp$tdwg_regions < 5) 

# find GBIF usageKeys for names
names <- count_wcvp_4$name_author
queries <- list()
for (i in seq_along(names)) {
  queries[[i]] <- name_backbone(name=names[i])
}
name_matching <- bind_rows(queries) 
name_matching <- name_matching %>%
  select(c(usageKey, verbatim_name, scientificName)) %>%
  left_join(count_wcvp_4[,c('taxon_name','name_author')], by = c('verbatim_name'='name_author')) %>%
  rename(tdwgFullName = verbatim_name, tdwgName = taxon_name, gbifName = scientificName)

# save list of gbif keys as csv file
write_csv(name_matching, 'compare_tdwg_gbif.csv')

#make vector of keys
keys <- as.vector(na.omit(name_matching$usageKey))

# download gbif records using taxon keys
#note: for large downloads, split into subsets

occ_download(
  pred_in("taxonKey", keys),
  pred_in("basisOfRecord", 'PRESERVED_SPECIMEN'),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email)

occ_download_wait('-')
occs <- occ_download_get('') %>%
  occ_download_import()


###############################################################################


#### cleaning gbif data ----

# add tdwg versions of names
occ_data <- occs %>% 
  merge(name_matching, by.x='speciesKey', by.y = 'usageKey') %>%
  select(-c(tdwgFullName, gbifName)) %>%
  filter(!is.na(speciesKey)) #removes erroneous records (e.g. genus-level)

#split data into georeferenced and non-georeferenced
geo_occs <- occ_data %>% filter(!is.na(decimalLatitude) & !is.na(decimalLongitude)) %>%
  filter(decimalLatitude != '' & decimalLongitude != '') 
non_occs <- occ_data %>% filter(!(gbifID %in% geo_occs$gbifID))


## clean non-georeferenced records

# remove records outside native range
#read in country mapping, distributions
countries <- read_csv('country_mapping.csv', na="") %>%
  select(c(LEVEL3_COD, ISO_code))
# get occupied countries per sp.
dist <- dist #get dist from wcvp_summary
dist_unq <- merge(x=dist, y = countries, by.x = 'area_code_l3', by.y = 'LEVEL3_COD') %>%
  distinct(taxon_name, ISO_code)
# filter occs outside native range
filtered_non_occs <- non_occs %>%
  semi_join(dist_unq, by = c("species"='taxon_name', 'countryCode'="ISO_code"))

# remove duplicated locality, or year & stateProvince 
filtered_non_occs <- filtered_non_occs %>%
  group_by(tdwgName) %>%
  filter( is.na(locality) | locality =='' | !duplicated(locality)) %>%
  filter((is.na(stateProvince) | is.na(year)) | year== '' | stateProvince == '' | !duplicated(paste(stateProvince, year))) %>%
  ungroup() %>%
  mutate(rounddecimalLatitude = decimalLatitude, rounddecimalLongitude = decimalLongitude) # to merge with georef data later

## clean georeferenced records

# remove duplicated georeferenced records (to 3 decimal places) & coord uncertainty >100 km
filtered_geo_occs <- geo_occs %>% 
  mutate(rounddecimalLatitude = round(as.numeric(decimalLatitude), digits=3),
         rounddecimalLongitude = round(as.numeric(decimalLongitude), digits=3)) %>%
  distinct(species,rounddecimalLatitude, rounddecimalLongitude, .keep_all = TRUE) %>%
  filter(as.numeric(coordinateUncertaintyInMeters) < 100000 | is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters == '')

# apply CoordinateCleaner filters
filtered_geo_occs <- clean_coordinates(x = filtered_geo_occs,
                                       lon = 'rounddecimalLongitude',
                                       lat = "rounddecimalLatitude",
                                       countries = "countryCode",
                                       species = "tdwgName",
                                       tests = c("equal","gbif", "institutions",
                                                 "zeros", "capitals", "centroids"),
                                       value="clean")

## clip to native ranges using rWCVP
# prep data for rWCVP
spatial_geo <- filtered_geo_occs %>% 
  select(c(tdwgName, decimalLatitude, decimalLongitude)) %>% 
  st_as_sf(coords=c(3,2), crs= st_crs(4326)) # make points spatial
# turns occ data into list grouped by species
split_data <- spatial_geo %>% group_by(tdwgName) %>% group_split()
# find native distributions
native_dist <- foreach(i=1:length(split_data), .packages = c("rWCVP", 'rWCVPdata'), .errorhandling = 'pass') %do% {
  native_dist[[i]] <- wcvp_distribution(split_data[[i]]$tdwgName[[1]], taxon_rank='species', introduced=FALSE, extinct=FALSE, location_doubtful=FALSE)
}
# intersect native distributions with occurrences 
sf::sf_use_s2(FALSE)
system.time(for (i in 1:length(native_dist)){
  tryCatch({
    split_data[[i]]$native <- factor(sf::st_intersects(split_data[[i]], st_union(st_buffer(native_dist[[i]], dist=.09)), sparse=FALSE))},
    error=function(e){cat('ERROR:',conditionMessage(e),'\n')})
})
# convert data back to df
unsplit_data <- data.table::rbindlist(split_data, fill=TRUE, use.names=TRUE)
# extract lat and long
coords <- do.call(rbind, st_geometry(unsplit_data$geometry)) %>% as_tibble() %>% setNames(c('lon','lat'))
unsplit_data$decimalLatitude <- coords$lat
unsplit_data$decimalLongitude <- coords$lon
# filter to native occurrences 
filtered_geo_occs <- filtered_geo_occs %>%
  merge(unsplit_data, by=c('tdwgName', 'decimalLatitude', 'decimalLongitude')) %>%
  filter(native==TRUE) %>% select(-c(native, geometry))


# bind cleaned occs back together 
clean_occs <- rbind(filtered_non_occs, filtered_geo_occs)


#### Calculate EOO and specimen count ----


eoo <- batchCon(taxa=filtered_geo_occs$tdwgName, lat=filtered_geo_occs$rounddecimalLatitude, lon = filtered_geo_occs$rounddecimalLongitude, cellsize=2000) %>%
  select(c(taxon, EOOkm2))

count_records <- count(clean_occs, tdwgName, name='recordCount')

#### gather all data into one df ----

ranges <- count_wcvp %>%
  left_join(count_records, by=c('taxon_name' ='tdwgName')) %>%
  left_join(eoo, by=c('taxon_name' = 'taxon')) %>%
  rename(accepted_name = taxon_name)

write_csv(ranges, 'range_sizes.csv')
