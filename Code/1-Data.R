# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                     LIBRARIES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(sf)
library(sp)
library(magrittr)
library(tidyverse)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                   DOWNLOAD DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The data used to characterize acidification comes from DFO
# We also use the benthic habitats description from Jean-Denis Dutil to get
# bathymetry data to interpolate dissolved oxygen.
# For more information read the repo's README.md document.

# Output location for downloaded data
output <- './Data/RawData'

# Will have to upload the data on zenodo and eventually get the data from SLGO.
# For now, I'm using the data downloaded manually from the website.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                   IMPORT DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# --------------- #
# Omega aragonite #
# --------------- #
# Import, spatial object &
arag <- read.csv('./data/rawdata/pHOx.csv') %>%
        st_as_sf(coords = c('Longitude','Latitude'), crs = 4326) %>%
        st_transform(crs = 32198)

# ---------- #
# Bathymetry #
# ---------- #
# Unzip file
unzip(zipfile = paste0(output, '/megahab.zip'), exdir = output)

# File name
fileName <- dir(output, pattern = 'shp$', full.names = T)

# Import and select bathymetry
bathy <- st_read(fileName) %>%
         select(Bathy_Mean) %>%
         mutate(Bathy_Mean = abs(round(Bathy_Mean))) %>%
         rename(DEPH = Bathy_Mean)

# ---------------- #
# Dissolved Oxygen #
# ---------------- #
# Import, spatial object & transform
o2 <- read.csv('./Data/RawData/Btm_O2_GSL_2013_2017.txt', sep = '\t') %>%
      st_as_sf(coords = c('LOND','LATD'), crs = 4326) %>%
      st_transform(crs = 32198)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 FORMAT DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# --------------- #
# Omega aragonite #
# --------------- #
arag <- arag %>%
        cbind(st_coordinates(.), .) %>% # Add coordinates
        .[!is.na(.$Aragonite), ] # Remove NA aragonite values

# Add depth to arag dataset (compare with pressure)
arag$DEPH <- st_intersects(arag, bathy) %>%
             unlist() %>%
             bathy$DEPH[.]


# ---------- #
# Bathymetry #
# ---------- #
bathy <- st_centroid(bathy) %>%
         st_coordinates() %>%
         cbind(bathy, .)


# ---------------- #
# Dissolved Oxygen #
# ---------------- #
# Divide by years and remove duplicated coordinates
ys <- unique(o2$YEAR) %>% sort()
o2Year <- vector('list', length(ys))

for(i in 1:length(ys)) {
  # Year ID
  id <- o2$YEAR == ys[i]

  # Select only year
  o2Year[[i]] <- o2[id, ]

  # Add coordinates
  o2Year[[i]] <- cbind(o2Year[[i]], st_coordinates(o2Year[[i]]))

  # Sort by coordinates and depth
  o2Year[[i]] <- o2Year[[i]][order(o2Year[[i]]$X, o2Year[[i]]$Y, o2Year[[i]]$DEPH), ]

  # Remove duplicated coordinates
  o2Year[[i]] <- o2Year[[i]][!duplicated(o2Year[[i]][, c('X','Y'), drop = T]), ]

  # Remove zero distance points
  temp <- as(o2Year[[i]], 'Spatial')
  temp <- zerodist(temp)[,1]
  if(length(temp) > 1) o2Year[[i]] <- o2Year[[i]][-temp,]
}

# Data point in 2014 messing up the interpolation, likely because it's at the maximum observed depth in the Dutil data
o2Year[[2]]$DEPH[o2Year[[2]]$DEPH == 521] <- 520

# Even though the interpolation works, the results are somehow very weird.
#I will remove 2014 from the data for now, but this should be further investigated
o2Year[[2]] <- NULL
ys <- ys[-2]

# Name hypoxia list
names(o2Year) <- ys

# Change name
doxy <- o2Year


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  EXPORT DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Export object as .RData
save(arag, file = './data/rawData/aragonite.RData')
save(bathy, file = './data/rawData/bathy.RData')
save(doxy, file = './data/rawData/doxy.RData')
