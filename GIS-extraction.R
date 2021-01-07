# GIS-extraction
# Used to extract and collate GIS variables for camera trap analysis
# M Fennell
# mitchfen@mail.ubc.ca
# Last updated: Jan 7, 2021

### 0. Setup ####

install.packages(c("raster","osmdata","sf","bcmaps"))

library(raster)
library(osmdata)
library(sf)
library(dplyr)
library(ggplot2)
library(bcmaps)

getwd()
setwd("D:/Mitch/Cathedral/3.Data/3.2 GIS")

### 1. Build workspace ####

lat_range = c(49.000608, 49.154864)
long_range = c(-120.347670, -120.007913) 

CATH_crs = CRS("+init=EPSG:32610")

convert_coords = function(lat,long, from = CRS("+init=epsg:4326"), to) {
  data = data.frame(long=long, lat=lat)
  coordinates(data) <- ~ long+lat
  proj4string(data) = from
  #Convert to coordinate system specified by EPSG code
  xy = data.frame(sp::spTransform(data, to))
  colnames(xy) = c("x","y")
  return(unlist(xy))
}

utm_bbox = convert_coords(lat = lat_range, long=long_range, to = crs(CATH_crs))
utm_bbox

extent_CATH = extent(utm_bbox[1], utm_bbox[2], utm_bbox[3], utm_bbox[4])

### 2. Extract OSM data ####

osm_bbox = c(long_range[1],lat_range[1], long_range[2],lat_range[2])

CATH_highway = opq(osm_bbox) %>%
  add_osm_feature("highway") %>%
  osmdata_sf()
CATH_highway

# plot to preview
CATH_lines = st_transform(CATH_highway$osm_lines, crs=crs(CATH_crs))

ggplot(CATH_lines, aes(color = osm_id)) +
  geom_sf() +
  theme(legend.position = "none")

# Look good? Now filter into desired categories
CATH_trails = CATH_lines %>%
  filter(highway %in% c("path","bridleway","footway"))

ggplot(CATH_trails, aes(color = osm_id)) +
  geom_sf() +
  theme(legend.position = "none")

# Save as shapefile
st_write(CATH_trails, paste0(getwd(),"/","CATH_trails.shp"))

### 3. Extract elevation for points ####
# Currently using 25m 1:250000 BC Gov DEM

CATH_raster<- cded_raster
