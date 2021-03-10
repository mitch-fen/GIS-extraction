# GIS-extraction
# Used to extract and collate GIS variables for camera trap analysis
# M Fennell
# mitchfen@mail.ubc.ca
# Last updated: Jan 7, 2021

### 0. Setup ####

#install.packages(c("raster","osmdata","sf","bcmaps","spex"))

library(raster)
library(osmdata)
library(sf)
library(spex)
library(dplyr)
library(ggplot2)
library(bcmaps)
library(rgeos)
library(spatialEco)
library(MODISTools)

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

CATH_roads = CATH_lines %>%
  filter(highway %in% c("track","unclassified","service"))

ggplot(CATH_roads, aes(color = osm_id)) + 
  geom_sf() +
  theme(legend.position = "none")

# Save as shapefile
st_write(CATH_trails, paste0(getwd(),"/","CATH_trails.shp"))
st_write(CATH_roads, paste0(getwd(),"/","CATH_roads.shp"))
st_write(CATH_lines, paste0(getwd(),"/","CATH_roads_trails.shp"))

### 3. Extract elevation for points ####
# Currently using 25m 1:250000 BC Gov CDED DEM 
# (https://www2.gov.bc.ca/gov/content/data/geographic-data-services/topographic-data/elevation/digital-elevation-model)

# Download DEM for AOI
extent_CATH_sf <- spex(extent_CATH, crs = CATH_crs)
CATH_raster<- cded_raster(aoi = extent_CATH_sf)

# Load camera location data
CATH_points <- read.csv("Cathedral_Camera_Deployments_August2020.csv")
colnames(CATH_points)[1] <- "Site"

# Convert site locations to spatial points
CATH_points_sp <- SpatialPoints(cbind.data.frame(CATH_points$Long, CATH_points$Lat))

# Extract value from raster at each point
CATH_points_el_sp <- extract(CATH_raster, CATH_points_sp, sp = T)

CATH_points_elev <- as.data.frame(CATH_points_el_sp)
colnames(CATH_points_elev) <- c("Elevation", "Long", "Lat")

# Add elevation data to master DF
CATH_points_master <- full_join(CATH_points,CATH_points_elev)

## Add aspect and slope using DEM ##

# Turn off scientific notation because it's annoying
options(scipen = 999)

aspect <- terrain(CATH_raster, "aspect", "degrees", 8)
CATH_points_aspct <- extract(aspect, CATH_points_sp, sp = T)
CATH_points_aspct <- as.data.frame(CATH_points_aspct)
colnames(CATH_points_aspct) <- c("Aspect", "Long", "Lat")

# Add aspect to master DF
CATH_points_master <- left_join(CATH_points_master, CATH_points_aspct)

slope <- terrain(CATH_raster, "slope", "degrees", 8)
CATH_points_slope <- extract(slope, CATH_points_sp, sp = T)
CATH_points_slope <- as.data.frame(CATH_points_slope)
colnames(CATH_points_slope) <- c("Slope", "Long", "Lat")

# Add slope to master DF
CATH_points_master <- left_join(CATH_points_master, CATH_points_slope)



### 4. Extract distance to water ####
# Using NRCAN CANVEC (https://maps.canada.ca/czs/index-en.html)

CATH_watercourse <- shapefile("D:/Mitch/Cathedral/3.Data/3.2 GIS/CANVEC_Water/watercourse_1.shp")
CATH_waterbody <- shapefile("D:/Mitch/Cathedral/3.Data/3.2 GIS/CANVEC_Water/waterbody_2.shp")
CATH_waterbody_line <- as(CATH_waterbody, "SpatialLinesDataFrame")

# Join waterbodies and watercourses (lakes and rivers/streams)
CATH_water_combi <- raster::union(CATH_watercourse, CATH_waterbody_line)
CATH_water_combi

# Re-project to UTM
CATH_water_combi_proj <- spTransform(CATH_water_combi, CATH_crs)
CATH_water_combi_proj

# Fix up points to match
proj4string(CATH_points_sp) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
CATH_points_sp_proj <- spTransform(CATH_points_sp, CATH_crs)
CATH_points_sp_proj

# Extract distance to water
CATH_dist_h2o <- cbind.data.frame(CATH_points$Site, CATH_points$Long, CATH_points$Lat)
colnames(CATH_dist_h2o) <- c("Site", "Long", "Lat")

for (i in 1:nrow(CATH_points)){
  tmp <- gDistance(CATH_points_sp_proj[i], CATH_water_combi_proj)
  CATH_dist_h2o$d.H2O[i] <- tmp
}

CATH_points_master <- left_join(CATH_points_master, CATH_dist_h2o)

### 5. Extract distance to trails and roads (and both) ####
#trails
CATH_trail_shp <- shapefile("D:/Mitch/Cathedral/3.Data/3.2 GIS/CATH_trails.shp")
plot(CATH_trail_shp)

CATH_dist_trails <- cbind.data.frame(CATH_points$Site, CATH_points$Long, CATH_points$Lat)
colnames(CATH_dist_trails) <- c("Site", "Long", "Lat")

for (i in 1:nrow(CATH_points)){
  tmp <- gDistance(CATH_points_sp_proj[i], CATH_trail_shp)
  CATH_dist_trails$d.TRL[i] <- tmp
}

CATH_points_master <- left_join(CATH_points_master, CATH_dist_trails)

# roads
CATH_road_shp <- shapefile("D:/Mitch/Cathedral/3.Data/3.2 GIS/CATH_roads.shp")
plot(CATH_road_shp)

CATH_dist_roads <- cbind.data.frame(CATH_points$Site, CATH_points$Long, CATH_points$Lat)
colnames(CATH_dist_roads) <- c("Site", "Long", "Lat")

for (i in 1:nrow(CATH_points)){
  tmp <- gDistance(CATH_points_sp_proj[i], CATH_road_shp)
  CATH_dist_roads$d.ROAD[i] <- tmp
}

CATH_points_master <- left_join(CATH_points_master, CATH_dist_roads)

#combined road and trails (linear features)
CATH_road_trail_shp <- shapefile("D:/Mitch/Cathedral/3.Data/3.2 GIS/CATH_roads_trails.shp")
plot(CATH_road_trail_shp)

CATH_dist_roads_trails <- cbind.data.frame(CATH_points$Site, CATH_points$Long, CATH_points$Lat)
colnames(CATH_dist_roads_trails) <- c("Site", "Long", "Lat")

for (i in 1:nrow(CATH_points)){
  tmp <- gDistance(CATH_points_sp_proj[i], CATH_road_trail_shp)
  CATH_dist_roads_trails$d.LIN[i] <- tmp
}

CATH_points_master <- left_join(CATH_points_master, CATH_dist_roads_trails)

### 6. Calculate terrain ruggedness ####
#VRM (Sappington et al 2007)

vrm3 <- vrm(CATH_raster, s=3) #Likely use this one
vrm5 <- vrm(CATH_raster, s=5)

# Extract value from raster at each point
CATH_points_vrm_sp <- extract(vrm3, CATH_points_sp, sp = T)

CATH_points_vrm <- as.data.frame(CATH_points_vrm_sp)
colnames(CATH_points_vrm) <- c("VRM", "Long", "Lat")

# Add VRM data to master DF
CATH_points_master <- left_join(CATH_points_master,CATH_points_vrm)


### ?. Save and export master ####

write.csv(CATH_points_master, "CATH_master_sites_Mar09.csv", row.names = F)

