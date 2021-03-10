###################################################################
# NDVI/EV/LAI/whatever you can get your hands on from MODIStools
library(MODISTools)
bands <- mt_bands(product = "MOD13Q1") #MOD13Q1
head(bands)

#https://lpdaac.usgs.gov/products/mod13q1v006/ ndv/evi
#https://lpdaac.usgs.gov/products/mod15a2hv006/ lai
# The algorithm chooses the best available pixel value from all the acquisitions from the 16 day period. 
#The criteria used is low clouds, low view angle, and the highest NDVI/EVI value.

# Import the station locations!
setwd("D:/Mitch/Cathedral/3.Data/3.2 GIS")

sta<- read.csv("./CATH_master_sites_Mar09.csv")

dates <- mt_dates(product = "MOD13Q1", lat = 49.05440, lon = -120.2136) #MOD13Q1
head(dates)
tail(dates)

tmp <- sta[,which(colnames(sta)%in%c("Site", "Long", "Lat"))]
colnames(tmp) <- c("site_name", "lat", "lon")
#tmp <- tmp[,c(1,3,2)]
str(tmp)


###################################################################
# GET lai 
algar_lai <- mt_batch_subset(product = "MOD15A2H", #MOD13Q1
                             df=tmp,
                             band = "Lai_500m",
                             start = "2015-01-11",
                             end = "2019-11-30",
                             km_lr = 0,
                             km_ab = 0,
                             internal = TRUE)

par(mar=c(5,4,1,1))
par(mfrow=c(1,1))
hist(algar_lai$value)


algar_lai$calendar_date <- strptime(algar_lai$calendar_date, "%Y-%m-%d")
i <- 1
plot(algar_lai$calendar_date, algar_lai$value, las=1, ylab="LAI", pch=19, type="n")
for(i in 1:length(unique(algar_lai$site)))
{
  points(algar_lai$calendar_date[algar_lai$site==unique(algar_lai$site)[i]], algar_lai$value[algar_lai$site==unique(algar_lai$site)[i]], pch=19, col=rgb(0,0,0,0.1))
}


###################################################################
# GET ndvi
CATH_NDVI_test<- mt_subset(product = "MOD13Q1",
          lat = 49.05440,
          lon = -120.2136,
          band = "250m_16_days_NDVI",
          start = "2019-06-26",
          end = "2019-09-30",
          km_lr = 0,
          km_ab = 0,
          internal = TRUE,
          progress = TRUE)
### Need to multiply values by 0.0001!

CATH_ndvi <- mt_batch_subset(product = "MOD13Q1",
                              df=tmp,
                              band = "250m_16_days_NDVI",
                              start = "2019-06-26",
                              end = "2020-10-03",
                              km_lr = 0.5,
                              km_ab = 0.5,
                              internal = TRUE,
                             out_dir = "D:/Mitch/Cathedral/3.Data/3.2 GIS/MODIS")

par(mar=c(5,4,1,1))
par(mfrow=c(1,1))
hist(CATH_ndvi$value)

# Remove stuff that is below zero... 
CATH_ndvi <- CATH_ndvi[CATH_ndvi$value>0,]

# correct NDVI values by scale factor (0.0001)
CATH_ndvi$value <- CATH_ndvi$value*0.0001

CATH_ndvi$calendar_date <- strptime(CATH_ndvi$calendar_date, "%Y-%m-%d")
i <- 1
plot(CATH_ndvi$calendar_date, CATH_ndvi$value, las=1, ylab="NDVI", pch=19, type="n")
for(i in 1:length(unique(CATH_ndvi$site)))
{
  points(CATH_ndvi$calendar_date[CATH_ndvi$site==unique(CATH_ndvi$site)[i]], CATH_ndvi$value[CATH_ndvi$site==unique(CATH_ndvi$site)[i]], pch=19, col=rgb(0,0,0,0.1))
}

# Load in detection data
setwd("D:/Mitch/Cathedral/3.Data/3.4 Data Analysis/3.4.3 Outputs")
site.week.vars <- read.csv("CATH_weekly_vars_dets.csv", header = T)

# Create site-year-week column for joins
#for (i in 1:length(site.week.vars$Date.Time.Captured)){
#    site.week.vars$Site.Year.Wk[i] <- paste0(site.week.vars$Deployment.Location.ID[i],"-",strftime(site.week.vars$Date.Time.Captured[i], format = "%Y-%V"))
#}

for (i in 1:length(CATH_ndvi$calendar_date)){
  CATH_ndvi$Site.Year.Wk[i] <- paste0(strftime(CATH_ndvi$calendar_date[i], format = "%Y-W%V"))
}

# filter out extra stuff in NDVI df
library(dplyr)
library(tidyr)
CATH_ndvi_short <- select(CATH_ndvi, c("value","Site.Year.Wk", "site"))
CATH_ndvi_short <- rename(CATH_ndvi_short,
                          NDVI = value,
                          Deployment.Location.ID = site,
                          Date = Site.Year.Wk,
                          )

# Join
site.week.vars.ndvi <- left_join(site.week.vars,CATH_ndvi_short,by = c("Deployment.Location.ID", "Date"))
site.week.vars.ndvi <- site.week.vars.ndvi[,c(1:13,44,14:ncol(site.week.vars.ndvi))]
site.week.vars.ndvi <- site.week.vars.ndvi[,-45]


#site.week.vars.ndvi <- site.week.vars.ndvi[,c(1,2,73,3:12,74,13:ncol(site.week.vars.ndvi))]
#site.week.vars.ndvi <- site.week.vars.ndvi[,-c(75,76)]

# Fill weeks with previous NDVI value until a new MODIS pass (16days)
site.week.vars.ndvi <- site.week.vars.ndvi %>% fill(NDVI, .direction=c("down"))
# For the first week, fill with the next value
#site.week.vars.ndvi <- site.week.vars.ndvi %>% fill(NDVI, .direction=c("up"))


setwd("D:/Mitch/Cathedral/3.Data/3.4 Data Analysis/3.4.3 Outputs")
write.csv(site.week.vars.ndvi, "CATH_Detections_Vars_Mar10.csv", row.names = F)


###################################################################
# GET evi 

algar_evi <- mt_batch_subset(product = "MOD13Q1",
                             df=tmp,
                             band = "250m_16_days_EVI",
                             start = "2015-01-11",  # for the purposes of corretly characterizing phenology and DHI, have the start date be as close to Jan 1 as possible even if the cameras didnt didnt operate then
                             end = "2019-11-30", # for the purposes of corretly characterizing phenology and DHI, have the end date be as close to Dec 31 as possible even if the cameras didnt didnt operate then
                             km_lr = 0,
                             km_ab = 0,
                             internal = TRUE)

hist(algar_evi$value)

# Remove stuff that is below zero... 
algar_evi <- algar_evi[algar_evi$value>0,]

algar_evi$calendar_date <- strptime(algar_evi$calendar_date, "%Y-%m-%d")
i <- 1
plot(algar_evi$calendar_date, algar_evi$value, las=1, ylab="evi", pch=19, type="n")
for(i in 1:length(unique(algar_evi$site)))
{
  points(algar_evi$calendar_date[algar_evi$site==unique(algar_evi$site)[i]], algar_evi$value[algar_evi$site==unique(algar_evi$site)[i]], pch=19, col=rgb(0,0,0,0.1))
}


#plitty plots

par(mfrow=c(2,1))
plot(algar_evi$calendar_date, algar_evi$value, las=1, ylab="evi", pch=19, type="n")
for(i in 1:length(unique(algar_evi$site)))
{
  points(algar_evi$calendar_date[algar_evi$site==unique(algar_evi$site)[i]], algar_evi$value[algar_evi$site==unique(algar_evi$site)[i]], pch=19, col=rgb(0,0,0,0.1))
}

i <- 1
plot(algar_ndvi$calendar_date, algar_ndvi$value, las=1, ylab="NDVI", pch=19, type="n")
for(i in 1:length(unique(algar_ndvi$site)))
{
  points(algar_ndvi$calendar_date[algar_ndvi$site==unique(algar_ndvi$site)[i]], algar_ndvi$value[algar_ndvi$site==unique(algar_ndvi$site)[i]], pch=19, col=rgb(0,0,0,0.1))
}
