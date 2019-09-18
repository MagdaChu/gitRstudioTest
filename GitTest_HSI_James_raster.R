Sys.setenv(LANG = "en")
rm(list=ls())
setwd("C:/Users/mec21/Documents/SMRU/Data/HSI_James/")

require(raster)
require(landscapeR)
library(sp)

# data are as csv file

da <- read.csv("Pv xyz predicted distributions_dist_included.csv")
da2 <- read.csv("Pv xyz predicted distributions no distance.csv")
land <- raster("C:/Users/mec21/Documents/SMRU/IBM/Model/Input/land.asc") # for East Coast

head(da)
head(da2)
DistFor <- da[,c(2,3,13)] # predictions based only on foraging points
DistAll <- da[,c(2,3,14)] # predictions based on all points (ask James if inlcuding haul-out)
NoDistFor <- da2[,c(2,3,13)] # predictions based only on foraging points, no distance to land included
NoDistAll <- da2[,c(2,3,14)] # predictions based on all points, no distance included (ask James if inlcuding haul-out)

coordinates(DistFor) <- ~x+y              
gridded(DistFor) <- TRUE
distFor <- raster(DistFor)
plot(distFor)
crs(distFor) <- "+proj=utm +zone=30 +datum=WGS84"

coordinates(DistAll) <- ~x+y              
gridded(DistAll) <- TRUE
distall <- raster(DistAll)
plot(distall)
crs(distall) <- "+proj=utm +zone=30 +datum=WGS84"

coordinates(NoDistFor) <- ~x+y              
gridded(NoDistFor) <- TRUE
NoDistFor <- raster(NoDistFor)
plot(NoDistFor)
crs(NoDistFor) <- "+proj=utm +zone=30 +datum=WGS84"

coordinates(NoDistAll) <- ~x+y              
gridded(NoDistAll) <- TRUE
NoDistall <- raster(NoDistAll)
plot(NoDistall)
crs(NoDistall) <- "+proj=utm +zone=30 +datum=WGS84"

# changing origin, and clipping to model domain, no need to change resolution as it is the same
# scaling from 0 to 1
# some of the land and James pixels do not overlap along the coastline hence values larger than 2. I change them all to land (so 2)

proj4string(land)=proj4string(distFor)
origin(distFor) <- origin(land)
distForModel <- crop(distFor, extent(land))
values(distForModel)[!is.na(values(distForModel))] <- round(distForModel[!is.na(distForModel)]/max(values(distForModel), na.rm=T),2)
plot(distForModel)
plot(land)
distForModel2 <- sum(land,distForModel, na.rm=T)
plot(distForModel2)
values(distForModel2)[values(distForModel2) >= 2] <- 2
plot(distForModel2)
distForModel2

proj4string(land)=proj4string(distall)
origin(distall) <- origin(land)
distAllModel <- crop(distall, extent(land))
values(distAllModel)[!is.na(values(distAllModel))] <- round(distAllModel[!is.na(distAllModel)]/max(values(distAllModel), na.rm=T),2)
plot(distAllModel)
plot(land)
distAllModel2 <- sum(land,distAllModel, na.rm=T)
plot(distAllModel2)
values(distAllModel2)[values(distAllModel2) >= 2] <- 2
plot(distAllModel2)
distAllModel2

proj4string(land)=proj4string(NoDistFor)
origin(NoDistFor) <- origin(land)
NoDistForModel <- crop(NoDistFor, extent(land))
values(NoDistForModel)[!is.na(values(NoDistForModel))] <- round(NoDistForModel[!is.na(NoDistForModel)]/max(values(NoDistForModel), na.rm=T),2)
plot(NoDistForModel)
plot(land)
NoDistForModel2 <- sum(land,NoDistForModel, na.rm=T)
plot(NoDistForModel2)
values(NoDistForModel2)[values(NoDistForModel2) >= 2] <- 2
plot(NoDistForModel2)
NoDistForModel2

proj4string(land)=proj4string(NoDistall)
origin(NoDistall) <- origin(land)
NoDistAllModel <- crop(NoDistall, extent(land))
values(NoDistAllModel)[!is.na(values(NoDistAllModel))] <- round(NoDistAllModel[!is.na(NoDistAllModel)]/max(values(NoDistAllModel), na.rm=T),2)
NoDistAllModelNA <- NoDistAllModel
values(NoDistAllModelNA)[is.na(values(NoDistAllModelNA))] <- 10
plot(NoDistAllModelNA) # so NA is land plus a lot of area in the center

par(mfrow=c(2,2))
plot(NoDistAllModel)
plot(land)
NoDistAllModel2 <- sum(land,NoDistAllModel, na.rm=T)
plot(NoDistAllModel2)
values(NoDistAllModel2)[values(NoDistAllModel2) >= 2] <- 2
plot(NoDistAllModel2)
NoDistAllModel2
summary(values(NoDistAllModel2))
NoDistAllModelZeros <- NoDistAllModel2
values(NoDistAllModelZeros)[values(NoDistAllModelZeros) == 0] <- 10
plot(NoDistAllModelZeros)


# James model does not reach small bays so I am going to calculate mean
# of adjacent cells
# I am going to turn zeros into NAs and than run running window to get mean
# I also have to turn land into zero so the land=2 does not cause large values next to shore

NoDistAllModelExtrap <- NoDistAllModel2
values(NoDistAllModelExtrap)[values(NoDistAllModelExtrap) == 0] <- NA
values(NoDistAllModelExtrap)[values(NoDistAllModelExtrap) == 2] <- 0
plot(NoDistAllModelExtrap)

# RUNNING MOVING WINDOW

NoDistAllModelExtrap2 <- focal(NoDistAllModelExtrap, matrix(c(1,1,1,1,1,1,1,1,1),ncol=3), function (i) {mean(i, na.rm=T)},pad = T, padValue = 0, NAonly=T)
plot(NoDistAllModelExtrap2)
summary(values(NoDistAllModelExtrap2))
summary(values(NoDistAllModelExtrap))

# now I turn the one which were NAs to 0.01 and zero back to 2 (land)

values(NoDistAllModelExtrap2)[values(NoDistAllModelExtrap) == 0] <- 2
values(NoDistAllModelExtrap2)[is.na(values(NoDistAllModelExtrap2))] <- 0.01

plot(NoDistAllModelExtrap2)
summary(values(NoDistAllModelExtrap2))
summary(values(NoDistAllModelExtrap2)[values(NoDistAllModelExtrap2)!=2])

par(mfrow=c(2,2))
plot(distForModel2)
plot(distAllModel2)
plot(NoDistForModel2)
plot(NoDistAllModel2)


# for East Coast
writeRaster(distForModel2,"HSI_James_dist_foragingOnly.asc",format="ascii",overwrite=TRUE)
writeRaster(distForModel2,"C:/Users/mec21/Documents/SMRU/IBM/Model/Input/HSI_James_dist_foragingOnly.asc",format="ascii",overwrite=TRUE)
writeRaster(distAllModel2,"HSI_James_dist_AllPoints.asc",format="ascii",overwrite=TRUE)
writeRaster(distAllModel2,"C:/Users/mec21/Documents/SMRU/IBM/Model/Input/HSI_James_dist_AllPoints.asc",format="ascii",overwrite=TRUE)
writeRaster(NoDistForModel2,"HSI_James_Nodist_foragingOnly.asc",format="ascii",overwrite=TRUE)
writeRaster(NoDistForModel2,"C:/Users/mec21/Documents/SMRU/IBM/Model/Input/HSI_James_Nodist_foragingOnly.asc",format="ascii",overwrite=TRUE)
writeRaster(NoDistAllModel2,"HSI_James_Nodist_AllPoints.asc",format="ascii",overwrite=TRUE)
writeRaster(NoDistAllModel2,"C:/Users/mec21/Documents/SMRU/IBM/Model/Input/HSI_James_Nodist_AllPoints.asc",format="ascii",overwrite=TRUE)
writeRaster(NoDistAllModelExtrap2,"C:/Users/mec21/Documents/SMRU/IBM/Model/Input/HSI_James_Nodist_AllPoints_Extrap_EastCoast.asc",format="ascii",overwrite=TRUE)

# making histogram with prportion of each hsi

hsi <- raster("C:/Users/mec21/Documents/SMRU/IBM/Model/Input/HSI_James_Nodist_AllPoints_Extrap_EastCoast.asc")
plot(hsi)
noland <- as.data.frame(hsi[values(hsi)<=1])
colnames(noland) <- c("hsi")
require(ggplot2)

ggplot(noland, aes(hsi)) + 
  geom_histogram(aes(y=c(..count../sum(..count..)*100)), binwidth = 0.01,show.legend=F)+
  xlab("Habitat suitability index (HSI)") +
  ylab("%")
  


# for Wash
writeRaster(distForModel2,"HSI_James_dist_foragingOnly_Wash.asc",format="ascii",overwrite=TRUE)
writeRaster(distForModel2,"C:/Users/mec21/Documents/SMRU/IBM/Model/Input/HSI_James_dist_foragingOnly_Wash.asc",format="ascii",overwrite=TRUE)
writeRaster(distAllModel2,"HSI_James_dist_AllPoints_Wash.asc",format="ascii",overwrite=TRUE)
writeRaster(distAllModel2,"C:/Users/mec21/Documents/SMRU/IBM/Model/Input/HSI_James_dist_AllPoints_Wash.asc",format="ascii",overwrite=TRUE)
writeRaster(NoDistForModel2,"HSI_James_Nodist_foragingOnly_Wash.asc",format="ascii",overwrite=TRUE)
writeRaster(NoDistForModel2,"C:/Users/mec21/Documents/SMRU/IBM/Model/Input/HSI_James_Nodist_foragingOnly_Wash.asc",format="ascii",overwrite=TRUE)
writeRaster(NoDistAllModel2,"HSI_James_Nodist_AllPoints_Wash.asc",format="ascii",overwrite=TRUE)
writeRaster(NoDistAllModel2,"C:/Users/mec21/Documents/SMRU/IBM/Model/Input/HSI_James_Nodist_AllPoints_Wash.asc",format="ascii",overwrite=TRUE)

# the model does reach far offshore so i am going to turn all zeros and NAs into 0.01 

distAllModel2_no0 <- distAllModel2
values(distAllModel2_no0)[values(distAllModel2_no0) == 0] <- 0.01
values(distAllModel2_no0)[is.na(values(distAllModel2_no0))] <- 0.01
plot(distAllModel2_no0)

NodistAllModel2_no0 <- NoDistAllModel2
values(NodistAllModel2_no0)[values(NodistAllModel2_no0) == 0] <- 0.01
values(NodistAllModel2_no0)[is.na(values(NodistAllModel2_no0))] <- 0.01
plot(NodistAllModel2_no0)

NodistForModel2_no0 <- NoDistForModel2
values(NodistForModel2_no0)[values(NodistForModel2_no0) == 0] <- 0.01
values(NodistForModel2_no0)[is.na(values(NodistForModel2_no0))] <- 0.01
plot(NodistForModel2_no0)





#East Coast
writeRaster(distAllModel2_no0,"HSI_James_dist_AllPoints_no0.asc",format="ascii",overwrite=TRUE)
writeRaster(distAllModel2_no0,"C:/Users/mec21/Documents/SMRU/IBM/Model/Input/HSI_James_dist_AllPoints_no0.asc",format="ascii",overwrite=TRUE)

