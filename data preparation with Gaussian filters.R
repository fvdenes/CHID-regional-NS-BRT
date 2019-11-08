library(raster)
library(dismo)
library(rpart)
library(maptools)
library(data.table)
library(rgdal)
library(dplyr)

load("D:/Avian data processed/BAM_data_package_August2019.RData")
#load("D:/Avian data processed/data_package_2016-04-18.Rdata")	
#load("D:/Avian data processed/offsets-v3_2016-04-18.Rdata")
#LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# coordinates(SS) <- c("X", "Y") 
# proj4string(SS) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# SSLCC <- spTransform(SScombo, LCC)
# 
# offl <- data.table(melt(OFF))
# names(offl) <- c("PKEY","SPECIES","logoffset")
# offl$SPECIES <- as.character(offl$SPECIES)
# offl$PKEY <-as.character(offl$PKEY)
# rm(OFF)

#eco<-shapefile("D:/CHID regional NS BRT/NNS_eco.shp")
#eco<-rasterize(eco,NNS)
#writeRaster(eco, filename="NNS_eco", format="GTiff",overwrite=TRUE)
eco <- raster("D:/CHID regional NS BRT/spatial/NNS_eco.tif") 
NNS <- raster("D:/CHID regional NS BRT/spatial/NNS.tif") 
proj4string(NNS)<-LCC

NNS_shp<-readOGR("D:/CHID regional NS BRT/spatial",layer="northernNSmngtunit")
NNS_shp<-spTransform(NNS_shp,LCC)
plot(NNS)

# Load Beaudoin layers (2011 and 2001), crop and mask for study region, save as rasters
b2011 <- list.files("D:/Beaudoin/2011/",pattern="tif$")
setwd("D:/Beaudoin/2011/")
bs2011 <- stack(raster(b2011[1]))
for (i in 2:length(b2011)) { bs2011 <- addLayer(bs2011, raster(b2011[i]))}
names(bs2011) <- gsub("NFI_MODIS250m_2011_kNN_","",names(bs2011))
abs2011 <- crop(bs2011,NNS)
abs2011<-mask(abs2011,NNS)
writeRaster(abs2011, filename="D:/Beaudoin/2011/Processed/NNS/abs2011_250m.grd", format="raster",overwrite=TRUE)
abs2011<-brick("D:/Beaudoin/2011/Processed/NNS/abs2011_250m.grd")

b2001 <- list.files("D:/Beaudoin/2001/",pattern="tif$")
setwd("D:/Beaudoin/2001/")
bs2001 <- stack(raster(b2001[1]))
for (i in 2:length(b2001)) { bs2001 <- addLayer(bs2001, raster(b2001[i]))}
names(bs2001) <- gsub("NFI_MODIS250m_2001_kNN_","",names(bs2001))
abs2001 <- crop(bs2001,NNS)
abs2001<-mask(abs2001,NNS)
writeRaster(abs2001, filename="D:/Beaudoin/2001/Processed/NNS/abs2001_250m.grd", format="raster",overwrite=TRUE)
abs2001<-brick("D:/Beaudoin/2001/Processed/NNS/abs2001_250m.grd")

# obtain weighted sums of neighourhood cells using Gaussian filter with sigma=250, and 750m for Beaudoin and CTI layers, save outputs as rasters
## sigma = 250m
fw250<-focalWeight(x=abs2011,d=250,type="Gauss")
abs2011_Gauss250<-stack(focal(abs2011[[1]],w=fw250,na.rm=TRUE))
names(abs2011_Gauss250)<-names(abs2011)[[1]]
for(i in 2:nlayers(abs2011)){
 abs2011_Gauss250<-addLayer(abs2011_Gauss250,focal(abs2011[[i]],w=fw250,na.rm=TRUE))
 names(abs2011_Gauss250)[i]<-names(abs2011)[[i]]
}
abs2011_Gauss250<-brick(abs2011_Gauss250)
writeRaster(abs2011_Gauss250, filename="D:/Beaudoin/2011/Processed/NNS/abs2011_250_Gauss250m.grd", format="raster",overwrite=TRUE)
abs2011_Gauss250<-brick("D:/Beaudoin/2011/Processed/NNS/abs2011_250_Gauss250m.grd")

## sigma = 750m
fw750<-focalWeight(x=abs2011,d=750,type="Gauss")
abs2011_Gauss750<-brick(focal(abs2011[[1]],w=fw750,na.rm=TRUE))
names(abs2011_Gauss750)<-names(abs2011)[[1]]
for(i in 2:nlayers(abs2011)){
 abs2011_Gauss750<-addLayer(abs2011_Gauss750,focal(abs2011[[i]],w=fw750,na.rm=TRUE))
 names(abs2011_Gauss750)[i]<-names(abs2011)[[i]]
}
abs2011_Gauss750<-brick(abs2011_Gauss750)
writeRaster(abs2011_Gauss750, filename="D:/Beaudoin/2011/Processed/NNS/abs2011_250_Gauss750m.grd", format="raster",overwrite=TRUE)
abs2011_Gauss750<-brick("D:/Beaudoin/2011/Processed/NNS/abs2011_250_Gauss750m.grd")



#Beaudoin 2001
## sigma = 250m
fw250<-focalWeight(x=abs2001,d=250,type="Gauss")
abs2001_Gauss250<-stack(focal(abs2001[[1]],w=fw250,na.rm=TRUE))
names(abs2001_Gauss250)<-names(abs2001)[[1]]
for(i in 2:nlayers(abs2001)){
 abs2001_Gauss250<-addLayer(abs2001_Gauss250,focal(abs2001[[i]],w=fw250,na.rm=TRUE))
 names(abs2001_Gauss250)[i]<-names(abs2001)[[i]]
}
abs2001_Gauss250<-brick(abs2001_Gauss250)
writeRaster(abs2001_Gauss250, filename="D:/Beaudoin/2001/Processed/NNS/abs2001_250_Gauss250m.grd", format="raster",overwrite=TRUE)
#abs2001_Gauss250<-brick("D:/Beaudoin/2001/Processed/NNS/abs2001_250_Gauss250m.grd")

## sigma = 750m
fw750<-focalWeight(x=abs2001,d=750,type="Gauss")
abs2001_Gauss750<-brick(focal(abs2001[[1]],w=fw750,na.rm=TRUE))
names(abs2001_Gauss750)<-names(abs2001)[[1]]
for(i in 2:nlayers(abs2001)){
 abs2001_Gauss750<-addLayer(abs2001_Gauss750,focal(abs2001[[i]],w=fw750,na.rm=TRUE))
 names(abs2001_Gauss750)[i]<-names(abs2001)[[i]]
}
abs2001_Gauss750<-brick(abs2001_Gauss750)
writeRaster(abs2001_Gauss750, filename="D:/Beaudoin/2001/Processed/NNS/abs2001_250_Gauss750m.grd", format="raster",overwrite=TRUE)
#abs2001_Gauss750<-brick("D:/Beaudoin/2001/Processed/NNS/abs2001_250_Gauss750m.grd")


# water table data. Water table data is a factor, so no sense in applying the Gaussian filter smoother. Instead, obtain the modal value over equal weight matrices of size 3x3 (250m GF analog) and 5x5 (750m GF analog) 
wtbl250 <- raster("D:/CHID regional NS BRT/spatial/watertableNNS_LCC.tif")
wtbl250_modal3x3<-focal(wtbl250,w=matrix(1,3,3),na.rm=TRUE, fun=modal)
wtbl250_modal5x5<-focal(wtbl250,w=matrix(1,5,5),na.rm=TRUE, fun=modal)


writeRaster(wtbl250_modal3x3, filename="D:/CHID regional NS BRT/spatial/wtbl250_modal3x3.asc", format="ascii",overwrite=TRUE)
writeRaster(wtbl250_modal5x5, filename="D:/CHID regional NS BRT/spatial/wtbl250_modal5x5.asc", format="ascii",overwrite=TRUE)

# CTI data
CTI<-raster("D:/CTI/CTI_NS.tif")
CTI<-projectRaster(CTI,wtbl250)
CTI<-crop(CTI,NNS)
CTI250<-resample(CTI,abs2011[[1]])
CTI250<-mask(CTI250,abs2011[[1]])

CTI250_Gauss250<-focal(CTI250,w=fw250,na.rm=TRUE)
CTI250_Gauss750<-focal(CTI250,w=fw750,na.rm=TRUE)

writeRaster(CTI250, filename="D:/CHID regional NS BRT/spatial/CTI250.asc", format="ascii",overwrite=TRUE)
CTI250<-raster("D:/CHID regional NS BRT/spatial/CTI250.asc")
writeRaster(CTI250_Gauss250, filename="D:/CHID regional NS BRT/spatial/CTI250_Gauss250.asc", format="ascii",overwrite=TRUE)
CTI250_Gauss250<-raster("D:/CHID regional NS BRT/spatial/CTI250_Gauss250.asc")
writeRaster(CTI250_Gauss750, filename="D:/CHID regional NS BRT/spatial/CTI250_Gauss750.asc", format="ascii",overwrite=TRUE)
CTI250_Gauss750<-raster("D:/CHID regional NS BRT/spatial/CTI250_Gauss750.asc")

# roads (Venter et al)
roadsNNS<-raster("D:/CHID regional NS BRT/spatial/roadsNNS.tif")
roadsNNS<-projectRaster(roadsNNS,wtbl250)
roadsNNS<-crop(roadsNNS,NNS)
roadsNNS<-resample(roadsNNS,abs2011[[1]])
roadsNNS<-mask(roadsNNS,abs2011[[1]])

roadsNNS_Gauss250<-focal(roadsNNS,w=fw250,na.rm=TRUE)
roadsNNS_Gauss750<-focal(roadsNNS,w=fw750,na.rm=TRUE)

# urban/agriculture
urbag <- raster("D:/urbag2011_lcc1/urbag2011_lcc1.tif")
urbagNNS<- crop(urbag,NNS)
urbagNNS<- mask(urbagNNS,abs2011[[1]])

urbagNNS_Gauss250<-focal(urbagNNS,w=fw250,na.rm=TRUE)
urbagNNS_Gauss750<-focal(urbagNNS,w=fw750,na.rm=TRUE)

# water and "lake edge density"
wat <- raster("D:/wat2011_lcc1/wat2011_lcc1.tif")
watNNS<- crop(wat,NNS)
watNNS<- mask(watNNS,abs2011[[1]])

watNNS_Gauss250<-focal(watNNS,w=fw250,na.rm=TRUE)
watNNS_Gauss750<-focal(watNNS,w=fw750,na.rm=TRUE)


# climate data- upload and resample to 250m resolution to match other layers, and attach to abs2011
#climateAW2010 <- list.files("D:/ClimateAdaptWest/baseline19812010/",pattern="asc$")
#setwd("D:/ClimateAdaptWest/baseline19812010/")
#clim2010 <- stack(raster(climateAW2010[1]))
#for (i in 2:length(climateAW2010)) { clim2010 <- addLayer(clim2010, raster(climateAW2010[i]))}
#proj4string(clim2010)<-LCC
#aclim2010 <- crop(clim2010,abs2011)
#aclim2010<-resample(clim2010,abs2011)
#aclim2010<-mask(aclim2010,abs2011$LandCover_NonVeg_v1)
# 
#writeRaster(aclim2010, filename="D:/ClimateAdaptWest/baseline19812010/Processed/NNS/aclim2010.grd", format="raster",overwrite=TRUE)
#aclim2010<-stack("D:/ClimateAdaptWest/baseline19812010/Processed/NNS/aclim2010.grd")

#for(i in 1:length(names(aclim2010))){ 
#  abs2011 <- addLayer(abs2011, aclim2010[[i]])
#  names(abs2011)[nlayers(abs2011)] <- names(aclim2010[[i]])
#}


# put together prediction rasterstack
## need to update names of layers with Gaussian filters first to differentiate them
for(i in 1:nlayers(abs2001_Gauss250)){
  names(abs2001_Gauss250)[i] <- paste(names(abs2001_Gauss250)[i],"_Gauss250",sep="")
  names(abs2001_Gauss750)[i] <- paste(names(abs2001_Gauss750)[i],"_Gauss750",sep="")
  names(abs2011_Gauss250)[i] <- paste(names(abs2011_Gauss250)[i],"_Gauss250",sep="")
  names(abs2011_Gauss750)[i] <- paste(names(abs2011_Gauss750)[i],"_Gauss750",sep="")
}

names(wtbl250) <- "wtbl250"
names(wtbl250_modal3x3) <- "wtbl250_modal3x3"
names(wtbl250_modal5x5) <- "wtbl250_modal5x5"

names(CTI250) <- "CTI250"
names(CTI250_Gauss250) <- "CTI250_Gauss250"
names(CTI250_Gauss750) <- "CTI250_Gauss750"

names(roadsNNS) <- "roadsNNS"
names(roadsNNS_Gauss250) <- "roadsNNS_Gauss250"
names(roadsNNS_Gauss750) <- "roadsNNS_Gauss750"

names(urbagNNS) <- "urbagNNS"
names(urbagNNS_Gauss250) <- "urbagNNS_Gauss250"
names(urbagNNS_Gauss750) <- "urbagNNS_Gauss750"

names(watNNS) <- "watNNS"
names(watNNS_Gauss250) <- "watNNS_Gauss250"
names(watNNS_Gauss750) <- "watNNS_Gauss750"

ARU<-raster(extent(pred_abs_2011),crs=LCC,resolution=res(pred_abs_2011),vals=0)
names(ARU)<-"ARU"

pred_abs_2011<-stack(abs2011,abs2011_Gauss250,abs2011_Gauss750,wtbl250,wtbl250_modal3x3,wtbl250_modal5x5,CTI250,CTI250_Gauss250,CTI250_Gauss750,roadsNNS,roadsNNS_Gauss250,roadsNNS_Gauss750,urbagNNS, urbagNNS_Gauss250, urbagNNS_Gauss750, watNNS, watNNS_Gauss250, watNNS_Gauss750,ARU, quick=TRUE)
names(pred_abs_2011)

writeRaster(pred_abs_2011, filename="D:/CHID regional NS BRT/prediction dataset/abs2011_250m.grd", format="raster",overwrite=TRUE)

pred_abs_2011<-stack("D:/CHID regional NS BRT/prediction dataset/abs2011_250m.grd")

# Extracting data from rasters for surveyed locations (ABSS)
## 2011

spSS<-SpatialPointsDataFrame(coords=SScombo[,2:3],data=SScombo,proj4string = LCC)
NSSS <- crop(spSS,NNS_shp)
NSSS <- NSSS@data
NSSS

dat2011 <- cbind(NSSS, extract(abs2011,as.matrix(cbind(NSSS$X,NSSS$Y)))) # includes Beaudoin layers (250m, no Gaussian filter) 

for(i in 1:nlayers(abs2011_Gauss250)){
  dat2011 <- cbind(dat2011, extract(abs2011_Gauss250[[i]],as.matrix(cbind(NSSS$X,NSSS$Y)))) # includes Beaudoin layers with Gaussian filter sigma=250m
  names(dat2011)[ncol(dat2011)] <- names(abs2011_Gauss250)[i]
}

for(i in 1:nlayers(abs2011_Gauss750)){
  dat2011 <- cbind(dat2011, extract(abs2011_Gauss750[[i]],as.matrix(cbind(NSSS$X,NSSS$Y)))) # includes Beaudoin layers with Gaussian filter sigma=750m
  names(dat2011)[ncol(dat2011)] <- names(abs2011_Gauss750)[i]
}

dat2011<-cbind(dat2011,extract(wtbl250,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes wtbl 250m resolution data
names(dat2011)[ncol(dat2011)] <- "wtbl250"

dat2011<-cbind(dat2011,extract(wtbl250_modal3x3,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes wtbl 250m resolution data, modal of 3x3 cell equal weight matrix
names(dat2011)[ncol(dat2011)] <- "wtbl250_modal3x3"

dat2011<-cbind(dat2011,extract(wtbl250_modal5x5,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes wtbl 250m resolution data, modal of 5x5 cell equal weight matrix
names(dat2011)[ncol(dat2011)] <- "wtbl250_modal5x5"

dat2011<-cbind(dat2011,extract(CTI250,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes CTI 250m resolution data
names(dat2011)[ncol(dat2011)] <- "CTI250"

dat2011<-cbind(dat2011,extract(CTI250_Gauss250,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes CTI 250m resolution data, Gaussian filter sigma=250m
names(dat2011)[ncol(dat2011)] <- "CTI250_Gauss250"

dat2011<-cbind(dat2011,extract(CTI250_Gauss750,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes CTI 250m resolution data, Gaussian filter sigma=750m
names(dat2011)[ncol(dat2011)] <- "CTI250_Gauss750"

dat2011<-cbind(dat2011,extract(roadsNNS,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes roads 250m resolution data
names(dat2011)[ncol(dat2011)] <- "roadsNNS"

dat2011<-cbind(dat2011,extract(roadsNNS_Gauss250,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes roads 250m resolution data, Gaussian filter sigma=250m
names(dat2011)[ncol(dat2011)] <- "roadsNNS_Gauss250"

dat2011<-cbind(dat2011,extract(roadsNNS_Gauss750,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes roads 250m resolution data, Gaussian filter sigma=750m
names(dat2011)[ncol(dat2011)] <- "roadsNNS_Gauss750"

dat2011<-cbind(dat2011,extract(urbagNNS,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes urbag 250m resolution data
names(dat2011)[ncol(dat2011)] <- "urbagNNS"

dat2011<-cbind(dat2011,extract(urbagNNS_Gauss250,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes urbag 250m resolution data, Gaussian filter sigma=250m
names(dat2011)[ncol(dat2011)] <- "urbagNNS_Gauss250"

dat2011<-cbind(dat2011,extract(urbagNNS_Gauss750,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes urbag 250m resolution data, Gaussian filter sigma=750m
names(dat2011)[ncol(dat2011)] <- "urbagNNS_Gauss750"

dat2011<-cbind(dat2011,extract(watNNS,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes wat 250m resolution data
names(dat2011)[ncol(dat2011)] <- "watNNS"

dat2011<-cbind(dat2011,extract(watNNS_Gauss250,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes wat 250m resolution data, Gaussian filter sigma=250m
names(dat2011)[ncol(dat2011)] <- "watNNS_Gauss250"

dat2011<-cbind(dat2011,extract(watNNS_Gauss750,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes wat 250m resolution data, Gaussian filter sigma=750m
names(dat2011)[ncol(dat2011)] <- "watNNS_Gauss750"


### set up weight matrix for SS, and calculate weight values for each row in dat2011
r2 <- abs2011[[1]]
samprast2011 <- rasterize(cbind(dat2011$X,dat2011$Y), r2, field=1)
sampsum25 <- focal(samprast2011, w=matrix(1/25, nc=5, nr=5), na.rm=TRUE)

dat2011 <- cbind(dat2011,extract(sampsum25,as.matrix(cbind(dat2011$X,dat2011$Y))))
names(dat2011)[ncol(dat2011)] <- "sampsum25"
dat2011$wt <- 1/dat2011$sampsum25

dat2011$SS <- as.character(dat2011$SS)
#dat2011$PCODE <- as.character(dat2011$PCODE)
#dat2011<-dat2011[,-c(25:43)] # remove columns with climate data (from avian dataset)

setwd("D:/CHID regional NS BRT/")
write.csv(dat2011,"NNSdat2011.csv")


## 2001
dat2001 <- cbind(NSSS, extract(abs2001,as.matrix(cbind(NSSS$X,NSSS$Y)))) # includes Beaudoin layers (250m, no Gaussian filter)


for(i in 1:nlayers(abs2001_Gauss250)){
  dat2001 <- cbind(dat2001, extract(abs2001_Gauss250[[i]],as.matrix(cbind(NSSS$X,NSSS$Y)))) # includes Beaudoin layers with Gaussian filter sigma=250m
  names(dat2001)[ncol(dat2001)] <- names(abs2001_Gauss250)[i]
}

for(i in 1:nlayers(abs2001_Gauss750)){
  dat2001 <- cbind(dat2001, extract(abs2001_Gauss750[[i]],as.matrix(cbind(NSSS$X,NSSS$Y)))) # includes Beaudoin layers with Gaussian filter sigma=750m
  names(dat2001)[ncol(dat2001)] <- names(abs2001_Gauss750)[i]
}

dat2001<-cbind(dat2001,extract(wtbl250,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes wtbl 250m resolution data
names(dat2001)[ncol(dat2001)] <- "wtbl250"

dat2001<-cbind(dat2001,extract(wtbl250_modal3x3,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes wtbl 250m resolution data, modal of 3x3 cell equal weight matrix
names(dat2001)[ncol(dat2001)] <- "wtbl250_modal3x3"

dat2001<-cbind(dat2001,extract(wtbl250_modal5x5,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes wtbl 250m resolution data, modal of 5x5 cell equal weight matrix
names(dat2001)[ncol(dat2001)] <- "wtbl250_modal5x5"

dat2001<-cbind(dat2001,extract(CTI250,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes CTI 250m resolution data
names(dat2001)[ncol(dat2001)] <- "CTI250"

dat2001<-cbind(dat2001,extract(CTI250_Gauss250,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes CTI 250m resolution data, Gaussian filter sigma=250m
names(dat2001)[ncol(dat2001)] <- "CTI250_Gauss250"

dat2001<-cbind(dat2001,extract(CTI250_Gauss750,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes CTI 250m resolution data, Gaussian filter sigma=750m
names(dat2001)[ncol(dat2001)] <- "CTI250_Gauss750"

dat2001<-cbind(dat2001,extract(roadsNNS,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes roads 250m resolution data
names(dat2001)[ncol(dat2001)] <- "roadsNNS"

dat2001<-cbind(dat2001,extract(roadsNNS_Gauss250,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes roads 250m resolution data, Gaussian filter sigma=250m
names(dat2001)[ncol(dat2001)] <- "roadsNNS_Gauss250"

dat2001<-cbind(dat2001,extract(roadsNNS_Gauss750,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes roads 250m resolution data, Gaussian filter sigma=750m
names(dat2001)[ncol(dat2001)] <- "roadsNNS_Gauss750"

dat2001<-cbind(dat2001,extract(urbagNNS,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes urbag 250m resolution data
names(dat2001)[ncol(dat2001)] <- "urbagNNS"

dat2001<-cbind(dat2001,extract(urbagNNS_Gauss250,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes urbag 250m resolution data, Gaussian filter sigma=250m
names(dat2001)[ncol(dat2001)] <- "urbagNNS_Gauss250"

dat2001<-cbind(dat2001,extract(urbagNNS_Gauss750,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes urbag 250m resolution data, Gaussian filter sigma=750m
names(dat2001)[ncol(dat2001)] <- "urbagNNS_Gauss750"

dat2001<-cbind(dat2001,extract(watNNS,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes wat 250m resolution data
names(dat2001)[ncol(dat2001)] <- "watNNS"

dat2001<-cbind(dat2001,extract(watNNS_Gauss250,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes wat 250m resolution data, Gaussian filter sigma=250m
names(dat2001)[ncol(dat2001)] <- "watNNS_Gauss250"

dat2001<-cbind(dat2001,extract(watNNS_Gauss750,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes wat 250m resolution data, Gaussian filter sigma=750m
names(dat2001)[ncol(dat2001)] <- "watNNS_Gauss750"



### set up weight matrix for SS, and calculate weight values for each row in dat2001
samprast2001 <- rasterize(cbind(dat2001$X,dat2001$Y), r2, field=1)
sampsum25 <- focal(samprast2001, w=matrix(1/25, nc=5, nr=5), na.rm=TRUE)
dat2001 <- cbind(dat2001,extract(sampsum25,as.matrix(cbind(dat2001$X,dat2001$Y))))
names(dat2001)[ncol(dat2001)] <- "sampsum25"
dat2001$wt <- 1/dat2001$sampsum25

dat2001$SS <- as.character(dat2001$SS)
# dat2001$PCODE <- as.character(dat2001$PCODE)
# dat2001<-dat2001[,-c(25:43)] # remove columns with climate data (from avian dataset, since using Adaptwest Climate data instead)
write.csv(dat2001,"NNSdat2001.csv")

# Prepare point count data for each SS and aggregate for 2001 and 2011.

# PC <- inner_join(PCTBL,PKEY[,1:8],by=c("PKEY","SS"))[,-1]
# colnames(PC)[10]<-"PCODE"
# PC <- inner_join(PC,SS@data[,c(2,5)],by="SS")
# NSPC <- PC[PC$JURS=="NS",]

PC<-inner_join(PCcombo,PKEYcombo,by=c("PKEY"))
PC<-inner_join(PC,SScombo,by="SS")
names(PC)

NSPC <- PC[which(PC$SS%in%NSSS$SS),]

NSPC$SS <- as.character(NSPC$SS)
NSPC$PKEY <- as.character(NSPC$PKEY)
#NSPC$PCODE <- as.character(NSPC$PCODE)
NSPC$SPECIES <- as.character(NSPC$SPECIES)
NSPC2001 <- NSPC[NSPC$YEAR < 2006,] 
NSPC2011 <- NSPC[NSPC$YEAR > 2005,] 
survey2001 <- aggregate(NSPC2001$ABUND, by=list("PKEY"=NSPC2001$PKEY,"SS"=NSPC2001$SS,"ARU"=NSPC2001$ARU), FUN=sum) #n=9865
survey2011 <- aggregate(NSPC2011$ABUND, by=list("PKEY"=NSPC2011$PKEY,"SS"=NSPC2011$SS,"ARU"=NSPC2011$ARU), FUN=sum) #n=10505

speclist<-levels(as.factor(offcombo$SPECIES))



rm(list=setdiff(ls(),c("pred_abs_2011","LCC","speclist","offcombo","NSPC2001","survey2001","dat2001","NSPC2011","survey2011","dat2011")))
gc()
save.image("D:/CHID regional NS BRT/NS_BRT_Rproject/data_pack.RData")
