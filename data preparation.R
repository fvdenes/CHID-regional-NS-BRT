library(raster)
library(dismo)
library(rpart)
library(maptools)
library(data.table)
library(rgdal)
library(dplyr)


load("D:/Avian data processed/data_package_2016-04-18.Rdata")	
load("D:/Avian data processed/offsets-v3_2016-04-18.Rdata")
LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
coordinates(SS) <- c("X", "Y") 
proj4string(SS) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
SSLCC <- spTransform(SS, LCC)

offl <- data.table(melt(OFF))
names(offl) <- c("PKEY","SPECIES","logoffset")
offl$SPECIES <- as.character(offl$SPECIES)
offl$PKEY <-as.character(offl$PKEY)
rm(OFF)

#eco<-shapefile("D:/CHID regional NS BRT/NNS_eco.shp")
#eco<-rasterize(eco,NNS)
#writeRaster(eco, filename="NNS_eco", format="GTiff",overwrite=TRUE)
eco <- raster("D:/CHID regional NS BRT/NNS_eco.tif") 
nalc <- raster("D:/CHID regional Alberta BRT/NA_LandCover_2005_LCC.img")
NNS <- raster("D:/CHID regional NS BRT/NNS.tif") 
proj4string(NNS)<-LCC


NNS_shp<-readOGR("D:/CHID regional NS BRT",layer="northernNSmngtunit")
NNS_shp<-spTransform(NNS_shp,LCC)
NSSS <- crop(SSLCC,NNS_shp)
NSSS <- as.data.frame(NSSS)

b2011 <- list.files("D:/Beaudoin/2011/",pattern="tif$")
setwd("D:/Beaudoin/2011/")
bs2011 <- stack(raster(b2011[1]))
for (i in 2:length(b2011)) { bs2011 <- addLayer(bs2011, raster(b2011[i]))}
names(bs2011) <- gsub("NFI_MODIS250m_2011_kNN_","",names(bs2011))

abs2011 <- crop(bs2011,NNS)
abs2011<-mask(abs2011,NNS)
abs2011_1km <- aggregate(abs2011, fact=4, fun=mean)
r2 <- abs2011_1km[[1]]

ecor1km <- resample(eco, abs2011_1km)
abs2011_1km <- addLayer(abs2011_1km, ecor1km)
names(abs2011_1km)[nlayers(abs2011_1km)] <- "eco"


wtbl <- raster("D:/CHID regional NS BRT/watertableNNS_LCC.tif")
wtbl1km<-aggregate(wtbl,fact=4,fun=modal)
abs2011_1km <- addLayer(abs2011_1km, wtbl1km)
names(abs2011_1km)[nlayers(abs2011_1km)] <- "wtbl"

writeRaster(abs2011_1km,"NS2011rasters",overwrite=TRUE)

roadsNNS<-raster("D:/CHID regional NS BRT/roadsNNS.tif")
roadsNNS<-crop(roadsNNS,abs2011_1km)
roadsNNS<-resample(roadsNNS,abs2011_1km)
roadsNNS<-mask(roadsNNS,abs2011_1km$LandCover_NonVeg_v1)
abs2011_1km <- addLayer(abs2011_1km, roadsNNS)
names(abs2011_1km)[nlayers(abs2011_1km)] <- "roads"

climateAW2010 <- list.files("D:/ClimateAdaptWest/baseline19812010/",pattern="asc$")
setwd("D:/ClimateAdaptWest/baseline19812010/")
clim2010 <- stack(raster(climateAW2010[1]))
for (i in 2:length(climateAW2010)) { clim2010 <- addLayer(clim2010, raster(climateAW2010[i]))}
proj4string(clim2010)<-LCC
aclim2010 <- crop(clim2010,abs2011_1km)
aclim2010<-resample(clim2010,abs2011_1km)
aclim2010<-mask(aclim2010,abs2011_1km$LandCover_NonVeg_v1)

for(i in 1:length(names(aclim2010))){ 
  abs2011_1km <- addLayer(abs2011_1km, aclim2010[[i]])
  names(abs2011_1km)[nlayers(abs2011_1km)] <- names(aclim2010[[i]])
}



b2001 <- list.files("D:/Beaudoin/2001/",pattern="tif$")
setwd("D:/Beaudoin/2001/")
bs2001 <- stack(raster(b2001[1]))
for (i in 2:length(b2001)) { bs2001 <- addLayer(bs2001, raster(b2001[i]))}
names(bs2001) <- gsub("NFI_MODIS250m_2001_kNN_","",names(bs2001))
abs2001 <- crop(bs2001,NNS)


dat2011 <- cbind(NSSS, extract(abs2011,as.matrix(cbind(NSSS$X,NSSS$Y))))
dat2011 <-cbind(dat2011,extract(nalc,as.matrix(cbind(dat2011$X,dat2011$Y)))) 
names(dat2011)[ncol(dat2011)] <- "LCC"
dat2011 <-cbind(dat2011,extract(eco,as.matrix(cbind(dat2011$X,dat2011$Y))))
names(dat2011)[ncol(dat2011)] <- "eco"

disaggregated_wtbl<-raster("D:/CHID regional NS BRT/disaggregated_25m_watertableNNS_LCC.tif")
dat2011<-cbind(dat2011,extract(disaggregated_wtbl,as.matrix(cbind(dat2011$X,dat2011$Y))))
names(dat2011)[ncol(dat2011)] <- "wtbl"

dat2011<-cbind(dat2011,extract(roadsNNS,as.matrix(cbind(dat2011$X,dat2011$Y))))
names(dat2011)[ncol(dat2011)] <- "roads"

for(i in which(names(abs2011_1km)=="AHM"):which(names(abs2011_1km)=="TD") ) {
  dat2011<-cbind(dat2011,extract(abs2011_1km[[i]],as.matrix(cbind(NSSS$X,NSSS$Y))))
  names(dat2011)[ncol(dat2011)] <- names(abs2011_1km)[[i]]
}


samprast2011 <- rasterize(cbind(dat2011$X,dat2011$Y), r2, field=1)
sampsum25 <- focal(samprast2011, w=matrix(1/25, nc=5, nr=5), na.rm=TRUE)

dat2011 <- cbind(dat2011,extract(sampsum25,as.matrix(cbind(dat2011$X,dat2011$Y))))
names(dat2011)[ncol(dat2011)] <- "sampsum25"
dat2011$wt <- 1/dat2011$sampsum25


dat2011$SS <- as.character(dat2011$SS)
dat2011$PCODE <- as.character(dat2011$PCODE)

setwd("D://CHID regional NS BRT/")
write.csv(dat2011,"NNSdat2011.csv")


dat2001 <- cbind(NSSS, extract(abs2001,as.matrix(cbind(NSSS$X,NSSS$Y))))
dat2001 <-cbind(dat2001,extract(nalc,as.matrix(cbind(dat2001$X,dat2001$Y)))) 
names(dat2001)[ncol(dat2001)] <- "LCC"
dat2001 <-cbind(dat2001,extract(eco,as.matrix(cbind(dat2001$X,dat2001$Y))))
names(dat2001)[ncol(dat2001)] <- "eco"



dat2001<-cbind(dat2001,extract(disaggregated_wtbl,as.matrix(cbind(dat2001$X,dat2001$Y))))
names(dat2001)[ncol(dat2001)] <- "wtbl"
dat2001<-cbind(dat2001,extract(roadsNNS,as.matrix(cbind(dat2001$X,dat2001$Y))))
names(dat2001)[ncol(dat2001)] <- "roads"

for(i in which(names(abs2011_1km)=="AHM"):which(names(abs2011_1km)=="TD") ) {
  dat2001<-cbind(dat2001,extract(abs2011_1km[[i]],as.matrix(cbind(NSSS$X,NSSS$Y))))
  names(dat2001)[ncol(dat2001)] <- names(abs2011_1km)[[i]]
}

samprast2001 <- rasterize(cbind(dat2001$X,dat2001$Y), r2, field=1)
sampsum25 <- focal(samprast2001, w=matrix(1/25, nc=5, nr=5), na.rm=TRUE) 
dat2001 <- cbind(dat2001,extract(sampsum25,as.matrix(cbind(dat2001$X,dat2001$Y))))
names(dat2001)[ncol(dat2001)] <- "sampsum25"
dat2001$wt <- 1/dat2001$sampsum25

dat2001<-dat2001[-which(is.na(dat2001$wt)),] 

dat2001$SS <- as.character(dat2001$SS)
dat2001$PCODE <- as.character(dat2001$PCODE)
write.csv(dat2001,"NNSdat2001.csv")

PC <- inner_join(PCTBL,PKEY[,1:8],by=c("PKEY","SS"))[,-1]
colnames(PC)[10]<-"PCODE"
PC <- inner_join(PC,SS@data[,c(2,5)],by="SS")


NSPC <- PC[which(PC$SS%in%NSSS$SS),]

NSPC$SS <- as.character(NSPC$SS)
NSPC$PKEY <- as.character(NSPC$PKEY)
NSPC$PCODE <- as.character(NSPC$PCODE)
NSPC$SPECIES <- as.character(NSPC$SPECIES)
NSPC2001 <- NSPC[NSPC$YEAR < 2006,] #n=83570
NSPC2011 <- NSPC[NSPC$YEAR > 2005,] #n=91656
survey2001 <- aggregate(NSPC2001$ABUND, by=list("PKEY"=NSPC2001$PKEY,"SS"=NSPC2001$SS,"PCODE"=NSPC2001$PCODE), FUN=sum) #n=9865
survey2011 <- aggregate(NSPC2011$ABUND, by=list("PKEY"=NSPC2011$PKEY,"SS"=NSPC2011$SS,"PCODE"=NSPC2011$PCODE), FUN=sum) #n=10505

speclist<-levels(as.factor(offl$SPECIES))

save.image("D:/CHID regional NS BRT/data_pack.RData")

