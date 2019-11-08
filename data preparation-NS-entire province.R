library(raster)
library(dismo)
library(rpart)
library(maptools)
library(data.table)
library(rgdal)
library(dplyr)

load("D:/Avian data processed/BAM_data_package_August2019.RData")

ltNS<-raster("D:/CHID regional NS BRT/spatial/landtypes_NovaScotia.tif")


NSshp<-readOGR(dsn="D:/basemaps",layer="bound_p")
NSshp<-subset(NSshp,STATEABB=="CA-NS")
plot(NSshp)
NSshp<-spTransform(NSshp,proj4string(ltNS))


NSras<-rasterize(NSshp,ltNS,field="COUNTRY")
plot(NSras)

NS2<-NSras
NS2<-projectRaster(NS2,crs=LCC,method="ngb")
plot(NS2)

# Load Beaudoin layers (2011 and 2001), crop and mask for study region, save as rasters

b2011 <- list.files("D:/Beaudoin/2011/",pattern="tif$")
setwd("D:/Beaudoin/2011/")
bs2011 <- stack(raster(b2011[1]))
for (i in 2:length(b2011)) { bs2011 <- addLayer(bs2011, raster(b2011[i]))}
names(bs2011) <- gsub("NFI_MODIS250m_2011_kNN_","",names(bs2011))
abs2011 <- crop(bs2011,NS2)
abs2011<-mask(abs2011,NS2)
writeRaster(abs2011, filename="D:/Beaudoin/2011/Processed/NovaScotia/abs2011_250m.grd", format="raster",overwrite=TRUE)
abs2011<-brick("D:/Beaudoin/2011/Processed/NovaScotia/abs2011_250m.grd")

b2001 <- list.files("D:/Beaudoin/2001/",pattern="tif$")
setwd("D:/Beaudoin/2001/")
bs2001 <- stack(raster(b2001[1]))
for (i in 2:length(b2001)) { bs2001 <- addLayer(bs2001, raster(b2001[i]))}
names(bs2001) <- gsub("NFI_MODIS250m_2001_kNN_","",names(bs2001))
abs2001 <- crop(bs2001,NS2)
abs2001<-mask(abs2001,NS2)
writeRaster(abs2001, filename="D:/Beaudoin/2001/Processed/NovaScotia/abs2001_250m.grd", format="raster",overwrite=TRUE)
abs2001<-brick("D:/Beaudoin/2001/Processed/NovaScotia/abs2001_250m.grd")

# obtain weighted sums of neighourhood cells using Gaussian filter with sigma=250, and 750m for Beaudoin and CTI layers, save outputs as rasters
#Beaudoin 2011
## sigma = 250m
fw250<-focalWeight(x=abs2011,d=250,type="Gauss")
abs2011_Gauss250<-stack(focal(abs2011[[1]],w=fw250,na.rm=TRUE))
names(abs2011_Gauss250)<-names(abs2011)[[1]]
for(i in 2:nlayers(abs2011)){
  abs2011_Gauss250<-addLayer(abs2011_Gauss250,focal(abs2011[[i]],w=fw250,na.rm=TRUE))
  names(abs2011_Gauss250)[i]<-names(abs2011)[[i]]
}
abs2011_Gauss250<-brick(abs2011_Gauss250)
writeRaster(abs2011_Gauss250, filename="D:/Beaudoin/2011/Processed/NovaScotia/abs2011_250_Gauss250m.grd", format="raster",overwrite=TRUE)
abs2011_Gauss250<-brick("D:/Beaudoin/2011/Processed/NovaScotia/abs2011_250_Gauss250m.grd")

## sigma = 750m
fw750<-focalWeight(x=abs2011,d=750,type="Gauss")
abs2011_Gauss750<-brick(focal(abs2011[[1]],w=fw750,na.rm=TRUE))
names(abs2011_Gauss750)<-names(abs2011)[[1]]
for(i in 2:nlayers(abs2011)){
  abs2011_Gauss750<-addLayer(abs2011_Gauss750,focal(abs2011[[i]],w=fw750,na.rm=TRUE))
  names(abs2011_Gauss750)[i]<-names(abs2011)[[i]]
}
abs2011_Gauss750<-brick(abs2011_Gauss750)
writeRaster(abs2011_Gauss750, filename="D:/Beaudoin/2011/Processed/NovaScotia/abs2011_250_Gauss750m.grd", format="raster",overwrite=TRUE)
abs2011_Gauss750<-brick("D:/Beaudoin/2011/Processed/NovaScotia/abs2011_250_Gauss750m.grd")



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
writeRaster(abs2001_Gauss250, filename="D:/Beaudoin/2001/Processed/NovaScotia/abs2001_250_Gauss250m.grd", format="raster",overwrite=TRUE)
#abs2001_Gauss250<-brick("D:/Beaudoin/2001/Processed/NovaScotia/abs2001_250_Gauss250m.grd")

## sigma = 750m
fw750<-focalWeight(x=abs2001,d=750,type="Gauss")
abs2001_Gauss750<-brick(focal(abs2001[[1]],w=fw750,na.rm=TRUE))
names(abs2001_Gauss750)<-names(abs2001)[[1]]
for(i in 2:nlayers(abs2001)){
  abs2001_Gauss750<-addLayer(abs2001_Gauss750,focal(abs2001[[i]],w=fw750,na.rm=TRUE))
  names(abs2001_Gauss750)[i]<-names(abs2001)[[i]]
}
abs2001_Gauss750<-brick(abs2001_Gauss750)
writeRaster(abs2001_Gauss750, filename="D:/Beaudoin/2001/Processed/NovaScotia/abs2001_250_Gauss750m.grd", format="raster",overwrite=TRUE)
#abs2001_Gauss750<-brick("D:/Beaudoin/2001/Processed/NovaScotia/abs2001_250_Gauss750m.grd")


# water table data. Water table data is a factor, so no sense in applying the Gaussian filter smoother. Instead, obtain the modal value over equal weight matrices of size 3x3 (250m GF analog) and 5x5 (750m GF analog) 
wtbl250 <- raster("D:/CHID regional NS BRT/spatial/watertableNovaScotia_LCC.tif")
wtbl250_modal3x3<-focal(wtbl250,w=matrix(1,3,3),na.rm=TRUE, fun=modal)
wtbl250_modal5x5<-focal(wtbl250,w=matrix(1,5,5),na.rm=TRUE, fun=modal)

writeRaster(wtbl250_modal3x3, filename="D:/CHID regional NS BRT/spatial/wtbl250_NS_modal3x3.asc", format="ascii",overwrite=TRUE)
writeRaster(wtbl250_modal5x5, filename="D:/CHID regional NS BRT/spatial/wtbl250_NS_modal5x5.asc", format="ascii",overwrite=TRUE)

# CTI data
CTI<-raster("D:/CTI/CTI_NS.tif")
CTI<-projectRaster(CTI,crs=LCC)
CTI<-crop(CTI,NS2)
CTI250<-resample(CTI,abs2011[[1]])
CTI250<-mask(CTI250,abs2011[[1]])

CTI250_Gauss250<-focal(CTI250,w=fw250,na.rm=TRUE)
CTI250_Gauss750<-focal(CTI250,w=fw750,na.rm=TRUE)

writeRaster(CTI250, filename="D:/CHID regional NS BRT/spatial/CTI250_NS.asc", format="ascii",overwrite=TRUE)
#CTI250<-raster("D:/CHID regional NS BRT/spatial/CTI250_NS.asc")
writeRaster(CTI250_Gauss250, filename="D:/CHID regional NS BRT/spatial/CTI250_NS_Gauss250.asc", format="ascii",overwrite=TRUE)
#CTI250_Gauss250<-raster("D:/CHID regional NS BRT/spatial/CTI250_NS_Gauss250.asc")
writeRaster(CTI250_Gauss750, filename="D:/CHID regional NS BRT/spatial/CTI250_NS_Gauss750.asc", format="ascii",overwrite=TRUE)
#CTI250_Gauss750<-raster("D:/CHID regional NS BRT/spatial/CTI250_NS_Gauss750.asc")

#roads (Venter et al)
roadsNS<-raster("D:/VenterEtAlFootprint/Maps/Roads.tif")
roadsNS<-projectRaster(roadsNS,NS2)
roadsNS<-crop(roadsNS,NS2)

roadsNS<-resample(roadsNS,abs2011[[1]])
roadsNS<-mask(roadsNS,abs2011[[1]])

roadsNS_Gauss250<-focal(roadsNS,w=fw250,na.rm=TRUE)
roadsNS_Gauss750<-focal(roadsNS,w=fw750,na.rm=TRUE)

# urban/agriculture
urbag <- raster("D:/urbag2011_lcc1/urbag2011_lcc1.tif")
urbagNS<- crop(urbag,NS2)
urbagNS<- mask(urbagNS,abs2011[[1]])

urbagNS_Gauss250<-focal(urbagNS,w=fw250,na.rm=TRUE)
urbagNS_Gauss750<-focal(urbagNS,w=fw750,na.rm=TRUE)

# water and "lake edge density"
wat <- raster("D:/wat2011_lcc1/wat2011_lcc1.tif")
watNS<- crop(wat,NS2)
watNS<- mask(watNS,abs2011[[1]])

watNS_Gauss250<-focal(watNS,w=fw250,na.rm=TRUE)
watNS_Gauss750<-focal(watNS,w=fw750,na.rm=TRUE)

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

names(roadsNS) <- "roadsNS"
names(roadsNS_Gauss250) <- "roadsNS_Gauss250"
names(roadsNS_Gauss750) <- "roadsNS_Gauss750"

names(urbagNS) <- "urbagNS"
names(urbagNS_Gauss250) <- "urbagNS_Gauss250"
names(urbagNS_Gauss750) <- "urbagNS_Gauss750"

names(watNS) <- "watNS"
names(watNS_Gauss250) <- "watNS_Gauss250"
names(watNS_Gauss750) <- "watNS_Gauss750"

ARU<-raster(extent(abs2011),crs=LCC,resolution=res(abs2011),vals=0)
names(ARU)<-"ARU"

pred_abs_2011<-stack(abs2011,abs2011_Gauss250,abs2011_Gauss750,wtbl250,wtbl250_modal3x3,wtbl250_modal5x5,CTI250,CTI250_Gauss250,CTI250_Gauss750,roadsNS,roadsNS_Gauss250,roadsNS_Gauss750,urbagNS, urbagNS_Gauss250, urbagNS_Gauss750, watNS, watNS_Gauss250, watNS_Gauss750,ARU, quick=TRUE)
names(pred_abs_2011)

writeRaster(pred_abs_2011, filename="D:/CHID regional NS BRT/prediction dataset/abs2011_250m_NS.grd", format="raster",overwrite=TRUE)
pred_abs_2011<-stack("D:/CHID regional NS BRT/prediction dataset/abs2011_250m_NS.grd")

# Extracting data from rasters for surveyed locations (ABSS)
## 2011

spSS<-SpatialPointsDataFrame(coords=SScombo[,2:3],data=SScombo,proj4string = LCC)
NSSS <- crop(spSS,ltNS)
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

dat2011<-cbind(dat2011,extract(roadsNS,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes roads 250m resolution data
names(dat2011)[ncol(dat2011)] <- "roadsNS"

dat2011<-cbind(dat2011,extract(roadsNS_Gauss250,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes roads 250m resolution data, Gaussian filter sigma=250m
names(dat2011)[ncol(dat2011)] <- "roadsNS_Gauss250"

dat2011<-cbind(dat2011,extract(roadsNS_Gauss750,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes roads 250m resolution data, Gaussian filter sigma=750m
names(dat2011)[ncol(dat2011)] <- "roadsNS_Gauss750"

dat2011<-cbind(dat2011,extract(urbagNS,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes urbag 250m resolution data
names(dat2011)[ncol(dat2011)] <- "urbagNS"

dat2011<-cbind(dat2011,extract(urbagNS_Gauss250,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes urbag 250m resolution data, Gaussian filter sigma=250m
names(dat2011)[ncol(dat2011)] <- "urbagNS_Gauss250"

dat2011<-cbind(dat2011,extract(urbagNS_Gauss750,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes urbag 250m resolution data, Gaussian filter sigma=750m
names(dat2011)[ncol(dat2011)] <- "urbagNS_Gauss750"

dat2011<-cbind(dat2011,extract(watNS,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes wat 250m resolution data
names(dat2011)[ncol(dat2011)] <- "watNS"

dat2011<-cbind(dat2011,extract(watNS_Gauss250,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes wat 250m resolution data, Gaussian filter sigma=250m
names(dat2011)[ncol(dat2011)] <- "watNS_Gauss250"

dat2011<-cbind(dat2011,extract(watNS_Gauss750,as.matrix(cbind(dat2011$X,dat2011$Y)))) # includes wat 250m resolution data, Gaussian filter sigma=750m
names(dat2011)[ncol(dat2011)] <- "watNS_Gauss750"

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
write.csv(dat2011,"NS_dat2011.csv")


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

dat2001<-cbind(dat2001,extract(roadsNS,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes roads 250m resolution data
names(dat2001)[ncol(dat2001)] <- "roadsNS"

dat2001<-cbind(dat2001,extract(roadsNS_Gauss250,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes roads 250m resolution data, Gaussian filter sigma=250m
names(dat2001)[ncol(dat2001)] <- "roadsNS_Gauss250"

dat2001<-cbind(dat2001,extract(roadsNS_Gauss750,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes roads 250m resolution data, Gaussian filter sigma=750m
names(dat2001)[ncol(dat2001)] <- "roadsNS_Gauss750"

dat2001<-cbind(dat2001,extract(urbagNS,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes urbag 250m resolution data
names(dat2001)[ncol(dat2001)] <- "urbagNS"

dat2001<-cbind(dat2001,extract(urbagNS_Gauss250,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes urbag 250m resolution data, Gaussian filter sigma=250m
names(dat2001)[ncol(dat2001)] <- "urbagNS_Gauss250"

dat2001<-cbind(dat2001,extract(urbagNS_Gauss750,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes urbag 250m resolution data, Gaussian filter sigma=750m
names(dat2001)[ncol(dat2001)] <- "urbagNS_Gauss750"

dat2001<-cbind(dat2001,extract(watNS,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes wat 250m resolution data
names(dat2001)[ncol(dat2001)] <- "watNS"

dat2001<-cbind(dat2001,extract(watNS_Gauss250,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes wat 250m resolution data, Gaussian filter sigma=250m
names(dat2001)[ncol(dat2001)] <- "watNS_Gauss250"

dat2001<-cbind(dat2001,extract(watNS_Gauss750,as.matrix(cbind(dat2001$X,dat2001$Y)))) # includes wat 250m resolution data, Gaussian filter sigma=750m
names(dat2001)[ncol(dat2001)] <- "watNS_Gauss750"


### set up weight matrix for SS, and calculate weight values for each row in dat2001
samprast2001 <- rasterize(cbind(dat2001$X,dat2001$Y), r2, field=1)
sampsum25 <- focal(samprast2001, w=matrix(1/25, nc=5, nr=5), na.rm=TRUE)
dat2001 <- cbind(dat2001,extract(sampsum25,as.matrix(cbind(dat2001$X,dat2001$Y))))
names(dat2001)[ncol(dat2001)] <- "sampsum25"
dat2001$wt <- 1/dat2001$sampsum25

dat2001$SS <- as.character(dat2001$SS)
# dat2001$PCODE <- as.character(dat2001$PCODE)
# dat2001<-dat2001[,-c(25:43)] # remove columns with climate data (from avian dataset, since using Adaptwest Climate data instead)
write.csv(dat2001,"NS_dat2001.csv")

# Prepare point count data for each SS and aggregate for 2001 and 2011.
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

# create layer that identifies management unit
NNS <- raster("D:/CHID regional NS BRT/spatial/NNS.tif") 
plot(NNS)
NNS<-projectRaster(NNS,pred_abs_2011,method="ngb")
values(NNS)[!is.na(values(NNS))]<-1


managt.unit<-pred_abs_2011[[1]]
values(managt.unit)[!is.na(values(managt.unit))]<-2
plot(managt.unit)

values(managt.unit)[base::intersect(which(values(NNS)==1),which(values(managt.unit)==2))]<-1
summary(values(managt.unit))
names(managt.unit)<-" management_unit"

pred_abs_2011<- addLayer(pred_abs_2011, managt.unit)

rm(list=setdiff(ls(),c("pred_abs_2011","LCC","speclist","offcombo","NSPC2001","survey2001","dat2001","NSPC2011","survey2011","dat2011")))
gc()
save.image("D:/CHID regional NS BRT/NS_BRT_Rproject/data_pack_NS.RData")
