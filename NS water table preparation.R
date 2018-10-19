library(raster)
library(dismo)
library(rpart)
library(maptools)
library(data.table)
library(rgdal)
library(dplyr)
library(sf)
library(fasterize)

LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
NNS <- raster("D:/CHID regional NS BRT/NNS.tif") 

wt1<-readOGR("D:/CHID regional NS BRT/water_table/wtbl_01/wtbl.shp")
wt1<-spTransform(wt1,LCC)
wt1_sf<-st_as_sf(wt1)
wt1r<-fasterize(wt1_sf,NNS,field="Depth")
rm(wt1)
plot(wt1r)

wt2<-readOGR("D:/CHID regional NS BRT/water_table/wtbl_02/wtbl.shp")
wt2<-spTransform(wt2,LCC)
wt2_sf<-st_as_sf(wt2)
wt2r<-fasterize(wt2_sf,NNS,field="Depth")
rm(wt2)
plot(wt2r)

wt3<-readOGR("D:/CHID regional NS BRT/water_table/wtbl_03/wtbl.shp")
wt3<-spTransform(wt3,LCC)
wt3_sf<-st_as_sf(wt3)
wt3r<-fasterize(wt3_sf,NNS,field="Depth")
rm(wt3)

wt4<-readOGR("D:/CHID regional NS BRT/water_table/wtbl_04/wtbl.shp")
wt4<-spTransform(wt4,LCC)
wt4_sf<-st_as_sf(wt4)
wt4r<-fasterize(wt4_sf,NNS,field="Depth")
rm(wt4)

wt5<-readOGR("D:/CHID regional NS BRT/water_table/wtbl_05/wtbl.shp")
wt5<-spTransform(wt5,LCC)
wt5_sf<-st_as_sf(wt5)
wt5r<-fasterize(wt5_sf,NNS,field="Depth")
rm(wt5)

wt6<-readOGR("D:/CHID regional NS BRT/water_table/wtbl_06/wtbl.shp")
wt6<-spTransform(wt6,LCC)
wt6_sf<-st_as_sf(wt6)
wt6r<-fasterize(wt6_sf,NNS,field="Depth")
rm(wt6)

wt7<-readOGR("D:/CHID regional NS BRT/water_table/wtbl_07/wtbl.shp")
wt7<-spTransform(wt7,LCC)
wt7_sf<-st_as_sf(wt7)
wt7r<-fasterize(wt7_sf,NNS,field="Depth")
rm(wt7)

wt8<-readOGR("D:/CHID regional NS BRT/water_table/wtbl_08/wtbl.shp")
wt8<-spTransform(wt8,LCC)
wt8_sf<-st_as_sf(wt8)
wt8r<-fasterize(wt8_sf,NNS,field="Depth")
rm(wt8)

wt9<-readOGR("D:/CHID regional NS BRT/water_table/wtbl_09/wtbl.shp")
wt9<-spTransform(wt9,LCC)
wt9_sf<-st_as_sf(wt9)
wt9r<-fasterize(wt9_sf,NNS,field="Depth")
rm(wt9)

wt10<-readOGR("D:/CHID regional NS BRT/water_table/wtbl_10/wtbl.shp")
wt10<-spTransform(wt10,LCC)
wt10_sf<-st_as_sf(wt10)
wt10r<-fasterize(wt10_sf,NNS,field="Depth")
rm(wt10)

wt11<-readOGR("D:/CHID regional NS BRT/water_table/wtbl_11/wtbl.shp")
wt11<-spTransform(wt11,LCC)
wt11_sf<-st_as_sf(wt11)
wt11r<-fasterize(wt11_sf,NNS,field="Depth")
rm(wt11)

wt12<-readOGR("D:/CHID regional NS BRT/water_table/wtbl_12/wtbl.shp")
wt12<-spTransform(wt12,LCC)
wt12_sf<-st_as_sf(wt12)
wt12r<-fasterize(wt12_sf,NNS,field="Depth")
rm(wt12)

wt13<-readOGR("D:/CHID regional NS BRT/water_table/wtbl_13/wtbl.shp")
wt13<-spTransform(wt13,LCC)
wt13_sf<-st_as_sf(wt13)
wt13r<-fasterize(wt13_sf,NNS,field="Depth")
rm(wt13)

wt14<-readOGR("D:/CHID regional NS BRT/water_table/wtbl_14/wtbl.shp")
wt14<-spTransform(wt14,LCC)
wt14_sf<-st_as_sf(wt14)
wt14r<-fasterize(wt14_sf,NNS,field="Depth")
rm(wt14)

wtmerge<-mosaic(wt1r,wt2r,wt3r,wt4r,wt5r,wt6r,wt7r,wt8r,wt9r,wt10r,wt11r,wt12r,wt13r,wt14r, fun=mean)
plot(wtmerge)

summary(getValues(wtmerge))

cc<-which(getValues(wtmerge)==1.5)
wtmerge@data@values[cc]<-1

nn<-which(getValues(wtmerge)==2.5)
wtmerge@data@values[nn]<-2

mm<-which(getValues(wtmerge)==3.5)
wtmerge@data@values[mm]<-3

plot(wtmerge)

wtmerge_mask<- mask(wtmerge,NNS)

plot(wtmerge_mask)

writeRaster(wtmerge_mask,"watertableNNS_LCC",format="GTiff",overwrite=TRUE)
