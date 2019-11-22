library(raster)
library(rpart)
library(maptools)
library(data.table)
library(rgdal)
library(dplyr)
library(blockCV)
library(sp)
library(dismo)
library(ggplot2)

# A function that fits the BRT model ('gbm.step' from dismo package) on pre-defined folds, and saves outputs ####
brt_blocks <- function(data = datcombo, pred.stack = pred_abs_2011, seed = 1222, pred.variables ,output.folder, blocks=NULL, keep.out = TRUE, tc=3,lr=0.001,bf=0.5, save.points.shp=FALSE){ 
  # Arguments for this function
  ## data: data.frame object containing data for model fitting
  ## pred.stack: the raster stack/brick used as prediction dataset
  ## pred.variables: a character vector giving the names of predictor variables that will be included in BRT models
  ## blocks: object resulting from 'spatialBlocks' function that contains classification of sample points into folds
  ## output.folder: path to output folder (if folder does not exist it is created)
  ## keep.out: logical, whether to keep the output in the workspace. both blocks object and the brt output are automatically saved in the output folder regardless
  ## tc: BRT tree complexity
  ## lr: BRT learning rate
  ## bf: BRT bag fraction
  ## save.points.shp: logical, whether to save survey points as a shapefile in output folder
  
  # fit BRT models using pre-determined folds for CV
  if (is.null(blocks)){
    folds<-NULL
    n.folds<-10
  }
  
  else {
    folds<-blocks$foldID
    n.folds<-blocks$k
  }
  
  x1 <-
    try(brt <-
          gbm.step(
            datcombo,
            gbm.y = "ABUND",
            gbm.x = pred.variables,
            fold.vector = folds,
            n.folds = n.folds,
            family = "poisson",
            tree.complexity = tc,
            learning.rate = lr,
            bag.fraction = bf,
            offset = datcombo$logoffset,
            site.weights = datcombo$wt,
            keep.fold.models = T,
            keep.fold.fit = T
          ))
  
  if(class(x1)=="NULL"){#retry models that didn't converge with smaller learning rate
    x1 <-
      try(brt <-
            gbm.step(
              datcombo,
              gbm.y = "ABUND",
              gbm.x = pred.variables,
              fold.vector = folds,
              n.folds = n.folds,
              family = "poisson",
              tree.complexity = tc,
              learning.rate = lr/10,
              bag.fraction = bf,
              offset = datcombo$logoffset,
              site.weights = datcombo$wt,
              keep.fold.models = T,
              keep.fold.fit = T
            ))
  }
  
  if(class(x1)=="NULL"){#retry models that didn't converge with smaller learning rate
    x1 <-
      try(brt <-
            gbm.step(
              datcombo,
              gbm.y = "ABUND",
              gbm.x = pred.variables,
              fold.vector = folds,
              n.folds = n.folds,
              family = "poisson",
              tree.complexity = tc,
              learning.rate = lr/100,
              bag.fraction = bf,
              offset = datcombo$logoffset,
              site.weights = datcombo$wt,
              keep.fold.models = T,
              keep.fold.fit = T
            ))
  }
  
  if(class(x1)=="NULL"){
    stop("Restart model with even smaller learning rate, or other predictors!")
  }
  
  # Define/create folders for storing outputs
  if (class(x1) != "try-error") {
    z <- output.folder
    
    if (file.exists(z) == FALSE) {
      dir.create(z)
    }
    
    if (is.null(blocks)){
      save(brt, file = paste(z, speclist[j], "brtAB.R", sep = ""))
    }
    
    else {
      save(blocks, file = paste(z, speclist[j], "blocks.R", sep = ""))
      save(brt, file = paste(z, speclist[j], "brtAB.R", sep = ""))
    }  
    
    
    ## Model evaluation
    varimp <- as.data.frame(brt$contributions)
    write.csv(varimp, file = paste(z, speclist[j], "varimp.csv", sep = ""))
    cvstats <- t(as.data.frame(brt$cv.statistics))
    write.csv(cvstats, file = paste(z, speclist[j], "cvstats.csv", sep =
                                      ""))
    pdf(paste(z, speclist[j], "_plot.pdf", sep = ""))
    gbm.plot(
      brt,
      n.plots = length(pred.variables),
      smooth = TRUE,
      plot.layout = c(3, 3),
      common.scale = T
    )
    dev.off()
    pdf(paste(z, speclist[j], "_plot.var-scale.pdf", sep = ""))
    gbm.plot(
      brt,
      n.plots = length(pred.variables),
      smooth = TRUE,
      plot.layout = c(3, 3),
      common.scale = F,
      write.title = F
    )
    dev.off()
    
    
    ## Model prediction
    
    rast <-
      predict(pred.stack,
              brt,
              type = "response",
              n.trees = brt$n.trees)
    writeRaster(
      rast,
      filename = paste(z, speclist[j], "_pred1km", sep = ""),
      format = "GTiff",
      overwrite = TRUE
    )
    
    data_sp <-SpatialPointsDataFrame(coords = data[,c("X","Y")], data = data, proj4string = LCC)
    png(paste(z, speclist[j], "_pred1km.png", sep = ""))
    plot(rast, zlim = c(0, 1))
    points(data_sp$X, data_sp$Y, cex = 0.05)
    dev.off()
    
    if(save.points.shp==T){
      writeOGR(
        data_sp,
        dsn = paste(z, "surveypoints.shp", sep = ""),
        layer = "data_sp",
        driver = "ESRI Shapefile"
      )
    }
    
    if(keep.out==T) {return(brt)}
  }
}

#load data and prepare objects ####
load("D:/CHID regional NS BRT/NS_BRT_Rproject/data_pack.RData")

#pred_abs_2011<- brick("D:/CHID regional NS BRT/prediction dataset/abs2011_250m.grd")

j<-which(speclist=="CAWA") 

specoff <- filter(offcombo, SPECIES==as.character(speclist[j]))
specoff <- distinct(specoff)

specdat2001 <- filter(NSPC2001, SPECIES == as.character(speclist[j]))
specdat2001x <- aggregate(specdat2001$ABUND,by=list("PKEY"=specdat2001$PKEY,"SS"=specdat2001$SS), FUN=sum)
names(specdat2001x)[3] <- "ABUND"
dat1 <- right_join(specdat2001x,survey2001[,1:3],by=c("SS","PKEY")) #n=31864
dat1$SPECIES <- as.character(speclist[j])
dat1$ABUND <- as.integer(ifelse(is.na(dat1$ABUND),0,dat1$ABUND)) 
s2001 <- left_join(dat1,specoff, by=c("SPECIES","PKEY"))
d2001 <- left_join(s2001, dat2001, by=c("SS")) 

specdat2011 <- filter(NSPC2011, SPECIES == as.character(speclist[j]))
specdat2011x <- aggregate(specdat2011$ABUND,by=list("PKEY"=specdat2011$PKEY,"SS"=specdat2011$SS), FUN=sum)
names(specdat2011x)[3] <- "ABUND"  
dat2 <- right_join(specdat2011x,survey2011[,1:3],by=c("SS","PKEY"))
dat2$SPECIES <- as.character(speclist[j])
dat2$ABUND <- as.integer(ifelse(is.na(dat2$ABUND),0,dat2$ABUND)) 
s2011 <- left_join(dat2,specoff, by=c("SPECIES","PKEY"))
d2011 <- left_join(s2011, dat2011, by=c("SS")) 


datcombo <- rbind(d2001,d2011)


#convert water table covariate to factor (the 5 levels are 1: 0-0.10m; 2: 0.11-0.50m; 30.51-2m; 4: 2.01-10m; 5: >10m)
datcombo$wtbl250 <- as.factor(datcombo$wtbl250)
datcombo$wtbl250_modal3x3 <- as.factor(datcombo$wtbl250_modal3x3)
datcombo$wtbl250_modal5x5 <- as.factor(datcombo$wtbl250_modal5x5)

#convert water covariate at cell level to factor. (0: no water; 1: water)
datcombo$watNNS <- as.factor(datcombo$watNNS)

rm(list=setdiff(ls(),c("datcombo","pred_abs_2011","w","LCC","speclist","randomCV_brt","j","brt_blocks")))
gc()

# create a column converting abundance to occupancy
datcombo$OCCU <- 0 
datcombo$OCCU[which(datcombo$ABUND > 0)] <- 1

# convert data into SpatialPointDataFrame class
datcombo_sp <-SpatialPointsDataFrame(coords = datcombo[,c("X","Y")], data = datcombo, proj4string = LCC)


# Full model ####
# list variables for full model
pred.variables<-c( 
  "Species_Abie_Bal_v1",#
  "Species_Abie_Bal_v1_Gauss250",#
  "Species_Abie_Bal_v1_Gauss750",#
  "Species_Acer_Sah_v1",#
  "Species_Acer_Sah_v1_Gauss250",#
  "Species_Acer_Sah_v1_Gauss750",#
  "Species_Acer_Rub_v1",#
  "Species_Acer_Rub_v1_Gauss250",#
  "Species_Acer_Rub_v1_Gauss750",#
  "Species_Acer_Spi_v1",#
  "Species_Acer_Spi_v1_Gauss250",#
  "Species_Acer_Spi_v1_Gauss750",#
  "Species_Betu_All_v1",#
  "Species_Betu_All_v1_Gauss250",#
  "Species_Betu_All_v1_Gauss750",#
  "Species_Betu_Pap_v1",#
  "Species_Betu_Pap_v1_Gauss250",#
  "Species_Betu_Pap_v1_Gauss750",#
  "Species_Fagu_Gra_v1",#
  "Species_Fagu_Gra_v1_Gauss250",#
  "Species_Fagu_Gra_v1_Gauss750",#
  "Species_Frax_Ame_v1",#
  "Species_Frax_Ame_v1_Gauss250",#
  "Species_Frax_Ame_v1_Gauss750",#
  "Species_Lari_Lar_v1",#
  "Species_Lari_Lar_v1_Gauss250",#
  "Species_Lari_Lar_v1_Gauss750",#
  "Species_Pice_Gla_v1",#
  "Species_Pice_Gla_v1_Gauss250",#
  "Species_Pice_Gla_v1_Gauss750",#
  "Species_Pice_Rub_v1",#
  "Species_Pice_Rub_v1_Gauss250",#
  "Species_Pice_Rub_v1_Gauss750",#
  "Species_Pice_Mar_v1",#
  "Species_Pice_Mar_v1_Gauss250",#
  "Species_Pice_Mar_v1_Gauss750",#
  "Species_Pinu_Ban_v1",#
  "Species_Pinu_Ban_v1_Gauss250",#
  "Species_Pinu_Ban_v1_Gauss750",#
  "Species_Pinu_Res_v1",#
  "Species_Pinu_Res_v1_Gauss250",#
  "Species_Pinu_Res_v1_Gauss750",#
  "Species_Pinu_Str_v1",#
  "Species_Pinu_Str_v1_Gauss250",#
  "Species_Pinu_Str_v1_Gauss750",#
  "Species_Popu_Bal_v1",#
  "Species_Popu_Bal_v1_Gauss250",#
  "Species_Popu_Bal_v1_Gauss750",#
   "Species_Quer_Rub_v1",#
  "Species_Quer_Rub_v1_Gauss250",#
  "Species_Quer_Rub_v1_Gauss750",#
   "Species_Tsug_Can_v1",#
  "Species_Tsug_Can_v1_Gauss250",#
  "Species_Tsug_Can_v1_Gauss750",#
  "Structure_Biomass_TotalLiveAboveGround_v1",
  "Structure_Biomass_TotalLiveAboveGround_v1_Gauss250",
  "Structure_Biomass_TotalLiveAboveGround_v1_Gauss750",
  "Structure_Stand_Age_v1",
  "Structure_Stand_Age_v1_Gauss250",
  "Structure_Stand_Age_v1_Gauss750",
  "wtbl250",
  "wtbl250_modal3x3",
  "wtbl250_modal5x5",
  "CTI250",
  "CTI250_Gauss250",
  "CTI250_Gauss750",
  "roadsNNS",
  "roadsNNS_Gauss250",
  "roadsNNS_Gauss750",
  "urbagNNS",
  "urbagNNS_Gauss250",
  "urbagNNS_Gauss750",
  "watNNS",
  "watNNS_Gauss250",
  "watNNS_Gauss750",
  "ARU"
)
  
  # create object storing the indices of predictor layers (from the prediction dataset) for the autocorrelation assessment that will inform block size. Only include continuous numeric variables here. Can use indexing "[-c(x:y)] to exclude them (e.g. exclude climate variables, HF, etc).
  ind <- which(names(pred_abs_2011) %in% pred.variables[-c(61:63,73)])
  
  # Calculate autocorrelation range
  start_time <- Sys.time()
  sp.auto.arr1 <- spatialAutoRange(pred_abs_2011[[ind]])
  end_time <- Sys.time()
  end_time - start_time

# Use spatial blocks to separate train and test folds
start_time <- Sys.time()
set.seed(123)
sp_block_full <-  spatialBlock(
  speciesData = datcombo_sp,
  species = "OCCU",
  rasterLayer = pred_abs_2011[[1]],
  iteration = 250,
  theRange = sp.auto.arr1$range,
  selection = "random",
  maskBySpecies = FALSE,
  k=5
)
end_time <- Sys.time()
end_time - start_time

start_time<-Sys.time()
brt1<- brt_blocks(data=datcombo,pred.variables = pred.variables, lr=0.01,tc=3, output.folder = "D://CHID regional NS BRT/BRT_outputs/full_model_new/", blocks=sp_block_full, save.points.shp = TRUE)  
end_time<-Sys.time()
end_time-start_time





# Simplified model (selecting most influential scales) ####
pred.variables2<-c( 
  #"Species_Abie_Bal_v1",#
  #"Species_Abie_Bal_v1_Gauss250",#
  "Species_Abie_Bal_v1_Gauss750",#
  "Species_Acer_Sah_v1",#
  #"Species_Acer_Sah_v1_Gauss250",#
  #"Species_Acer_Sah_v1_Gauss750",#
  "Species_Acer_Rub_v1",#
  #"Species_Acer_Rub_v1_Gauss250",#
  #"Species_Acer_Rub_v1_Gauss750",#
  #"Species_Acer_Spi_v1",#
  #"Species_Acer_Spi_v1_Gauss250",#
  "Species_Acer_Spi_v1_Gauss750",#
  "Species_Betu_All_v1",#
  #"Species_Betu_All_v1_Gauss250",#
  #"Species_Betu_All_v1_Gauss750",#
  #"Species_Betu_Pap_v1",#
  #"Species_Betu_Pap_v1_Gauss250",#
  "Species_Betu_Pap_v1_Gauss750",#
  #"Species_Fagu_Gra_v1",#
  "Species_Fagu_Gra_v1_Gauss250",#
  #"Species_Fagu_Gra_v1_Gauss750",#
  #"Species_Frax_Ame_v1",#
  #"Species_Frax_Ame_v1_Gauss250",#
  "Species_Frax_Ame_v1_Gauss750",#
  #"Species_Lari_Lar_v1",#
  #"Species_Lari_Lar_v1_Gauss250",#
  "Species_Lari_Lar_v1_Gauss750",#
  "Species_Pice_Gla_v1",#
  #"Species_Pice_Gla_v1_Gauss250",#
  #"Species_Pice_Gla_v1_Gauss750",#
  #"Species_Pice_Rub_v1",#
  "Species_Pice_Rub_v1_Gauss250",#
  #"Species_Pice_Rub_v1_Gauss750",#
  #"Species_Pice_Mar_v1",#
  "Species_Pice_Mar_v1_Gauss250",#
  #"Species_Pice_Mar_v1_Gauss750",#
  #"Species_Pinu_Ban_v1",#
  #"Species_Pinu_Ban_v1_Gauss250",#
  "Species_Pinu_Ban_v1_Gauss750",#
  "Species_Pinu_Res_v1",#
  #"Species_Pinu_Res_v1_Gauss250",#
  #"Species_Pinu_Res_v1_Gauss750",#
  #"Species_Pinu_Str_v1",#
  #"Species_Pinu_Str_v1_Gauss250",#
  "Species_Pinu_Str_v1_Gauss750",#
  "Species_Popu_Bal_v1",#
  #"Species_Popu_Bal_v1_Gauss250",#
  #"Species_Popu_Bal_v1_Gauss750",#
  #"Species_Quer_Rub_v1",#
  #"Species_Quer_Rub_v1_Gauss250",#
  "Species_Quer_Rub_v1_Gauss750",#
  #"Species_Tsug_Can_v1",#
  #"Species_Tsug_Can_v1_Gauss250",#
  "Species_Tsug_Can_v1_Gauss750",#
  #"Structure_Biomass_TotalLiveAboveGround_v1",
  "Structure_Biomass_TotalLiveAboveGround_v1_Gauss250",
  #"Structure_Biomass_TotalLiveAboveGround_v1_Gauss750",
  #"Structure_Stand_Age_v1",
  "Structure_Stand_Age_v1_Gauss250",
  #"Structure_Stand_Age_v1_Gauss750",
  #"wtbl250",
  "wtbl250_modal3x3",
  #"wtbl250_modal5x5",
  #"CTI250",
  "CTI250_Gauss250",
  #"CTI250_Gauss750",
  "roadsNNS",
  #"roadsNNS_Gauss250",
  #"roadsNNS_Gauss750",
  #"urbagNNS",
  #"urbagNNS_Gauss250",
  "urbagNNS_Gauss750",
  #"watNNS",
  "watNNS_Gauss250",
  #,"watNNS_Gauss750",
  "ARU"
)

# create object storing the indices of predictor layers (from the prediction dataset) for the autocorrelation assessment that will inform block size. Only include continuous numeric variables here. Can use indexing "[-c(x:y)] to exclude them (e.g. exclude climate variables, HF, etc).
ind2 <- which(names(pred_abs_2011) %in% pred.variables2[-21])

# Calculate autocorrelation range
start_time <- Sys.time()
sp.auto.arr2 <- spatialAutoRange(pred_abs_2011[[ind2]])
end_time <- Sys.time()
end_time - start_time

# Use spatial blocks to separate train and test folds
start_time <- Sys.time()
set.seed(123)
sp_block_select <-  spatialBlock(
  speciesData = datcombo_sp,
  species = "OCCU",
  rasterLayer = pred_abs_2011[[1]],
  iteration = 250,
  theRange = sp.auto.arr2$range,
  selection = "random",
  maskBySpecies = FALSE,
  k=5
)
end_time <- Sys.time()
end_time - start_time

start_time<-Sys.time()
brt2<- brt_blocks(data=datcombo,pred.variables = pred.variables2, lr=0.01,tc=3, output.folder = "D://CHID regional NS BRT/BRT_outputs/selected_scales_new/", blocks=sp_block_select, save.points.shp = TRUE)  
end_time<-Sys.time()
end_time-start_time


# Cell level (250m, no Gaussian filter) ####
pred.variables3<-c( 
  "Species_Abie_Bal_v1",#
  #"Species_Abie_Bal_v1_Gauss250",#
  #"Species_Abie_Bal_v1_Gauss750",#
  "Species_Acer_Sah_v1",#
  #"Species_Acer_Sah_v1_Gauss250",#
  #"Species_Acer_Sah_v1_Gauss750",#
  "Species_Acer_Rub_v1",#
  #"Species_Acer_Rub_v1_Gauss250",#
  #"Species_Acer_Rub_v1_Gauss750",#
  "Species_Acer_Spi_v1",#
  #"Species_Acer_Spi_v1_Gauss250",#
  #"Species_Acer_Spi_v1_Gauss750",#
  "Species_Betu_All_v1",#
  #"Species_Betu_All_v1_Gauss250",#
  #"Species_Betu_All_v1_Gauss750",#
  "Species_Betu_Pap_v1",#
  #"Species_Betu_Pap_v1_Gauss250",#
  #"Species_Betu_Pap_v1_Gauss750",#
  "Species_Fagu_Gra_v1",#
  #"Species_Fagu_Gra_v1_Gauss250",#
  #"Species_Fagu_Gra_v1_Gauss750",#
  "Species_Frax_Ame_v1",#
  #"Species_Frax_Ame_v1_Gauss250",#
  #"Species_Frax_Ame_v1_Gauss750",#
  "Species_Lari_Lar_v1",#
  #"Species_Lari_Lar_v1_Gauss250",#
  #"Species_Lari_Lar_v1_Gauss750",#
  "Species_Pice_Gla_v1",#
  #"Species_Pice_Gla_v1_Gauss250",#
  #"Species_Pice_Gla_v1_Gauss750",#
  "Species_Pice_Rub_v1",#
  #"Species_Pice_Rub_v1_Gauss250",#
  #"Species_Pice_Rub_v1_Gauss750",#
  "Species_Pice_Mar_v1",#
  #"Species_Pice_Mar_v1_Gauss250",#
  #"Species_Pice_Mar_v1_Gauss750",#
  "Species_Pinu_Ban_v1",#
  #"Species_Pinu_Ban_v1_Gauss250",#
  #"Species_Pinu_Ban_v1_Gauss750",#
  "Species_Pinu_Res_v1",#
  #"Species_Pinu_Res_v1_Gauss250",#
  #"Species_Pinu_Res_v1_Gauss750",#
  "Species_Pinu_Str_v1",#
  #"Species_Pinu_Str_v1_Gauss250",#
  #"Species_Pinu_Str_v1_Gauss750",#
  "Species_Popu_Bal_v1",#
  #"Species_Popu_Bal_v1_Gauss250",#
  #"Species_Popu_Bal_v1_Gauss750",#
  "Species_Quer_Rub_v1",#
  #"Species_Quer_Rub_v1_Gauss250",#
  #"Species_Quer_Rub_v1_Gauss750",#
  "Species_Tsug_Can_v1",#
  #"Species_Tsug_Can_v1_Gauss250",#
  #"Species_Tsug_Can_v1_Gauss750",#
  "Structure_Biomass_TotalLiveAboveGround_v1",
  #"Structure_Biomass_TotalLiveAboveGround_v1_Gauss250",
  #"Structure_Biomass_TotalLiveAboveGround_v1_Gauss750",
  "Structure_Stand_Age_v1",
  #"Structure_Stand_Age_v1_Gauss250",
  #"Structure_Stand_Age_v1_Gauss750",
  "wtbl250",
  #"wtbl250_modal3x3",
  #"wtbl250_modal5x5",
  "CTI250",
  #"CTI250_Gauss250",
  #"CTI250_Gauss750",
  "roadsNNS",
  #"roadsNNS_Gauss250",
  #"roadsNNS_Gauss750",
  "urbagNNS",
  #"urbagNNS_Gauss250",
  #"urbagNNS_Gauss750",
  "watNNS",
  #"watNNS_Gauss250",
  #"watNNS_Gauss750",
  "ARU"
)

# create object storing the indices of predictor layers (from the prediction dataset) for the autocorrelation assessment that will inform block size. Only include continuous numeric variables here. Can use indexing "[-c(x:y)] to exclude them (e.g. exclude climate variables, HF, etc).
ind3 <- which(names(pred_abs_2011) %in% pred.variables3[-21])

# Calculate autocorrelation range
start_time <- Sys.time()
sp.auto.arr3 <- spatialAutoRange(pred_abs_2011[[ind3]])
end_time <- Sys.time()
end_time - start_time

# Use spatial blocks to separate train and test folds
start_time <- Sys.time()
set.seed(123)
sp_block_cell <-  spatialBlock(
  speciesData = datcombo_sp,
  species = "OCCU",
  rasterLayer = pred_abs_2011[[1]],
  iteration = 250,
  theRange = sp.auto.arr3$range,
  selection = "random",
  maskBySpecies = FALSE,
  k=5
)
end_time <- Sys.time()
end_time - start_time

start_time<-Sys.time()
brt3<- brt_blocks(data=datcombo,pred.variables = pred.variables3, lr=0.00001,tc=2, output.folder = "D://CHID regional NS BRT/BRT_outputs/cell_level_new/", blocks=sp_block_cell, save.points.shp = FALSE)  
end_time<-Sys.time()
end_time-start_time

# using these predictors gbm.step does not run even with lr=0.000001 and tc=2

# Only Gaussian filter with sigma = 250m ####
pred.variables4<-c( 
  #"Species_Abie_Bal_v1",#
  "Species_Abie_Bal_v1_Gauss250",#
  #"Species_Abie_Bal_v1_Gauss750",#
  #"Species_Acer_Sah_v1",#
  "Species_Acer_Sah_v1_Gauss250",#
  #"Species_Acer_Sah_v1_Gauss750",#
  #"Species_Acer_Rub_v1",#
  "Species_Acer_Rub_v1_Gauss250",#
  #"Species_Acer_Rub_v1_Gauss750",#
  #"Species_Acer_Spi_v1",#
  "Species_Acer_Spi_v1_Gauss250",#
  #"Species_Acer_Spi_v1_Gauss750",#
  #"Species_Betu_All_v1",#
  "Species_Betu_All_v1_Gauss250",#
  #"Species_Betu_All_v1_Gauss750",#
  #"Species_Betu_Pap_v1",#
  "Species_Betu_Pap_v1_Gauss250",#
  #"Species_Betu_Pap_v1_Gauss750",#
  #"Species_Fagu_Gra_v1",#
  "Species_Fagu_Gra_v1_Gauss250",#
  #"Species_Fagu_Gra_v1_Gauss750",#
  #"Species_Frax_Ame_v1",#
  "Species_Frax_Ame_v1_Gauss250",#
  #"Species_Frax_Ame_v1_Gauss750",#
  #"Species_Lari_Lar_v1",#
  "Species_Lari_Lar_v1_Gauss250",#
  #"Species_Lari_Lar_v1_Gauss750",#
  #"Species_Pice_Gla_v1",#
  "Species_Pice_Gla_v1_Gauss250",#
  #"Species_Pice_Gla_v1_Gauss750",#
  #"Species_Pice_Rub_v1",#
  "Species_Pice_Rub_v1_Gauss250",#
  #"Species_Pice_Rub_v1_Gauss750",#
  #"Species_Pice_Mar_v1",#
  "Species_Pice_Mar_v1_Gauss250",#
  #"Species_Pice_Mar_v1_Gauss750",#
  #"Species_Pinu_Ban_v1",#
  "Species_Pinu_Ban_v1_Gauss250",#
  #"Species_Pinu_Ban_v1_Gauss750",#
  #"Species_Pinu_Res_v1",#
  "Species_Pinu_Res_v1_Gauss250",#
  #"Species_Pinu_Res_v1_Gauss750",#
  #"Species_Pinu_Str_v1",#
  "Species_Pinu_Str_v1_Gauss250",#
  #"Species_Pinu_Str_v1_Gauss750",#
  #"Species_Popu_Bal_v1",#
  "Species_Popu_Bal_v1_Gauss250",#
  #"Species_Popu_Bal_v1_Gauss750",#
  #"Species_Quer_Rub_v1",#
  "Species_Quer_Rub_v1_Gauss250",#
  #"Species_Quer_Rub_v1_Gauss750",#
  #"Species_Tsug_Can_v1",#
  "Species_Tsug_Can_v1_Gauss250",#
  #"Species_Tsug_Can_v1_Gauss750",#
  #"Structure_Biomass_TotalLiveAboveGround_v1",
  "Structure_Biomass_TotalLiveAboveGround_v1_Gauss250",
  #"Structure_Biomass_TotalLiveAboveGround_v1_Gauss750",
  #"Structure_Stand_Age_v1",
  "Structure_Stand_Age_v1_Gauss250",
  #"Structure_Stand_Age_v1_Gauss750",
  #"wtbl250",
  "wtbl250_modal3x3",
  #"wtbl250_modal5x5",
  "CTI250",
  #"CTI250_Gauss250",
  #"CTI250_Gauss750",
  #"roadsNNS",
  "roadsNNS_Gauss250",
  #"roadsNNS_Gauss750",
  #"urbagNNS",
  "urbagNNS_Gauss250",
  #"urbagNNS_Gauss750",
  #"watNNS",
  "watNNS_Gauss250",
  #"watNNS_Gauss750",
  "ARU"
)

# create object storing the indices of predictor layers (from the prediction dataset) for the autocorrelation assessment that will inform block size. Only include continuous numeric variables here. Can use indexing "[-c(x:y)] to exclude them (e.g. exclude climate variables, HF, etc).
ind4 <- which(names(pred_abs_2011) %in% pred.variables4[-21])

# Calculate autocorrelation range
start_time <- Sys.time()
sp.auto.arr4 <- spatialAutoRange(pred_abs_2011[[ind4]])
end_time <- Sys.time()
end_time - start_time

# Use spatial blocks to separate train and test folds
start_time <- Sys.time()
set.seed(123)
sp_block_GF250 <-  spatialBlock(
  speciesData = datcombo_sp,
  species = "OCCU",
  rasterLayer = pred_abs_2011[[1]],
  iteration = 250,
  theRange = sp.auto.arr4$range,
  selection = "random",
  maskBySpecies = FALSE,
  k=5
)
end_time <- Sys.time()
end_time - start_time

start_time<-Sys.time()
brt4<- brt_blocks(data=datcombo,pred.variables = pred.variables4, lr=0.01,tc=3, output.folder = "D://CHID regional NS BRT/BRT_outputs/GFsigma250m_new/", blocks=sp_block_GF250, save.points.shp = TRUE)  
end_time<-Sys.time()
end_time-start_time


# Only Gaussian filter with sigma = 750m ####
pred.variables5<-c( 
  #"Species_Abie_Bal_v1",#
  #"Species_Abie_Bal_v1_Gauss250",#
  "Species_Abie_Bal_v1_Gauss750",#
  #"Species_Acer_Sah_v1",#
  #"Species_Acer_Sah_v1_Gauss250",#
  "Species_Acer_Sah_v1_Gauss750",#
  #"Species_Acer_Rub_v1",#
  #"Species_Acer_Rub_v1_Gauss250",#
  "Species_Acer_Rub_v1_Gauss750",#
  #"Species_Acer_Spi_v1",#
  #"Species_Acer_Spi_v1_Gauss250",#
  "Species_Acer_Spi_v1_Gauss750",#
  #"Species_Betu_All_v1",#
  #"Species_Betu_All_v1_Gauss250",#
  "Species_Betu_All_v1_Gauss750",#
  #"Species_Betu_Pap_v1",#
  #"Species_Betu_Pap_v1_Gauss250",#
  "Species_Betu_Pap_v1_Gauss750",#
  #"Species_Fagu_Gra_v1",#
  #"Species_Fagu_Gra_v1_Gauss250",#
  "Species_Fagu_Gra_v1_Gauss750",#
  #"Species_Frax_Ame_v1",#
  #"Species_Frax_Ame_v1_Gauss250",#
  "Species_Frax_Ame_v1_Gauss750",#
  #"Species_Lari_Lar_v1",#
  #"Species_Lari_Lar_v1_Gauss250",#
  "Species_Lari_Lar_v1_Gauss750",#
  #"Species_Pice_Gla_v1",#
  #"Species_Pice_Gla_v1_Gauss250",#
  "Species_Pice_Gla_v1_Gauss750",#
  #"Species_Pice_Rub_v1",#
  #"Species_Pice_Rub_v1_Gauss250",#
  "Species_Pice_Rub_v1_Gauss750",#
  #"Species_Pice_Mar_v1",#
  #"Species_Pice_Mar_v1_Gauss250",#
  "Species_Pice_Mar_v1_Gauss750",#
  #"Species_Pinu_Ban_v1",#
  #"Species_Pinu_Ban_v1_Gauss250",#
  "Species_Pinu_Ban_v1_Gauss750",#
  #"Species_Pinu_Res_v1",#
  #"Species_Pinu_Res_v1_Gauss250",#
  "Species_Pinu_Res_v1_Gauss750",#
  #"Species_Pinu_Str_v1",#
  #"Species_Pinu_Str_v1_Gauss250",#
  "Species_Pinu_Str_v1_Gauss750",#
  #"Species_Popu_Bal_v1",#
  #"Species_Popu_Bal_v1_Gauss250",#
  "Species_Popu_Bal_v1_Gauss750",#
  #"Species_Quer_Rub_v1",#
  #"Species_Quer_Rub_v1_Gauss250",#
  "Species_Quer_Rub_v1_Gauss750",#
  #"Species_Tsug_Can_v1",#
  #"Species_Tsug_Can_v1_Gauss250",#
  "Species_Tsug_Can_v1_Gauss750",#
  #"Structure_Biomass_TotalLiveAboveGround_v1",
  #"Structure_Biomass_TotalLiveAboveGround_v1_Gauss250",
  "Structure_Biomass_TotalLiveAboveGround_v1_Gauss750",
  #"Structure_Stand_Age_v1",
  #"Structure_Stand_Age_v1_Gauss250",
  "Structure_Stand_Age_v1_Gauss750",
  #"wtbl250",
  #"wtbl250_modal3x3",
  "wtbl250_modal5x5",
  #"CTI250",
  #"CTI250_Gauss250",
  "CTI250_Gauss750",
  #"roadsNNS",
  #"roadsNNS_Gauss250",
  "roadsNNS_Gauss750",
  #"urbagNNS",
  #"urbagNNS_Gauss250",
  "urbagNNS_Gauss750",
  #"watNNS",
  #"watNNS_Gauss250",
  "watNNS_Gauss750",
  "ARU"
)

# create object storing the indices of predictor layers (from the prediction dataset) for the autocorrelation assessment that will inform block size. Only include continuous numeric variables here. Can use indexing "[-c(x:y)] to exclude them (e.g. exclude climate variables, HF, etc).
ind5 <- which(names(pred_abs_2011) %in% pred.variables5[-21])

# Calculate autocorrelation range
start_time <- Sys.time()
sp.auto.arr5 <- spatialAutoRange(pred_abs_2011[[ind5]])
end_time <- Sys.time()
end_time - start_time

# Use spatial blocks to separate train and test folds
start_time <- Sys.time()
set.seed(123)
sp_block_GF750 <-  spatialBlock(
  speciesData = datcombo_sp,
  species = "OCCU",
  rasterLayer = pred_abs_2011[[1]],
  iteration = 250,
  theRange = sp.auto.arr5$range,
  selection = "random",
  maskBySpecies = FALSE,
  k=5
)
end_time <- Sys.time()
end_time - start_time

start_time<-Sys.time()
brt5<- brt_blocks(data=datcombo,pred.variables = pred.variables5, lr=0.01,tc=3, output.folder = "D://CHID regional NS BRT/BRT_outputs/GFsigma750m_new/", blocks=sp_block_GF750, save.points.shp = TRUE)  
end_time<-Sys.time()
end_time-start_time

save.image("D:/CHID regional NS BRT/BRT_outputs/outputs.RData")

# Comparison of the different models ####
library(ggplot2)


# Model deviance

df_devplot<- data.frame(grp=c("Full model",
                      "Selected scales",
                      #"cell level",
                      "250m Gaussian smooth",
                      "750m Gaussian smooth"),
                deviance=c(brt1$cv.statistics$deviance.mean,
                           brt2$cv.statistics$deviance.mean,
                           #brt3$cv.statistics$deviance.mean,
                           brt4$cv.statistics$deviance.mean,
                           brt5$cv.statistics$deviance.mean),
                se=c(brt1$cv.statistics$deviance.se,
                     brt2$cv.statistics$deviance.se,
                     #brt3$cv.statistics$deviance.se,
                     brt4$cv.statistics$deviance.se,
                     brt5$cv.statistics$deviance.se)
)
df_devplot
k<-ggplot(df_devplot, aes(grp,deviance,ymin=deviance-se,ymax=deviance+se))
k+geom_pointrange()

# correlation
df_corrplot<- data.frame(grp=c("Full model",
                       "Selected scales",
                       #"cell level",
                       "250m Gaussian smooth",
                       "750m Gaussian smooth"),
                 correlation=c(brt1$cv.statistics$correlation.mean,
                               brt2$cv.statistics$correlation.mean,
                               #brt3$cv.statistics$correlation.mean,
                               brt4$cv.statistics$correlation.mean,
                               brt5$cv.statistics$correlation.mean),
                 se=c(brt1$cv.statistics$correlation.se,
                      brt2$cv.statistics$correlation.se,
                      #brt3$cv.statistics$correlation.se,
                      brt4$cv.statistics$correlation.se,
                      brt5$cv.statistics$correlation.se)
)
df_corrplot
k2<-ggplot(df_corrplot, aes(grp,correlation,ymin=correlation-se,ymax=correlation+se))
k2+geom_pointrange()

# Calibration:  The first two statistics were the estimated intercepts and slopes of linear regression models of predictions against observations. The intercept measures the magnitude and direction of bias, with values close to 0 indicating low or no bias. The slope yields information about the consistency in the bias as a function of the mean, with a value of 1 indicating a consistent bias if the intercept is a nonzero value.
## intercept
df_calli_int_plot<- data.frame(grp=c("Full model",
                       "Selected scales",
                       #"cell level",
                       "250m Gaussian smooth",
                       "750m Gaussian smooth"),
                 calibration.intercept=c(brt1$cv.statistics$calibration.mean[1],
                                         brt2$cv.statistics$calibration.mean[1],
                                         #brt3$cv.statistics$calibration.mean[1],
                                         brt4$cv.statistics$calibration.mean[1],
                                         brt5$cv.statistics$calibration.mean[1]),
                 se=c(brt1$cv.statistics$calibration.se[1],
                      brt2$cv.statistics$calibration.se[1],
                      #brt3$cv.statistics$calibration.se[1],
                      brt4$cv.statistics$calibration.se[1],
                      brt5$cv.statistics$calibration.se[1])
)
df_calli_int_plot
k3<-ggplot(df_calli_int_plot, aes(grp,calibration.intercept,ymin=calibration.intercept-se,ymax=calibration.intercept+se))
k3+geom_pointrange()

## slope
df_calli_slope_plot<- data.frame(grp=c("Full model",
                       "Selected scales",
                       #"cell level",
                       "250m Gaussian smooth",
                       "750m Gaussian smooth"),
                 calibration.slope=c(brt1$cv.statistics$calibration.mean[2],
                                     brt2$cv.statistics$calibration.mean[2],
                                     #brt3$cv.statistics$calibration.mean[2],
                                     brt4$cv.statistics$calibration.mean[2],
                                     brt5$cv.statistics$calibration.mean[2]),
                 se=c(brt1$cv.statistics$calibration.se[2],
                      brt2$cv.statistics$calibration.se[2],
                      #brt3$cv.statistics$calibration.se[2],
                      brt4$cv.statistics$calibration.se[2],
                      brt5$cv.statistics$calibration.se[2])
)
df_calli_slope_plot
k4<-ggplot(df_calli_slope_plot, aes(grp,calibration.slope,ymin=calibration.slope-se,ymax=calibration.slope+se))
k4+geom_pointrange()


# Plotting deviance for each brt model ####
dev_plot<-function(brt){
  m<-brt
  y.bar <- min(m$cv.values) 
  y.min <- min(m$cv.values - m$cv.loss.ses)
  y.max <- max(m$cv.values + m$cv.loss.ses)
  
  plot(m$trees.fitted, m$cv.values, type = 'l', ylab = "Holdout deviance", xlab = "no. of trees", ylim = c(y.min,y.max))
  abline(h = y.bar, col = 3)
  
  lines(m$trees.fitted, m$cv.values + m$cv.loss.ses, lty=2)  
  lines(m$trees.fitted, m$cv.values - m$cv.loss.ses, lty=2)  
  
  target.trees <- m$trees.fitted[match(TRUE,m$cv.values == y.bar)]
  abline(v = target.trees, col=4)
}

dev_plot(brt5)


####
# A function to obtain confidence intervals of model predictions ####
boot_brt<-function(data,brtmodel,blocks,pred.data,nsamples=100,output.folder){
  rast <-  predict(pred.data,
                   brtmodel,
                   type = "response",
                   n.trees = brtmodel$n.trees)
  stack<-stack(rast)
  names(stack)[1]<-"model estimate"
  
  z <- output.folder
  
  if (file.exists(z) == FALSE) {
    dir.create(z)
  }
  
  for(i in 1:nsamples){
    
    
    cat("loop",i,"\n") # this prints the loop number on console to track function progress
    
   
    brt <- NULL
    attempt <- 0
    while( is.null(brt) && attempt <= 20 ) {
      attempt <- attempt + 1
      
      sample<-sample(1:nrow(data),size=nrow(data),replace=T)
      data2<-data[sample,]
      
      try(brt <-
            gbm.step(
              data = data2,
              gbm.y = brtmodel$gbm.call$gbm.y,
              gbm.x = brtmodel$gbm.call$gbm.x,
              fold.vector = blocks$foldID[sample],
              n.folds = brtmodel$gbm.call$cv.folds,
              family = brtmodel$gbm.call$family,
              tree.complexity = brtmodel$gbm.call$tree.complexity,
              learning.rate = brtmodel$gbm.call$learning.rate,
              bag.fraction = brtmodel$gbm.call$bag.fraction,
              offset = brtmodel$gbm.call$offset[sample],
              site.weights = data2$wt
            ))
    } 
    
    

    
    
    rast <- predict(pred.data,
                    brt,
                    type = "response",
                    n.trees = brt$n.trees)
    
    
    
    stack <- addLayer(stack, rast)
    names(stack)[i+1]<-paste0("sample",i) 
    
    
    saveRDS(brt, file=paste0(z,"bootstrap_model",i,".R"))
    
    
  }
  
  fun0.05 <- function(x) {quantile(x, probs = 0.05, na.rm = TRUE)}
  lower<- calc(stack[[-1]],fun0.05)
  fun0.95 <- function(x) {quantile(x, probs = 0.95, na.rm = TRUE)}
  upper<- calc(stack[[-1]],fun0.95)
  
  writeRaster(lower, filename=paste0(z, "confint_lower.tif"), format="GTiff",overwrite=TRUE)
  writeRaster(upper, filename=paste0(z, "confint_upper.tif"), format="GTiff",overwrite=TRUE)
  writeRaster(stack, filename=paste0(z, "samples.grd"), format="raster",overwrite=TRUE)
  
  return(stack)
}

start_time <- Sys.time()
confintbrt5<-boot_brt(datcombo,brt5,sp_block_GF750,pred_abs_2011,nsamples=250,output.folder = "D://CHID regional NS BRT/BRT_outputs/GFsigma750m_new/confint/")
end_time <- Sys.time()
end_time - start_time


save.image("D:/CHID regional NS BRT/BRT_outputs/outputs.RData")

writeRaster(confintbrt5,file="D://CHID regional NS BRT/BRT_outputs/GFsigma750m_new/confint/samples.grd", format="raster",overwrite=T)

names(confintbrt5)





# Population size from density predictions ####
#pred_abs_2011<- brick("D:/CHID regional NS BRT/prediction dataset/abs2011_250m.grd")
load("D:/CHID regional NS BRT/BRT_outputs/outputs.RData")

library(gbm)
start_time <- Sys.time()
rast <-  predict(pred_abs_2011,
                 brt5,
                 type = "response",
                 n.trees = brt5$n.trees)
end_time <- Sys.time()
end_time - start_time
plot(rast) # density is singing males/ha

# need to consider only the raster cells that went into the simulations
load("D:/CHID regional NS BRT/Landscape simulation results/Prediction raster/birddensityRasters_NovaScotia_baseline_BudwormBaselineFire_1_CAWA.RData")
plot(PredictBirdsp.l[[1]])

rast2<-mask(rast,PredictBirdsp.l[[1]])

# Multiplying density values by 6.25 (hectars in a 250-250 pixel) to get abundance for each cell and overall population size

density_brt5_250m<-rast2*6.25
popsize<-cellStats(density_brt5_250m,stat=sum,na.rm=T) # THIS IS ESTIMATED POPULATION SIZE


density_brt5_250m_full<-rast*6.25
popsize_full<-cellStats(density_brt5_250m_full,stat=sum,na.rm=T)


# alternative is to obtain 100 random Poisson draws using cell density as rate parameter:
v1<-values(density_brt5_250m)[!is.na(values(density_brt5_250m))]
abundance_brt5<-matrix(NA,length(v1),100)
for(i in 1:nrow(abundance_brt5)){
  abundance_brt5[i,]<- rpois(100,v1[i])
}
## mean of sums across the 100 draws
mean(colSums(abundance_brt5)) 


# Same procedure for upper and lower CI predictions  
upper<-raster("D:/CHID regional NS BRT/BRT_outputs/GFsigma750m_new/confint/confint_upper.tif") 
upper2<-mask(upper,PredictBirdsp.l[[1]])
upper_brt5_250m<-upper2*6.25
popsize.upp<-cellStats(upper_brt5_250m,stat=sum,na.rm=T)

upper_brt5_250m_full<-upper*6.25
popsize.upp_full<-cellStats(upper_brt5_250m_full,stat=sum,na.rm=T)

# v2<-values(upper_brt5_250m)[!is.na(values(upper_brt5_250m))]
# abundance_upper_brt5<-matrix(NA,length(v2),100)
# for(i in 1:nrow(abundance_upper_brt5)){
#   abundance_upper_brt5[i,]<- rpois(100,v2[i])
# }
# mean(colSums(abundance_upper_brt5))

lower<-raster("D://CHID regional NS BRT/BRT_outputs/GFsigma750m_new/confint/confint_lower.tif")
lower2<-mask(lower,PredictBirdsp.l[[1]])
lower_brt5_250m<-lower2*6.25
popsize.low<-cellStats(lower_brt5_250m,stat=sum,na.rm=T)

lower_brt5_250m_full<-lower*6.25
popsize.low_full<-cellStats(lower_brt5_250m_full,stat=sum,na.rm=T)
# v3<-values(lower_brt5_250m)[!is.na(values(lower_brt5_250m))]
# abundance_lower_brt5<-matrix(NA,length(v3),100)
# for(i in 1:nrow(abundance_lower_brt5)){
#   abundance_lower_brt5[i,]<- rpois(100,v3[i])
# }
# mean(colSums(abundance_lower_brt5))

# to obtain custom CI (other than 90%):
upper3<- calc(confintbrt5[[-1]],function(x) {quantile(x, probs = 0.95, na.rm = TRUE)})
upper3<-mask(upper3,PredictBirdsp.l[[1]])
upper3_250m<-upper3*6.25
popsize.upp<-cellStats(upper3_250m,stat=sum,na.rm=T)

lower3<- calc(confintbrt5[[-1]],function(x) {quantile(x, probs = 0.05, na.rm = TRUE)})
lower3<-mask(lower3,PredictBirdsp.l[[1]])
lower3_250m<-lower3*6.25
popsize.low<-cellStats(lower3_250m,stat=sum,na.rm=T)



# cell level SDs
bootsamples<-confintbrt5[[-1]]
names(bootsamples)
bootsamples_masked<-mask(bootsamples,PredictBirdsp.l[[1]])
SDpreds<-calc(bootsamples_masked,fun=sd,na.rm=T)


# there are some outliers with very high SD. Set those values to be equal to the 99% quantile
#values(SDpreds)[which(values(SDpreds)>quantile(values(SDpreds),0.9995,na.rm=T))]<-quantile(values(SDpreds),0.9995,na.rm=T)

plot(SDpreds)

popsizeCIs<-c(estimate=popsize,lower=popsize.low,upper=popsize.upp)

# Current population size probability  ####

bootpopsizes<-rep(NA,250)
for(i in 1:250){
  bootpopsizes[i]<-cellStats(bootsamples_masked[[i]]*6.25,stat=sum,na.rm=T)
}


bootpopsizes

mean(bootpopsizes)
median(bootpopsizes)
sd(bootpopsizes)  
###

#with normal distribution:
poprange<-seq(15000,70000,by=1000)
id.dist<-pnorm(poprange,mean=mean(bootpopsizes),sd=sd(bootpopsizes),lower.tail = F)
id.df<-data.frame("Pop.range"=poprange,"Probability"=id.dist,Distribution="normal")
library(ggplot2)
ggplot(id.df, aes(x = Pop.range, y = Probability)) + geom_point()+ geom_vline(xintercept = popsize,  color = "blue", size=0.8)

#with gamma distribution:
findGamma1<-function(scale.init,shape.init,low,upp){
  mat<-expand.grid(scale=runif(1000,min=scale.init[1],max=scale.init[2]),shape=runif(1000,min=shape.init[1],max=shape.init[2]))
  mat[,3:4]<-t(apply(mat,1,FUN = function(x){qgamma(p=c(0.05,0.95),shape=x[2],scale=x[1])}))
  mat$distances<-pointDistance(p1=mat[,3:4],p2=c(low,upp),lonlat = F)
  out1<-mat[which.min(mat$distances),]
  
  scale.init2<-unlist(c(out1[1]-(out1[1]*0.25),out1[1]+(out1[1]*0.25)))
  shape.init2<-unlist(c(out1[2]-(out1[2]*0.25),out1[2]+(out1[2]*0.25)))
  
  mat2<-expand.grid(scale=runif(1000,min=scale.init2[1],max=scale.init2[2]),shape=runif(1000,min=shape.init2[2],max=shape.init2[2]))
  mat2[,3:4]<-t(apply(mat2,1,FUN = function(x){qgamma(p=c(0.05,0.95),shape=x[2],scale=x[1])}))
  mat2$distances<-pointDistance(p1=mat2[,3:4],p2=c(low,upp),lonlat = F)
  
  mat3<-rbind(mat,mat2)
  colnames(mat3)[3:4]<-c("lowerCL_0.5","upperCL_0.95")
  return(mat3[which.min(mat3$distances),])
}


GammaParamsCurrent<-findGamma1(scale.init=c(2000,60000),shape.init = c(0.25,15),upp=popsize.upp,low=popsize.low)
GammaParamsCurrent


id.dist.gamma<-pgamma(poprange,shape=GammaParamsCurrent$shape,scale=GammaParamsCurrent$scale,lower.tail = F)
id.df.gamma<-data.frame("Pop.range"=poprange,"Probability"=id.dist.gamma,Distribution="gamma")


ycurrentprednormal<-pnorm(popsize,mean=mean(bootpopsizes),sd=sd(bootpopsizes),lower.tail = F)
ycurrentpredgamma<-pgamma(popsize,shape=GammaParamsCurrent$shape,scale=GammaParamsCurrent$scale,lower.tail = F)


library(ggplot2)
ggplot(rbind(id.df,id.df.gamma)) + 
  geom_point(aes(x = Pop.range, y = Probability,color=Distribution))+ 
  geom_errorbarh(aes(xmax=popsize.upp,xmin=popsize.low,y=ycurrentpredgamma,height=0.03),color="#00BFC4")+
  geom_point(x = popsize,y=ycurrentpredgamma)+
  geom_errorbarh(aes(xmax=popsize.upp,xmin=popsize.low,y=ycurrentprednormal,height=0.03),color="#F8766D")+
  geom_point(x = popsize,y=ycurrentprednormal)+  
  labs(x="N",y="Probability (population size >= N)")+
  annotate("text", x= 38000, y = 0.57,hjust=0, label = paste0("Model prediction (90% CI) = ",round(popsize,0)," (",round(popsize.low,0)," - ",round(popsize.upp,0),")"),size=4)

ggplot(rbind(id.df)) + 
  geom_point(aes(x = Pop.range, y = Probability),color="#F8766D")+ 
  #geom_errorbarh(aes(xmax=popsize.upp,xmin=popsize.low,y=ycurrentpredgamma,height=0.03),color="#00BFC4")+
  #geom_point(x = popsize,y=ycurrentpredgamma)+
  geom_errorbarh(aes(xmax=popsize.upp,xmin=popsize.low,y=ycurrentprednormal,height=0.03),color="#F8766D")+
  geom_point(x = popsize,y=ycurrentprednormal)+  
  labs(x="N",y="Probability (population size >= N)")+
  annotate("text", x= 38000, y = 0.57,hjust=0, label = paste0("Model prediction (90% CI) = ",round(popsize,0)," (",round(popsize.low,0)," - ",round(popsize.upp,0),")"),size=4)

ggplot(rbind(id.df.gamma)) + 
  geom_point(aes(x = Pop.range, y = Probability),color="#00BFC4")+ 
  geom_errorbarh(aes(xmax=popsize.upp,xmin=popsize.low,y=ycurrentpredgamma,height=0.03),color="#00BFC4")+
  geom_point(x = popsize,y=ycurrentpredgamma)+
  #geom_errorbarh(aes(xmax=popsize.upp,xmin=popsize.low,y=ycurrentprednormal,height=0.03),color="#F8766D")+
  #geom_point(x = popsize,y=ycurrentprednormal)+  
  labs(x="N",y="Probability (population size >= N)")+
  annotate("text", x= 38000, y = 0.57,hjust=0, label = paste0("Model prediction (90% CI) = ",round(popsize,0)," (",round(popsize.low,0)," - ",round(popsize.upp,0),")"),size=4)

#### Future Population size - current conditions ####
# trend from BBS (Smith et al 2014), CAWA in Nova Scotia
trends<- c(long_trend = -1.78, # long term trend (>40 years)
                    long_trend_lower= -3.56,
                    long_trend_upper= 0.124,
                    short_trend = 5.46, # short term trend (10 years) 
                    short_trend_lower =-3.31,
                    short_trend_upper = 16)

# Population growth model - also from Smith et al 2014

popgrowth<-function(data){
  years<-data[1]
  N0<-data[2]
  trend<-data[3]
  lower<-data[4]
  upper<-data[5]

  est<-N0*(1+trend/100)^years
  low<-N0*(1+lower/100)^years
  upp<-N0*(1+upper/100)^years
  out<-c(est=est,low=low,upp=upp)
  out
  }

year.range<-23

df<-data.frame(years=rep(seq(1:year.range),6),N0=rep(popsizeCIs,each=2*year.range),trend=c(rep(trends[1],year.range),rep(trends[4],year.range)),lower=c(rep(trends[2],year.range),rep(trends[5],year.range)),upper=c(rep(trends[3],year.range),rep(trends[6],year.range)))
df[,6:8]<-t(apply(df,1,popgrowth))
names(df)[6:8]<-c("Ni","Ni_lower", "Ni_upper")
df$group<-c(rep("Estimate",2*year.range),rep("2.5%CI",2*year.range),rep("97.5%CI",2*year.range))
df$group2<-rep(c(rep("long",year.range),rep("short",year.range)),3)  
df

df_long<-subset(df,group2=="long")
row.names(df_long)<-seq(1:nrow(df_long))
df_short<-subset(df,group2=="short")
row.names(df_short)<-seq(1:nrow(df_short))

ggplot(data=df_long)+
  geom_ribbon(mapping=aes(x=years,ymin=Ni_lower,ymax=Ni_upper,fill=group), alpha=0.2)+
  geom_line(mapping=aes(x=years,y=Ni,group=group,color=group),size=1)+ 
  #theme(legend.title=element_blank())+ 
  guides(fill = guide_legend(reverse=TRUE),color=guide_legend(reverse=TRUE))+
  labs(x = "Years", y="Population size", fill="Current population",color="Current population")+
  xlim(1,year.range)+
  ylim(0,100000)+
  annotate("text", x= 5, y = 100000, label = ("Long term trend = -1.78 (-3.56, 0.174)"),size=5)

ggplot(data=df_short)+
  geom_ribbon(mapping=aes(x=years,ymin=Ni_lower,ymax=Ni_upper,fill=group), alpha=0.2)+
  geom_line(mapping=aes(x=years,y=Ni,group=group,color=group),size=1)+ 
  #theme(legend.title=element_blank())+ 
  guides(fill = guide_legend(reverse=TRUE),color=guide_legend(reverse=TRUE))+
  labs(x = "Years", y="Population size", fill="Current population",color="Current population")+
  #xlim(1,year.range)+
  #ylim(0,100000)+
  annotate("text", x= 5, y = 1000000, label = ("Short term trend = 5.46 (-3.31, 16)"),size=5)

# # find out what year the current estimates most relate to, based on year of avian data:
load("D:/Avian data processed/BAM_data_package_August2019.RData")
yearsdata<-PKEYcombo$YEAR[match(datcombo$PKEY,PKEYcombo$PKEY)]
yearsdata<-as.data.frame(t(table(yearsdata)))[,2:3]
names(yearsdata)<-c("Year","Freq")
yearsdata$Year<-as.numeric(as.character(yearsdata$Year))
str(yearsdata)

median(PKEYcombo$YEAR[match(datcombo$PKEY,PKEYcombo$PKEY)])
mean(PKEYcombo$YEAR[match(datcombo$PKEY,PKEYcombo$PKEY)])
plot(table(PKEYcombo$YEAR[match(datcombo$PKEY,PKEYcombo$PKEY)]),ylab="Frequency")

# 2006 is median, but use 2007 to include in short term trend
BBS_AI <- read.csv("D:/CHID regional NS BRT/BBS_annual_indices_CAWA_NS.csv")

AI2007<-BBS_AI[which(BBS_AI$year==2007&BBS_AI$timeFrame=="Long-term"),7:9]

df_long

newdf<-data.frame(year=c(1970:2040),annualIndex=NA,lower=NA,upper=NA,psize=NA,psizeLow=NA,psizeUpp=NA)
newdf[1:48,2:4]<-BBS_AI[1:48,7:9]

for(i in c(1:48)){
  newdf[i,5]<-newdf[i,2]*popsizeCIs[1]/newdf[38,2]
  newdf[i,6]<-newdf[i,2]*popsizeCIs[2]/newdf[38,2]
  newdf[i,7]<-newdf[i,2]*popsizeCIs[3]/newdf[38,2]
}

newdf

for(i in 49:71){
  newdf[i,5]<-newdf[48,5]*(1+trends[1]/100)^(i-48)
  newdf[i,6]<-newdf[48,6]*(1+trends[1]/100)^(i-48)
  newdf[i,7]<-newdf[48,7]*(1+trends[1]/100)^(i-48)

  newdf[i,2]<-round(AI2007[1]*newdf[i,5]/newdf[38,5],3)
  newdf[i,3]<-round(AI2007[1]*newdf[i,6]/newdf[38,5],3)
  newdf[i,4]<-round(AI2007[1]*newdf[i,7]/newdf[38,5],3)


}

ggplot(data=newdf[1:71,])+
  geom_ribbon(aes(x=year,ymin=upper,ymax=lower),alpha=0.3)+
  geom_line(aes(x=year,y=annualIndex))+
  geom_vline(xintercept = 2007,  color = "blue", size=0.5)+
  geom_vline(xintercept = 2018,  color = "red", size=0.5)


ggplot(data=newdf[1:71,])+
  geom_ribbon(aes(x=year,ymin=psizeLow,ymax=psizeUpp),alpha=0.3)+
  geom_line(aes(x=year,y=psize))+
  geom_vline(xintercept = 2007,  color = "blue", size=0.5)+
  geom_vline(xintercept = 2018,  color = "red", size=0.5)


newdf[72:117,]<-NA

newdf$group<-c(rep("Trend estimate",71),rep("2.5%CI",23),rep("97.5%CI",23))

newdf[72:117,1]<-rep(2018:2040,times=2)

for(i in 72:94){
  newdf[i,5]<-newdf[48,5]*(1+trends[2]/100)^(i-68)
  newdf[i,6]<-newdf[48,6]*(1+trends[2]/100)^(i-68)
  newdf[i,7]<-newdf[48,7]*(1+trends[2]/100)^(i-68)
  
  newdf[i,2]<-round(AI2007[1]*newdf[i,5]/newdf[38,5],3)
  newdf[i,3]<-round(AI2007[1]*newdf[i,6]/newdf[38,5],3)
  newdf[i,4]<-round(AI2007[1]*newdf[i,7]/newdf[38,5],3)
}
newdf

for(i in 95:117){
  newdf[i,5]<-newdf[48,5]*(1+trends[3]/100)^(i-68)
  newdf[i,6]<-newdf[48,6]*(1+trends[3]/100)^(i-68)
  newdf[i,7]<-newdf[48,7]*(1+trends[3]/100)^(i-68)
  
  newdf[i,2]<-round(AI2007[1]*newdf[i,5]/newdf[38,5],3)
  newdf[i,3]<-round(AI2007[1]*newdf[i,6]/newdf[38,5],3)
  newdf[i,4]<-round(AI2007[1]*newdf[i,7]/newdf[38,5],3)
}
newdf


ggplot(data=newdf)+
  geom_ribbon(aes(x=year,ymin=upper,ymax=lower,fill=group),alpha=0.2)+
  geom_line(aes(x=year,y=annualIndex,color=group))+
  geom_vline(xintercept = 2007,  size=0.5)+
  geom_segment(x=2018,y=0.7,xend=2037,yend=0.7,size=0.5,arrow = arrow(angle=90,ends="both",length=unit(0.08,"inches")))+
  guides(fill = guide_legend(reverse=TRUE),color=guide_legend(reverse=TRUE))+
  labs(x = "Years", y="Annual Index", fill="",color="")+
  annotate("text", x= 2006, y = 2,angle=90, label = ("BRT model reference year"),size=4)+
  annotate("text", x= 2028, y = 0.8, label = ("Predictions with BBS trend (2017)"),size=4)

library(scales)

ggplot(data=newdf)+
  geom_ribbon(aes(x=year,ymin=psizeLow,ymax=psizeUpp,fill=group),alpha=0.2)+
  geom_line(aes(x=year,y=psize,color=group))+
  geom_vline(xintercept = 2007,  size=0.5)+
  geom_segment(x=2018,y=175000,xend=2040,yend=175000,size=0.5,arrow = arrow(angle=90,ends="both",length=unit(0.08,"inches")))+
  guides(fill = guide_legend(reverse=TRUE),color=guide_legend(reverse=TRUE))+
  labs(x = "Years", y="Population size", fill="",color="")+
  scale_y_continuous(labels = comma)+
  annotate("text", x= 2006, y = 500000,angle=90, label = ("BRT model reference year (2007)"),size=4)+
  annotate("text", x= 2028, y = 200000, label = ("Predictions with BBS trend (2017)"),size=4)+
  geom_rug(data=yearsdata, aes(x=Year,y=Freq),sides="b", alpha=0.5, position="jitter")

# Population size trajectory and future predictions with BBS long-term trend

newdf2<-subset(newdf,group=="Trend estimate")
newdf2[49:71,"psizeLow"]<-newdf[72:94,"psizeLow"]
newdf2[49:71,"psizeUpp"]<-newdf[95:117,"psizeUpp"]

jitter <- position_jitter(width = 0.05, height = 1)

p1<-ggplot(data=newdf2)+
  geom_ribbon(aes(x=year,ymin=psizeLow,ymax=psizeUpp),fill="blue",alpha=0.2)+
  geom_line(aes(x=year,y=psize),color="blue")+
  geom_vline(xintercept = 2007,  size=0.5,linetype="dashed")+
  geom_errorbarh(aes(xmax=2040,xmin=2018,y=126000,height=20000))+
  guides(fill = guide_legend(reverse=TRUE),color=guide_legend(reverse=TRUE))+
  labs(x = "Years", y="Population size (males)", fill="",color="")+
  scale_y_continuous(labels = comma)+
  annotate("text", x= 2008, y = 400000,hjust=0, label = paste0("BRT model reference year (2007) \n", paste0("Estimate (90% CI) = ",round(popsize,0)," (",round(popsize.low,0)," - ",round(popsize.upp,0),")")),size=4)+
  annotate("text", x= 2029, y = 155000, label = ("Predictions with BBS long-term trend (2017): \n -1.78 (-3.560, 0.124)"),size=4)+
  geom_rug(data=yearsdata, aes(x=Year,y=Freq),length=unit(yearsdata$Freq/50000,"npc"),sides="b", alpha=0.5, position=jitter)

p1


# Population size trajectory and future predictions with BBS short-term trend (10 years)
newdf3<-data.frame(year=c(2018:2040),psize=NA,psizeLow=NA,psizeUpp=NA)
newdf3<-rbind(newdf3,newdf3,newdf3)
newdf3$group<-c(rep("Trend estimate",23),rep("2.5%CI",23),rep("97.5%CI",23))

for(i in 1:23){
  newdf3[i,2]<-newdf[48,5]*(1+trends[4]/100)^(newdf3$year[i]-2017)
  newdf3[i,3]<-newdf[48,6]*(1+trends[4]/100)^(newdf3$year[i]-2017)
  newdf3[i,4]<-newdf[48,7]*(1+trends[4]/100)^(newdf3$year[i]-2017)
  }

for(i in 24:46){
  newdf3[i,2]<-newdf[48,5]*(1+trends[5]/100)^(newdf3$year[i]-2017)
  newdf3[i,3]<-newdf[48,6]*(1+trends[5]/100)^(newdf3$year[i]-2017)
  newdf3[i,4]<-newdf[48,7]*(1+trends[5]/100)^(newdf3$year[i]-2017)
}

for(i in 47:69){
  newdf3[i,2]<-newdf[48,5]*(1+trends[6]/100)^(newdf3$year[i]-2017)
  newdf3[i,3]<-newdf[48,6]*(1+trends[6]/100)^(newdf3$year[i]-2017)
  newdf3[i,4]<-newdf[48,7]*(1+trends[6]/100)^(newdf3$year[i]-2017)
}
newdf3

newdf4<-newdf2
newdf4[49:71,"psize"]<-newdf3[1:23,"psize"]
newdf4[49:71,"psizeLow"]<-newdf3[24:46,"psizeLow"]
newdf4[49:71,"psizeUpp"]<-newdf3[47:69,"psizeUpp"]

p2<- ggplot(data=newdf4)+
  geom_ribbon(aes(x=year,ymin=psizeLow,ymax=psizeUpp),fill="darkgreen",alpha=0.2)+
  geom_line(aes(x=year,y=psize),color="darkgreen")+
  geom_vline(xintercept = 2007,  size=0.5,linetype="dashed")+
  geom_errorbarh(aes(xmax=2040,xmin=2018,y=210000,height=20000))+
  guides(fill = guide_legend(reverse=TRUE),color=guide_legend(reverse=TRUE))+
  labs(x = "Years", y="Population size (males)", fill="",color="")+
  scale_y_continuous(labels = comma)+
  annotate("text", x= 2008, y = 650000,hjust=0, label = paste0("BRT model reference year (2007) \n", paste0("Estimate (90% CI) = ",round(popsize,0)," (",round(popsize.low,0)," - ",round(popsize.upp,0),")")),size=4)+
  annotate("text", x= 2029, y = 240000, label = ("Predictions with BBS short-term trend (2017): \n 5.460 (-3.310, 16.000)"),size=4)+
  geom_rug(data=yearsdata, aes(x=Year,y=Freq),length=unit(yearsdata$Freq/50000,"npc"),sides="b", alpha=0.5, position=jitter)+
  coord_cartesian(ylim = c(0, 675000)) 
p2  

library(gridExtra)
grid.arrange(p1,p2,nrow=2)

# Compare future pop.size based on trend with that from simulations in LANDIS ####

getLandisPreds2040<-function(scenario){
  filelist<-list.files("D:/CHID regional NS BRT/Landscape simulation results/Prediction raster/",pattern=scenario)
  preds<-rep(NA,5)
  for(i in 1:length(filelist)){
    load(paste0("D:/CHID regional NS BRT/Landscape simulation results/Prediction raster/",filelist[i]))
    landispreds2040<-PredictBirdsp.l[[3]] # 3 is for 2040
    
    landispreds2040<-landispreds2040*6.25  # Multiplying density values by 6.25 (hectars in a 250-250 pixel) to get abundance for each cell and overall population size
    preds[i]<-cellStats(landispreds2040,stat=sum,na.rm=T)
  }
  preds.out<-data.frame(mean=mean(preds),sd=sd(preds))
  
  return(preds.out)
  
}

scenarios<-getLandisPreds2040("baseline_BudwormBaselineFire_")
scenarios[2,]<-getLandisPreds2040("baseline_BudwormBaselineFireEBFMHarvest_")
scenarios[3,]<-getLandisPreds2040("baseline_BudwormBaselineFireHighCPRSHarvest_")
scenarios[4,]<-getLandisPreds2040("baseline_BudwormBaselineFireHighPartialHarvest_")
scenarios[5,]<-getLandisPreds2040("baseline_BudwormBaselineFireLowCPRSHarvest_")
scenarios[6,]<-getLandisPreds2040("baseline_BudwormBaselineFireLowPartialHarvest_")

scenarios[7,]<-getLandisPreds2040("RCP26_GrowthBudwormProjectedFire_")
scenarios[8,]<-getLandisPreds2040("RCP26_GrowthBudwormProjectedFireEBFMHarvest_")
scenarios[9,]<-getLandisPreds2040("RCP26_GrowthBudwormProjectedFireHighCPRSHarvest_")
scenarios[10,]<-getLandisPreds2040("RCP26_GrowthBudwormProjectedFireHighPartialHarvest_")
scenarios[11,]<-getLandisPreds2040("RCP26_GrowthBudwormProjectedFireLowCPRSHarvest_")
scenarios[12,]<-getLandisPreds2040("RCP26_GrowthBudwormProjectedFireLowPartialHarvest_")

scenarios[13,]<-getLandisPreds2040("RCP45_GrowthBudwormProjectedFire_")
scenarios[14,]<-getLandisPreds2040("RCP45_GrowthBudwormProjectedFireEBFMHarvest_")
scenarios[15,]<-getLandisPreds2040("RCP45_GrowthBudwormProjectedFireHighCPRSHarvest_")
scenarios[16,]<-getLandisPreds2040("RCP45_GrowthBudwormProjectedFireHighPartialHarvest_")
scenarios[17,]<-getLandisPreds2040("RCP45_GrowthBudwormProjectedFireLowCPRSHarvest_")
scenarios[18,]<-getLandisPreds2040("RCP45_GrowthBudwormProjectedFireLowPartialHarvest_")

scenarios[19,]<-getLandisPreds2040("RCP85_GrowthBudwormProjectedFire_")
scenarios[20,]<-getLandisPreds2040("RCP85_GrowthBudwormProjectedFireEBFMHarvest_")
scenarios[21,]<-getLandisPreds2040("RCP85_GrowthBudwormProjectedFireHighCPRSHarvest_")
scenarios[22,]<-getLandisPreds2040("RCP85_GrowthBudwormProjectedFireHighPartialHarvest_")
scenarios[23,]<-getLandisPreds2040("RCP85_GrowthBudwormProjectedFireLowCPRSHarvest_")
scenarios[24,]<-getLandisPreds2040("RCP85_GrowthBudwormProjectedFireLowPartialHarvest_")

# adding columns with scenario names
filelist<-list.files("D:/CHID regional NS BRT/Landscape simulation results/Prediction raster/")

scenario.names<-rep(NA,24)
for(i in 1:24 ){
  splits<-unlist(strsplit(filelist[seq(1,length(filelist),by=5)[i]], split="_"))
  name.sce<-paste(splits[3],splits[4],sep="_")
  scenario.names[i]<-name.sce
}



scenarios$scenario<-scenario.names
scenarios$Climate<-rep(c("Baseline","RCP26","RCP45","RCP85"),each=6)


scenarios

p_inset<-
  ggplot(data=newdf2)+
  geom_line(aes(x=year,y=psize),color="blue",size=1)+
  guides(fill = guide_legend(reverse=TRUE),color=guide_legend(reverse=TRUE))+
  scale_y_continuous(labels = comma)+
  geom_pointrange(data=scenarios,mapping=aes(x=2040,y=mean,ymin=mean-1.96*sd,ymax=mean+1.96*sd,color=Climate),position=position_dodge2(width=0.2),alpha=0.8)+ 
  coord_cartesian(xlim=c(2039,2042),ylim=c(39000,58000) )+
  scale_x_continuous(labels=c("","2040","",""))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.75, 0.5)
        )

p_inset

p1+
  geom_pointrange(data=scenarios,mapping=aes(x=2040,y=mean,ymin=mean-1.96*sd,ymax=mean+1.96*sd))+ theme(legend.position = "none")+
  annotation_custom(grob = ggplotGrob(p_inset), xmin=2020, xmax=2040, ymin=200000,ymax=350000)


p_inset2<-
  ggplot(data=newdf4)+
  geom_line(aes(x=year,y=psize),color="darkgreen")+
  geom_vline(xintercept = 2007,  size=0.5,linetype="dashed")+
  geom_segment(x=2018,y=275000,xend=2040,yend=275000,size=0.5,arrow = arrow(angle=90,ends="both",length=unit(0.08,"inches")))+
  guides(fill = guide_legend(reverse=TRUE),color=guide_legend(reverse=TRUE))+
  labs(x = "Years", y="Population size", fill="",color="")+
  scale_y_continuous(labels = comma)+
  geom_pointrange(data=scenarios,mapping=aes(x=2040,y=mean,ymin=mean-1.96*sd,ymax=mean+1.96*sd,color=Climate),position=position_dodge2(width=0.2),alpha=0.8)+
  coord_cartesian(xlim=c(2039,2042),ylim=c(39000,200000))+
  scale_x_continuous(labels=c("","2040","",""))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.75, 0.5)) 
#p_inset2  
p2+
  geom_pointrange(data=scenarios,mapping=aes(x=2040,y=mean,ymin=mean-1.96*sd,ymax=mean+1.96*sd))+ theme(legend.position = "none")+
  annotation_custom(grob = ggplotGrob(p_inset2), xmin=2020, xmax=2040, ymin=400000,ymax=600000)


# Future population size probability ####
# long trend #
future.pop<-newdf2[which(newdf2$year=="2040"),"psize"] ###

future.pop.upp<-newdf2[which(newdf2$year=="2040"),"psizeUpp"]
future.pop.low<-newdf2[which(newdf2$year=="2040"),"psizeLow"]


GammaParamsFuture<-findGamma1(scale.init=c(2000,80000),shape.init = c(0.25,15),upp=future.pop.upp,low=future.pop.low)
GammaParamsFuture


futurepoprange<-seq(1000,200000,by=1000)
id.dist.future<-pgamma(futurepoprange,shape=GammaParamsFuture$shape,scale=GammaParamsFuture$scale,lower.tail = F)
id.df.future<-data.frame("Pop.range"=futurepoprange,"Probability"=id.dist.future)
yfuturepred<-pgamma(future.pop,shape=GammaParamsFuture$shape,scale=GammaParamsFuture$scale,lower.tail = F)

ggplot(id.df.future, aes(x = Pop.range, y = Probability)) + 
  geom_point()+ 
  geom_errorbarh(aes(xmax=future.pop.upp,xmin=future.pop.low,y=yfuturepred,height=0.05))+
  geom_point(x = future.pop,y=yfuturepred,col="red")+  
  labs(x="N",y="Probability (population size >= N)")+
  annotate("text", x= 55000, y = 0.625,angle=0, label = paste0("Model prediction (90% CI) = ",round(future.pop,0)," (",round(future.pop.low,0)," - ",round(future.pop.upp,0),")"),size=4, hjust=0)


# If using log-normal distribution
# mu<-log(future.pop)
# 
# phi<-0.73737
# qlnorm(p=c(0.025,0.975),meanlog = mu,sdlog=phi)
# 
# 
# futurepoprange<-seq(1000,200000,by=1000)
# id.dist.future<-plnorm(futurepoprange,meanlog=mu,sdlog=phi,lower.tail = F)
# id.df.future<-data.frame("Pop.range"=futurepoprange,"Probability"=id.dist.future)
# 
# ggplot(id.df.future, aes(x = Pop.range, y = Probability)) + geom_point()+ geom_vline(xintercept = future.pop,  color = "blue", size=0.8)+
#   labs(x="Population size")


# short trend #
# future.pop2<-newdf4[which(newdf4$year=="2040"),"psize"] ###
# 
# upp2<-newdf4[which(newdf4$year=="2040"),"psizeUpp"]
# low2<-newdf4[which(newdf4$year=="2040"),"psizeLow"]
# 
# mat2<-expand.grid(scale=runif(1000,min=1000,max=100000),shape=runif(1000,min=0.2,max=30))
# mat2[,3:4]<-t(apply(mat2,1,FUN = function(x){qgamma(p=c(0.025,0.975),shape=x[2],scale=x[1])}))
# mat2$distances<-pointDistance(p1=mat2[,3:4],p2=c(low2,upp2),lonlat = F)
# mat2[which.min(mat2$distances),]
# 
# low2
# upp2
# qgamma(p=c(0.025,0.975),shape = 15,scale=220000)
# 
# futurepoprange<-seq(1000,200000,by=1000)
# id.dist.future<-pgamma(futurepoprange,shape=2.705717,scale=18677.86,lower.tail = F)
# id.df.future<-data.frame("Pop.range"=futurepoprange,"Probability"=id.dist.future)



# Future conditions population size - LANDIS simulations ####

#1. Predictions for 2100 
getLandisPreds2100<-function(scenario){
  filelist<-list.files("D:/CHID regional NS BRT/Landscape simulation results/Prediction raster/",pattern=scenario)
  preds<-rep(NA,5)
  for(i in 1:length(filelist)){
    load(paste0("D:/CHID regional NS BRT/Landscape simulation results/Prediction raster/",filelist[i]))
    landispreds2040<-PredictBirdsp.l[[6]] # 6 is for 2100
    
    landispreds2040<-landispreds2040*6.25  # Multiplying density values by 6.25 (hectars in a 250-250 pixel) to get abundance for each cell and overall population size
    preds[i]<-cellStats(landispreds2040,stat=sum,na.rm=T)
  }
  preds.out<-data.frame(mean=mean(preds),sd=sd(preds))
  
  return(preds.out)
  }

scenarios2100<-getLandisPreds2100("baseline_BudwormBaselineFire_")
scenarios2100[2,]<-getLandisPreds2100("baseline_BudwormBaselineFireEBFMHarvest_")
scenarios2100[3,]<-getLandisPreds2100("baseline_BudwormBaselineFireHighCPRSHarvest_")
scenarios2100[4,]<-getLandisPreds2100("baseline_BudwormBaselineFireHighPartialHarvest_")
scenarios2100[5,]<-getLandisPreds2100("baseline_BudwormBaselineFireLowCPRSHarvest_")
scenarios2100[6,]<-getLandisPreds2100("baseline_BudwormBaselineFireLowPartialHarvest_")

scenarios2100[7,]<-getLandisPreds2100("RCP26_GrowthBudwormProjectedFire_")
scenarios2100[8,]<-getLandisPreds2100("RCP26_GrowthBudwormProjectedFireEBFMHarvest_")
scenarios2100[9,]<-getLandisPreds2100("RCP26_GrowthBudwormProjectedFireHighCPRSHarvest_")
scenarios2100[10,]<-getLandisPreds2100("RCP26_GrowthBudwormProjectedFireHighPartialHarvest_")
scenarios2100[11,]<-getLandisPreds2100("RCP26_GrowthBudwormProjectedFireLowCPRSHarvest_")
scenarios2100[12,]<-getLandisPreds2100("RCP26_GrowthBudwormProjectedFireLowPartialHarvest_")

scenarios2100[13,]<-getLandisPreds2100("RCP45_GrowthBudwormProjectedFire_")
scenarios2100[14,]<-getLandisPreds2100("RCP45_GrowthBudwormProjectedFireEBFMHarvest_")
scenarios2100[15,]<-getLandisPreds2100("RCP45_GrowthBudwormProjectedFireHighCPRSHarvest_")
scenarios2100[16,]<-getLandisPreds2100("RCP45_GrowthBudwormProjectedFireHighPartialHarvest_")
scenarios2100[17,]<-getLandisPreds2100("RCP45_GrowthBudwormProjectedFireLowCPRSHarvest_")
scenarios2100[18,]<-getLandisPreds2100("RCP45_GrowthBudwormProjectedFireLowPartialHarvest_")

scenarios2100[19,]<-getLandisPreds2100("RCP85_GrowthBudwormProjectedFire_")
scenarios2100[20,]<-getLandisPreds2100("RCP85_GrowthBudwormProjectedFireEBFMHarvest_")
scenarios2100[21,]<-getLandisPreds2100("RCP85_GrowthBudwormProjectedFireHighCPRSHarvest_")
scenarios2100[22,]<-getLandisPreds2100("RCP85_GrowthBudwormProjectedFireHighPartialHarvest_")
scenarios2100[23,]<-getLandisPreds2100("RCP85_GrowthBudwormProjectedFireLowCPRSHarvest_")
scenarios2100[24,]<-getLandisPreds2100("RCP85_GrowthBudwormProjectedFireLowPartialHarvest_")

scenarios2100$scenario<-scenario.names
scenarios2100$Climate<-rep(c("Baseline","RCP26","RCP45","RCP85"),each=6)

scenarios2100
#CONTINUE - LANDIS BOOTSTRAP PREDICTIONS ####
#2. Load bootstrap predictions tables
bootpreds<-readRDS("D:/CHID regional NS BRT/Landscape simulation results/bootstrap_predictions2100.R")
str(bootpreds)
bootpreds$treatment<-as.factor(bootpreds$treatment)
bootpreds$Year<-as.factor(bootpreds$Year)
bootpreds$FireTreatment<-as.factor(bootpreds$FireTreatment)
bootpreds$HarvestTreatment<-as.factor(bootpreds$HarvestTreatment)

bootpreds$scenarioname<-rep(rep(scenario.names,each=250),5)


#3. Quantify uncertainty from bootstrap predictions, for each scenario. This is achieved by aproximating the parameters of a gamma distribution whose quantiles fit those of the bootstrap distribution
bootpreds2100<-subset(bootpreds,Year=="2100")
summary(bootpreds2100)
levels(bootpreds$treatment)

# each of the objects below, representing the different scenarios, has 250 (bootstrao samples) x 5 (simulation replicates) = 1250 population size estimates
baseline_BudwormBaselineFire_2100 <- bootpreds2100[which(bootpreds2100$treatment=="BudwormBaselineFire"),] 
baseline_BudwormBaselineFireEBFMHarvest_2100 <- bootpreds2100[which(bootpreds2100$treatment=="BudwormBaselineFireEBFMHarvest"),]
baseline_BudwormBaselineFireHighCPRSHarvest_2100 <- bootpreds2100[which(bootpreds2100$treatment=="BudwormBaselineFireHighCPRSHarvest"),]
baseline_BudwormBaselineFireHighPartialHarvest_2100 <- bootpreds2100[which(bootpreds2100$treatment=="BudwormBaselineFireHighPartialHarvest"),]
baseline_BudwormBaselineFireLowCPRSHarvest_2100 <- bootpreds2100[which(bootpreds2100$treatment=="BudwormBaselineFireLowCPRSHarvest"),]
baseline_BudwormBaselineFireLowPartialHarvest_2100 <- bootpreds2100[which(bootpreds2100$treatment=="BudwormBaselineFireLowPartialHarvest"),]

RCP26_GrowthBudwormProjectedFire_2100 <- bootpreds2100[which(bootpreds2100$treatment=="GrowthBudwormProjectedFire"&bootpreds2100$scenario=="RCP26"),]
RCP26_GrowthBudwormProjectedFireEBFMHarvest_2100 <- bootpreds2100[which(bootpreds2100$treatment=="GrowthBudwormProjectedFireEBFMHarvest"&bootpreds2100$scenario=="RCP26"),]
RCP26_GrowthBudwormProjectedFireHighCPRSHarvest_2100 <- bootpreds2100[which(bootpreds2100$treatment=="GrowthBudwormProjectedFireHighCPRSHarvest"&bootpreds2100$scenario=="RCP26"),]
RCP26_GrowthBudwormProjectedFireHighPartialHarvest_2100 <- bootpreds2100[which(bootpreds2100$treatment=="GrowthBudwormProjectedFireHighPartialHarvest"&bootpreds2100$scenario=="RCP26"),]
RCP26_GrowthBudwormProjectedFireLowCPRSHarvest_2100 <- bootpreds2100[which(bootpreds2100$treatment=="GrowthBudwormProjectedFireLowCPRSHarvest"&bootpreds2100$scenario=="RCP26"),]
RCP26_GrowthBudwormProjectedFireLowPartialHarvest_2100 <- bootpreds2100[which(bootpreds2100$treatment=="GrowthBudwormProjectedFireLowPartialHarvest"&bootpreds2100$scenario=="RCP26"),]

RCP45_GrowthBudwormProjectedFire_2100 <- bootpreds2100[which(bootpreds2100$treatment=="GrowthBudwormProjectedFire"&bootpreds2100$scenario=="RCP45"),]
RCP45_GrowthBudwormProjectedFireEBFMHarvest_2100 <- bootpreds2100[which(bootpreds2100$treatment=="GrowthBudwormProjectedFireEBFMHarvest"&bootpreds2100$scenario=="RCP45"),]
RCP45_GrowthBudwormProjectedFireHighCPRSHarvest_2100 <- bootpreds2100[which(bootpreds2100$treatment=="GrowthBudwormProjectedFireHighCPRSHarvest"&bootpreds2100$scenario=="RCP45"),]
RCP45_GrowthBudwormProjectedFireHighPartialHarvest_2100 <- bootpreds2100[which(bootpreds2100$treatment=="GrowthBudwormProjectedFireHighPartialHarvest"&bootpreds2100$scenario=="RCP45"),]
RCP45_GrowthBudwormProjectedFireLowCPRSHarvest_2100 <- bootpreds2100[which(bootpreds2100$treatment=="GrowthBudwormProjectedFireLowCPRSHarvest"&bootpreds2100$scenario=="RCP45"),]
RCP45_GrowthBudwormProjectedFireLowPartialHarvest_2100 <- bootpreds2100[which(bootpreds2100$treatment=="GrowthBudwormProjectedFireLowPartialHarvest"&bootpreds2100$scenario=="RCP45"),]

RCP85_GrowthBudwormProjectedFire_2100 <- bootpreds2100[which(bootpreds2100$treatment=="GrowthBudwormProjectedFire"&bootpreds2100$scenario=="RCP85"),]
RCP85_GrowthBudwormProjectedFireEBFMHarvest_2100 <- bootpreds2100[which(bootpreds2100$treatment=="GrowthBudwormProjectedFireEBFMHarvest"&bootpreds2100$scenario=="RCP85"),]
RCP85_GrowthBudwormProjectedFireHighCPRSHarvest_2100 <- bootpreds2100[which(bootpreds2100$treatment=="GrowthBudwormProjectedFireHighCPRSHarvest"&bootpreds2100$scenario=="RCP85"),]
RCP85_GrowthBudwormProjectedFireHighPartialHarvest_2100 <- bootpreds2100[which(bootpreds2100$treatment=="GrowthBudwormProjectedFireHighPartialHarvest"&bootpreds2100$scenario=="RCP85"),]
RCP85_GrowthBudwormProjectedFireLowCPRSHarvest_2100 <- bootpreds2100[which(bootpreds2100$treatment=="GrowthBudwormProjectedFireLowCPRSHarvest"&bootpreds2100$scenario=="RCP85"),]
RCP85_GrowthBudwormProjectedFireLowPartialHarvest_2100 <- bootpreds2100[which(bootpreds2100$treatment=="GrowthBudwormProjectedFireLowPartialHarvest"&bootpreds2100$scenario=="RCP85"),]

boot_preds_2100_list<-list(baseline_BudwormBaselineFire_2100,
                           baseline_BudwormBaselineFireEBFMHarvest_2100,
                           baseline_BudwormBaselineFireHighCPRSHarvest_2100,
                           baseline_BudwormBaselineFireHighPartialHarvest_2100,
                           baseline_BudwormBaselineFireLowCPRSHarvest_2100,
                           baseline_BudwormBaselineFireLowPartialHarvest_2100,
                           RCP26_GrowthBudwormProjectedFire_2100,
                           RCP26_GrowthBudwormProjectedFireEBFMHarvest_2100,
                           RCP26_GrowthBudwormProjectedFireHighCPRSHarvest_2100,
                           RCP26_GrowthBudwormProjectedFireHighPartialHarvest_2100,
                           RCP26_GrowthBudwormProjectedFireLowCPRSHarvest_2100,
                           RCP26_GrowthBudwormProjectedFireLowPartialHarvest_2100,
                           RCP45_GrowthBudwormProjectedFire_2100,
                           RCP45_GrowthBudwormProjectedFireEBFMHarvest_2100,
                           RCP45_GrowthBudwormProjectedFireHighCPRSHarvest_2100,
                           RCP45_GrowthBudwormProjectedFireHighPartialHarvest_2100,
                           RCP45_GrowthBudwormProjectedFireLowCPRSHarvest_2100,
                           RCP45_GrowthBudwormProjectedFireLowPartialHarvest_2100,
                           RCP85_GrowthBudwormProjectedFire_2100,
                           RCP85_GrowthBudwormProjectedFireEBFMHarvest_2100,
                           RCP85_GrowthBudwormProjectedFireHighCPRSHarvest_2100,
                           RCP85_GrowthBudwormProjectedFireHighPartialHarvest_2100,
                           RCP85_GrowthBudwormProjectedFireLowCPRSHarvest_2100,
                           RCP85_GrowthBudwormProjectedFireLowPartialHarvest_2100) 

probsCI<-c(0.05,0.95) # can change these probabilities for more or less conservative intervals

boot_preds_2100_CIs<-matrix(unlist(lapply(boot_preds_2100_list,FUN=function(x){quantile(x$Popsize, probs=probsCI)})),24,2,byrow = T)
boot_preds_2100_CIs

scenarios2100$bootLower<-boot_preds_2100_CIs[,1]
scenarios2100$bootUpper<-boot_preds_2100_CIs[,2]
scenarios2100 

findGamma2<-function(x,scale.init,shape.init){
 low<-as.numeric(x[5])
 upp<-as.numeric(x[6])
 mat<-expand.grid(scale=runif(1000,min=scale.init[1],max=scale.init[2]),shape=runif(1000,min=shape.init[1],max=shape.init[2]))
 mat[,3:4]<-t(apply(mat,1,FUN = function(x){qgamma(p=c(0.025,0.975),shape=x[2],scale=x[1])}))
 mat$distances<-pointDistance(p1=mat[,3:4],p2=c(low,upp),lonlat = F)
 out1<-mat[which.min(mat$distances),]
 
 scale.init2<-unlist(c(out1[1]-(out1[1]*0.25),out1[1]+(out1[1]*0.25)))
 shape.init2<-unlist(c(out1[2]-(out1[2]*0.25),out1[2]+(out1[2]*0.25)))
 
 mat2<-expand.grid(scale=runif(1000,min=scale.init2[1],max=scale.init2[2]),shape=runif(1000,min=shape.init2[2],max=shape.init2[2]))
 mat2[,3:4]<-t(apply(mat2,1,FUN = function(x){qgamma(p=c(0.025,0.975),shape=x[2],scale=x[1])}))
 mat2$distances<-pointDistance(p1=mat2[,3:4],p2=c(low,upp),lonlat = F)
 
 mat3<-rbind(mat,mat2)
 mat3[which.min(mat3$distances),]
 }

parametersGamma2100<-apply(scenarios2100,MARGIN = 1,FUN=findGamma2,scale.init=c(20000,80000),shape.init = c(0.25,15))
parametersGamma2100<-as.data.frame(matrix(unlist(parametersGamma2100),ncol=5,nrow=24,byrow=T, dimnames=list(1:24,names(parametersGamma2100[[1]]))))
parametersGamma2100$scenario<-scenario.names
parametersGamma2100
names(parametersGamma2100)[c(3,4)]<-c("lowerCL","upperCL")

# Scenario-specific future pop size probability plots ####

pop.prob.2100<-as.list(scenarios2100$scenario)
names(pop.prob.2100)<-scenarios2100$scenario

id.df.future.2100<-data.frame("Pop.range"=NA,"Probability"=NA,"scenario"=NA)
for(i in 1:length(pop.prob.2100)){
  futurepoprange<-seq(1000,350000,by=1000)
  id.dist.future.2100<-pgamma(futurepoprange,shape=parametersGamma2100$shape[i],scale=parametersGamma2100$scale[i],lower.tail = F)
  
  id.df.future.2100[(1:length(futurepoprange))+length(futurepoprange)*(i-1),]<-cbind(futurepoprange,id.dist.future.2100,names(pop.prob.2100)[i])
}

id.df.future.2100$Pop.range<-as.numeric(id.df.future.2100$Pop.range)
id.df.future.2100$Probability<-as.numeric(id.df.future.2100$Probability)
str(id.df.future.2100)

pop.prob.2100.plot<-ggplot(id.df.future.2100, aes(x = Pop.range, y = Probability,color=scenario)) +
  geom_point()+
  geom_vline(mapping=aes(xintercept = mean,color=scenario),data=scenarios2100, size=0.8,alpha=0.4)+
  labs(x="Population size")+
  annotate("text", x= 55000, y = 0.93,angle=0, label = paste0("Future population size predictions \n (mean = ",round(mean(scenarios2100$mean),0),")"),size=4, hjust=0)+
  #annotate("text", x= 150000, y = 0.93,angle=0, label = paste0("Scenario: ",scenarios2100$scenario[i]),size=4)+
  ylim(0,1)+
  scale_color_discrete(guide=guide_legend(ncol=1))+
  scale_x_continuous(labels = scales::comma)+
  labs(x="N",y="Probability (population size >= N)")
pop.prob.2100.plot

# population trajectories for each scenario
getLandisTrajectory<-function(scenario){
  filelist<-list.files("D:/CHID regional NS BRT/Landscape simulation results/Prediction raster/",pattern=scenario)
  preds<-matrix(NA,11,5)
  for(i in 1:length(filelist)){
    load(paste0("D:/CHID regional NS BRT/Landscape simulation results/Prediction raster/",filelist[i]))
    preds[,i]<-unlist(lapply(PredictBirdsp.l,FUN=function(x){cellStats(x*6.25,stat=sum,na.rm=T)}))
    
  }
  preds.out<-data.frame(mean=apply(preds,1,mean),year=seq(2000,2200,by=20))

  return(preds.out)
}

trajectories<-as.list(rep(NA,24))
for(i in 1:length(scenario.names)){
  trajectories[[i]]<-getLandisTrajectory(paste0(scenario.names,"_")[i])
}
trajectories

trajectories<-bind_rows(trajectories)
trajectories$scenario<-rep(scenario.names,each=11)

trajectories<-subset(trajectories,year>2010)
trajectories<-subset(trajectories,year<2120)

traj_boot<-aggregate(bootpreds$Popsize,by=list(bootpreds$scenarioname,bootpreds$Year),FUN=quantile,probs=c(0.05,0.95))

traj_boot<-data.frame(scenario=traj_boot$Group.1,year=traj_boot$Group.2,lower=traj_boot$x[,1],upper=traj_boot$x[,2])
traj_boot$year<-as.numeric(as.character(traj_boot$year))


p3<-ggplot(data=newdf2)+
  geom_ribbon(aes(x=year,ymin=psizeLow,ymax=psizeUpp),fill="blue",alpha=0.2)+
  geom_line(aes(x=year,y=psize),color="blue")+
  geom_vline(xintercept = 2007,  size=0.5,linetype="dashed")+
  geom_errorbarh(aes(xmax=2040,xmin=2018,y=126000,height=20000))+
  guides(fill = guide_legend(reverse=TRUE),color=guide_legend(reverse=TRUE))+
  labs(x = "Years", y="Population size (males)", fill="",color="")+
  scale_y_continuous(labels = comma)+
  annotate("text", x= 2008, y = 400000,hjust=0, label = paste0("BRT model reference year (2007) \n", paste0("Estimate (90% CI) = ",round(popsize,0)," (",round(popsize.low,0)," - ",round(popsize.upp,0),")")),size=4)+
  annotate("text", x= 2029, y = 155000, label = ("Predictions with BBS long-term trend (2017): \n -1.78 (-3.560, 0.124)"),size=4)+
  geom_rug(data=yearsdata, aes(x=Year,y=Freq),length=unit(yearsdata$Freq/50000,"npc"),sides="b", alpha=0.5, position=jitter)+
  geom_errorbar(data=traj_boot,aes(x=year,ymin=lower,ymax=upper,color=scenario),position=position_dodge(width=3),width=20,alpha=0.5)+
  geom_errorbarh(aes(xmax=2100,xmin=2020,y=68000,height=20000))+
  annotate("text", x= 2050, y = 80000,angle=0, label = "Predictions with LANDIS-II simulations",size=4, hjust=0)+
  geom_line(data=subset(trajectories,year<2101),aes(x=year,y=mean,color=scenario))+
  guides(col=guide_legend(ncol=1))+
  theme(legend.text=element_text(size=8))+
  annotate("text", x= 2100, y = 300000,angle=0, label = "Bootstrap 90%CI (2100)",size=4, hjust=1)
p3

#save.image("D:/CHID regional NS BRT/BRT_outputs/outputs2.RData")

#load("D:/CHID regional NS BRT/BRT_outputs/outputs.RData")

#6. Calculate future trend 
futuretrend<-function(datafuture,datacurrent=newdf2,baseyear=2017){
  years<-2100-baseyear
  
  N0<-datacurrent[which(datacurrent$year==baseyear),"psize"]
  lower<-datacurrent[which(datacurrent$year==baseyear),"psizeLow"]
  upper<-datacurrent[which(datacurrent$year==baseyear),"psizeUpp"]
  
  N1<-datafuture
  
  trend<-100*((N1/N0)^(1/years)-1)
  trendLow<-100*((N1/upper)^(1/years)-1)
  trendUpp<-100*((N1/lower)^(1/years)-1)
  
  return(data.frame(trend=trend,lowerCL=trendLow,upperCL=trendUpp))
}

## Estimates
futuretrendsSims<-futuretrend(datafuture = scenarios2100$mean)
futuretrendsSims$scenario<-as.factor(scenario.names)
futuretrendsSims

## bootstrap samples
boot_preds_2100_df<-rbind(baseline_BudwormBaselineFire_2100,
                           baseline_BudwormBaselineFireEBFMHarvest_2100,
                           baseline_BudwormBaselineFireHighCPRSHarvest_2100,
                           baseline_BudwormBaselineFireHighPartialHarvest_2100,
                           baseline_BudwormBaselineFireLowCPRSHarvest_2100,
                           baseline_BudwormBaselineFireLowPartialHarvest_2100,
                           RCP26_GrowthBudwormProjectedFire_2100,
                           RCP26_GrowthBudwormProjectedFireEBFMHarvest_2100,
                           RCP26_GrowthBudwormProjectedFireHighCPRSHarvest_2100,
                           RCP26_GrowthBudwormProjectedFireHighPartialHarvest_2100,
                           RCP26_GrowthBudwormProjectedFireLowCPRSHarvest_2100,
                           RCP26_GrowthBudwormProjectedFireLowPartialHarvest_2100,
                           RCP45_GrowthBudwormProjectedFire_2100,
                           RCP45_GrowthBudwormProjectedFireEBFMHarvest_2100,
                           RCP45_GrowthBudwormProjectedFireHighCPRSHarvest_2100,
                           RCP45_GrowthBudwormProjectedFireHighPartialHarvest_2100,
                           RCP45_GrowthBudwormProjectedFireLowCPRSHarvest_2100,
                           RCP45_GrowthBudwormProjectedFireLowPartialHarvest_2100,
                           RCP85_GrowthBudwormProjectedFire_2100,
                           RCP85_GrowthBudwormProjectedFireEBFMHarvest_2100,
                           RCP85_GrowthBudwormProjectedFireHighCPRSHarvest_2100,
                           RCP85_GrowthBudwormProjectedFireHighPartialHarvest_2100,
                           RCP85_GrowthBudwormProjectedFireLowCPRSHarvest_2100,
                           RCP85_GrowthBudwormProjectedFireLowPartialHarvest_2100) 





fut.trend.df<-futuretrend(datafuture=boot_preds_2100_df$Popsize)
fut.trend.df

fut.trend.df$scenario<-as.factor(rep(scenario.names,each=1250))

fut.trends.plot<-ggplot(fut.trend.df,aes(x=scenario,y=trend))+geom_boxplot()+coord_flip()+geom_hline(yintercept=0, colour="blue")+geom_point(data=futuretrendsSims,aes(x=scenario,y=trend),shape=3,color="red")
fut.trends.plot


futuretrend2<-function(scenarioindex,datafuture=scenarios2100,datacurrent=newdf2,baseyear=2017){
  years<-2100-baseyear

  N0<-datacurrent[which(datacurrent$year==baseyear),"psize"]
  lower<-datacurrent[which(datacurrent$year==baseyear),"psizeLow"]
  upper<-datacurrent[which(datacurrent$year==baseyear),"psizeUpp"]

  N1<-datafuture[scenarioindex,"mean"]
  N1low<-datafuture[scenarioindex,"bootLower"]
  N1upp<-datafuture[scenarioindex,"bootUpper"]

  trend<-100*((N1/N0)^(1/years)-1)
  trendLow<-100*((N1low/upper)^(1/years)-1)
  trendUpp<-100*((N1upp/lower)^(1/years)-1)

  return(data.frame(trend=trend,lowerCL=trendLow,upperCL=trendUpp))
}


fut.trend.df2<-as.data.frame(matrix(unlist(lapply(1:24,FUN=futuretrend2)),ncol=3,nrow=24,byrow=T, dimnames=list(1:24,c("trend","lowerCL","upperCL"))))

fut.trend.df2$scenario<-as.factor(scenario.names)
fut.trend.df2

write.csv(fut.trend.df2,"D:/CHID regional NS BRT/November 2019 meetings/future_trends.csv")

fut.trend.df2$maxmin<-3
fut.trend.df2$maxmin[(which.min(fut.trend.df2$trend))]<-1
fut.trend.df2$maxmin[(which.max(fut.trend.df2$trend))]<-2


fut.trends.plot2<-ggplot(fut.trend.df2,aes(x=trend,y=scenario,color=factor(maxmin)))+
  geom_point()+xlim(-4,4)+
  geom_errorbarh(data=fut.trend.df2,aes(xmax=upperCL,xmin=lowerCL,color=factor(maxmin)))+
  labs(x="Future trend (2017-2100)",y="Scenario")+
  scale_color_manual(labels=c("Minimum","Maximum","Others"),values=c("Red","Blue", "Black")) +
  theme(legend.title = element_blank())
fut.trends.plot2


# Estimate distribution function for future trend ####
# As we can see from plot above, CL are also not symetrical with the mean
# So need to use a distribution that can allow asymetrical shapes AND negative values
# enter the Exponentially Modified Gaussian (EMG) Distribution

library(emg)

findEMG<-function(x,mu.init,sigma.init,lambda.init,N=10){
  low<-as.numeric(x[2])
  upp<-as.numeric(x[3])
  mat<-expand.grid(mu=runif(N,min=min(mu.init),max=max(mu.init)),sigma=runif(N,min=sigma.init[1],max=sigma.init[2]),lambda=runif(N,min=lambda.init[1],max=lambda.init[2]))
  mat[,4:5]<-t(apply(mat,1,FUN = function(x){qemg(p=c(0.05,0.95),mu=x[1],sigma=x[2],lambda=x[3])}))
  mat$distances<-pointDistance(p1=mat[,4:5],p2=c(low,upp),lonlat = F)
  out1<-mat[which.min(mat$distances),]
  
  mu.init2<-unlist(c(out1[1]-(out1[1]*0.25),out1[1]+(out1[1]*0.25)))
  sigma.init2<-unlist(c(out1[2]-(out1[2]*0.25),out1[2]+(out1[2]*0.25)))
  lambda.init2<-unlist(c(out1[3]-(out1[3]*0.25),out1[3]+(out1[3]*0.25)))
  
  mat2<-expand.grid(mu=runif(N,min=min(mu.init2),max=max(mu.init2)),sigma=runif(N,min=min(sigma.init2),max=max(sigma.init2)),lambda=runif(N,min=lambda.init2[1],max=lambda.init2[2]))
  mat2[,4:5]<-t(apply(mat,1,FUN = function(x){qemg(p=c(0.05,0.95),mu=x[1],sigma=x[2],lambda=x[3])}))
  mat2$distances<-pointDistance(p1=mat2[,4:5],p2=c(low,upp),lonlat = F)
  
  mat3<-rbind(mat,mat2)
  colnames(mat3)[4:5]<-c("lowerCL_0.5","upperCL_0.95")
  return(mat3[which.min(mat3$distances),])
}

library(pbapply)
EMGparams2100<-pbapply(fut.trend.df2, MARGIN=1, FUN=findEMG,mu.init=c(-2,2),sigma.init=c(0.5,5),lambda.init=c(0.1,1),N=25)
EMGparams2100<-as.data.frame(matrix(unlist(EMGparams2100),ncol=6,nrow=24,byrow=T, dimnames=list(1:24,names(EMGparams2100[[1]]))))

#8. Generate plot for future trend distribution function
pop.prob.2100<-as.list(scenarios2100$scenario)
names(pop.prob.2100)<-scenarios2100$scenario

id.df.fut.trend.2100<-data.frame("Trend"=NA,"Probability"=NA,"scenario"=NA)
for(i in 1:length(pop.prob.2100)){
  fut.trend.range<-seq(-4,5,length.out=350)
  id.dist.fut.trend.2100<-pemg(fut.trend.range,mu=EMGparams2100$mu[i],sigma=EMGparams2100$sigma[i],lambda=EMGparams2100$lambda[i],lower.tail = F)
  
  id.df.fut.trend.2100[(1:length(fut.trend.range))+length(fut.trend.range)*(i-1),]<-cbind(fut.trend.range,id.dist.fut.trend.2100,names(pop.prob.2100)[i])
}

id.df.fut.trend.2100$Trend<-as.numeric(id.df.fut.trend.2100$Trend)
id.df.fut.trend.2100$Probability<-as.numeric(id.df.fut.trend.2100$Probability)
str(id.df.fut.trend.2100)

trend.prob.2100.plot<-ggplot(id.df.fut.trend.2100, aes(x = Trend, y = Probability,color=scenario)) +
  geom_point(alpha=0.6)+
  geom_vline(mapping=aes(xintercept = trend,color=scenario),data=fut.trend.df2, size=0.8,alpha=0.4)+
  #annotate("text", x= 55000, y = 0.93,angle=0, label = paste0("Future trend predictions \n (mean = ",round(mean(fut.trend.df2$trend),0),")"),size=4, hjust=0)+
  #annotate("text", x= 150000, y = 0.93,angle=0, label = paste0("Scenario: ",scenarios2100$scenario[i]),size=4)+
  xlim(-4,5)+
  ylim(0,1)+
  scale_color_discrete(guide=guide_legend(ncol=1))+
  #scale_x_continuous(labels = scales::comma)+
  labs(x="T",y="Probability (Trend >= T)")
trend.prob.2100.plot

### Summarize current predictions by land cover type and CTI ####

##MODIS-based landcover (250-m)
nalc2005<-raster("M:/DataStuff/SpatialData/LCC05_NALCMS2010/Land_Cover_MXD/NA_LandCover_2005/data/NA_LandCover_2005/NA_LandCover_2005_LCC.img")
nalc2005<-projectRaster(nalc2005,rast2, method="ngb")
lc<-crop(nalc2005,rast2)
lc<-mask(nalc2005,rast2)

pred<-getValues(rast2)
predUpp<-getValues(upper2)
predLow<-getValues(lower2)

lcv<-getValues(lc)
lcv<-as.factor(as.character(lcv))
table(lcv)

pred1<-data.frame(pred=pred,upper=predUpp,lower=predLow,class=lcv)
pred1$count<-1

lcdens<-aggregate(pred1[,1:3], by=list(pred1$class),FUN="mean", na.rm=T)
names(lcdens)<-c("class", "mean_density","upper", "lower")

ltnalc <- read.csv("C:/Users/voeroesd/Documents/Repos/bamanalytics/lookup/nalcms.csv")
lcdens$LCclass<-ltnalc$Label[match(lcdens$class, ltnalc$Value)]
lcdens$LCclass<-factor(as.character(lcdens$LCclass))
lcdens$LCclassName<-ltnalc$Class_name[match(lcdens$class, ltnalc$Value)]
lcdens$LCclassName<-factor(as.character(lcdens$LCclassName))
lcdens


lc1<-aggregate(pred1$count, by=list(pred1$class),FUN="sum", na.rm=T)
names(lc1)<-c("class","area")
lc1$areaprop<-lc1$area/sum(lc1$area) 
lc1

lcdens$area<-lc1$area
lcdens$areaprop<-lc1$areaprop
lcdens

p<-ggplot(lcdens,aes(x=LCclassName,y=mean_density))+
  geom_bar(aes(fill=LCclassName),width=lcdens$areaprop,alpha=0.3, stat="identity")+
  geom_pointrange(aes(ymin=lower,ymax=upper,color=LCclassName))+
  xlab("Land cover class")+
  ylab("Mean density (males/ha)")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  labs(fill = "Land cover class",color="Land cover class")
p



# CTI
CTI<-raster("D:/CTI/CTI_NS.tif")
CTI<-projectRaster(CTI,rast2)
CTI<-mask(CTI,rast2)

CTIv<-getValues(CTI)

pred1$CTIq<-NA
pred1$CTIq[which(CTIv<6.9)]<-"Q1"
pred1$CTIq[which(6.9<CTIv&CTIv<7)]<-"Q2"
pred1$CTIq[which(7<CTIv&CTIv<8)]<-"Q3"
pred1$CTIq[which(CTIv>8)]<-"Q4"

CTI_habdens<-aggregate(pred1[,1:3], by=list(pred1$class,pred1$CTIq),FUN="mean", na.rm=T)
names(CTI_habdens)<-c("class","CTIq" ,"mean_density","upper", "lower")
CTI_habdens$LCclass<-ltnalc$Label[match(CTI_habdens$class, ltnalc$Value)]
CTI_habdens$LCclass<-factor(as.character(CTI_habdens$LCclass))
CTI_habdens$LCclassName<-ltnalc$Class_name[match(CTI_habdens$class, ltnalc$Value)]
CTI_habdens$LCclassName<-factor(as.character(CTI_habdens$LCclassName))
CTI_habdens

lc2<-aggregate(pred1$count, by=list(pred1$class,pred1$CTIq),FUN="sum", na.rm=T)
names(lc2)<-c("class","CTIq","area")
lc2$areaprop<-lc2$area/sum(lc2$area) 
lc2

CTI_habdens$area<-lc2$area
CTI_habdens$areaprop<-lc2$areaprop
CTI_habdens$id<-paste0(CTI_habdens$LCclassName,CTI_habdens$CTIq)

pCTI<-ggplot(CTI_habdens,aes(x=id,y=mean_density))+
  geom_bar(aes(fill=LCclassName),alpha=0.3, stat="identity")+
  geom_errorbar(aes(ymin=lower,ymax=upper,color=LCclassName,linetype=CTIq),size=0.8,position = position_dodge(width = 1))+
  #geom_point(aes(color=LCclassName))+
  xlab("Land cover class")+
  ylab("Mean density (males/ha)")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  scale_x_discrete(limits = rev(levels(CTI_habdens$id)))+
  labs(fill = "Land cover class",color="Land cover class",linetype="CTI quartiles")
pCTI


# Assessing sampling representativeness ####
# The method measures the similarity of any given point to a reference set ofpoints, with respect to the chosen predictor variables. It reports the closeness of the point to the distribution of reference points, gives negative values for dissimilar points and maps these values across the whole prediction region.



mess1<-mess(pred_abs_2011[[brt5$var.names]], extract(pred_abs_2011[[brt5$var.names]],cbind(datcombo$X,datcombo$Y)))

plot(mess1)

writeRaster(mess1, filename="D://CHID regional NS BRT/BRT_outputs/GFsigma250m/mess_brt_predsGF250m.tif", format="GTiff",overwrite=TRUE)


# Plot current model predictions and CIs ####
png("D://CHID regional NS BRT/November 2019 meetings/preds_allpoints.png",width=1080,height=700)
par(mar=c(0.5, 0.5, 0.5, 0.5))
plot(rast,  zlim = c(0, (0.15)),xaxt="n",yaxt="n",legend.width=1,legend.args = list(text = 'Density (males/ha)',side=4,line=2.5,cex=1,font=2))
points(datcombo_sp,pch=16,col=alpha(1,0.5), cex = 0.8)
dev.off()

pointsplot<-subset(datcombo_sp,ABUND>0)

png("D://CHID regional NS BRT/November 2019 meetings/preds_allpresences.png",width=1080,height=700)
par(mar=c(0.5, 0.5, 0.5, 0.5))
plot(rast, zlim = c(0, (0.15)),xaxt="n",yaxt="n",legend.width=1,legend.args = list(text = 'Density (males/ha)',side=4,line=2.5,cex=1,font=2))
points(pointsplot,pch=16,col=alpha(1,0.5), cex = 0.8)
dev.off()

png("D://CHID regional NS BRT/November 2019 meetings/preds_nopoints.png",width=1080,height=700)
par(mar=c(0.5, 0.5, 0.5, 0.5))
plot(rast,  zlim = c(0, (0.15)),xaxt="n",yaxt="n",legend.width=1,legend.args = list(text = 'Density (males/ha)',side=4,line=2.5,cex=1,font=2))
dev.off()

png("D://CHID regional NS BRT/November 2019 meetings/predsLower_nopoints.png",width=1080,height=700)
par(mar=c(0.5, 0.5, 0.5, 0.5))
plot(lower,  zlim = c(0, (0.15)),xaxt="n",yaxt="n",legend.width=1,legend.args = list(text = 'Density (males/ha)',side=4,line=2.5,cex=1,font=2))
dev.off()

png("D://CHID regional NS BRT/November 2019 meetings/predsUpper_nopoints.png",width=1080,height=700)
par(mar=c(0.5, 0.5, 0.5, 0.5))
plot(upper,  zlim = c(0, (0.15)),xaxt="n",yaxt="n",legend.width=1,legend.args = list(text = 'Density (males/ha)',side=4,line=2.5,cex=1,font=2))
dev.off()

png("D://CHID regional NS BRT/November 2019 meetings/predsLANDISarea_nopoints.png",width=1080,height=700)
par(mar=c(0.5, 0.5, 0.5, 0.5))
plot(rast2,  zlim = c(0, (0.15)),xaxt="n",yaxt="n",legend.width=1,legend.args = list(text = 'Density (males/ha)',side=4,line=2.5,cex=1,font=2))
dev.off()