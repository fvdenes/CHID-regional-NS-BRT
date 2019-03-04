library(raster)
library(rpart)
library(maptools)
library(data.table)
library(rgdal)
library(dplyr)
library(blockCV)
library(sp)
library(dismo)

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
    
    data_sp <-SpatialPointsDataFrame(coords = data[, 34:35], data = data, proj4string = LCC)
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

specoff <- filter(offl, SPECIES==as.character(speclist[j]))
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
datcombo_sp <-SpatialPointsDataFrame(coords = datcombo[, 34:35], data = datcombo, proj4string = LCC)


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
  "watNNS_Gauss750"
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
brt1<- brt_blocks(data=datcombo,pred.variables = pred.variables, lr=0.01,tc=3, output.folder = "D://CHID regional NS BRT/BRT_outputs/full_model/", blocks=sp_block_full, save.points.shp = TRUE)  
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
  "Species_Lari_Lar_v1_Gauss250",#
  #"Species_Lari_Lar_v1_Gauss750",#
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
  "Structure_Stand_Age_v1",
  #"Structure_Stand_Age_v1_Gauss250",
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
  "watNNS_Gauss250"
  #,"watNNS_Gauss750"
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
brt2<- brt_blocks(data=datcombo,pred.variables = pred.variables2, lr=0.01,tc=3, output.folder = "D://CHID regional NS BRT/BRT_outputs/selected_scales/", blocks=sp_block_select, save.points.shp = TRUE)  
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
  "watNNS"
  #,
  #"watNNS_Gauss250",
  #"watNNS_Gauss750"
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
brt3<- brt_blocks(data=datcombo,pred.variables = pred.variables3, lr=0.00001,tc=2, output.folder = "D://CHID regional NS BRT/BRT_outputs/cell_level/", blocks=sp_block_cell, save.points.shp = FALSE)  
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
  #"roadsNNS",
  "roadsNNS_Gauss250",
  #"roadsNNS_Gauss750",
  #"urbagNNS",
  "urbagNNS_Gauss250",
  #"urbagNNS_Gauss750",
  #"watNNS",
  "watNNS_Gauss250"
  #,
  #"watNNS_Gauss750"
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
brt4<- brt_blocks(data=datcombo,pred.variables = pred.variables4, lr=0.01,tc=3, output.folder = "D://CHID regional NS BRT/BRT_outputs/GFsigma250m/", blocks=sp_block_GF250, save.points.shp = TRUE)  
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
  #"roadsNNS",
  #"roadsNNS_Gauss250",
  "roadsNNS_Gauss750",
  #"urbagNNS",
  #"urbagNNS_Gauss250",
  "urbagNNS_Gauss750",
  #"watNNS",
  #"watNNS_Gauss250",
  "watNNS_Gauss750"
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
brt5<- brt_blocks(data=datcombo,pred.variables = pred.variables5, lr=0.01,tc=3, output.folder = "D://CHID regional NS BRT/BRT_outputs/GFsigma750m/", blocks=sp_block_GF750, save.points.shp = TRUE)  
end_time<-Sys.time()
end_time-start_time

save.image("D:/CHID regional NS BRT/BRT_outputs/outputs.RData")

# Comparison of the different models ####
library(ggplot2)


# Model deviance

df<- data.frame(grp=c("Full model",
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
df
k<-ggplot(df, aes(grp,deviance,ymin=deviance-se,ymax=deviance+se))
k+geom_pointrange()

# correlation
df2<- data.frame(grp=c("Full model",
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
df2
k2<-ggplot(df2, aes(grp,correlation,ymin=correlation-se,ymax=correlation+se))
k2+geom_pointrange()

# Calibration:  The first two statistics were the estimated intercepts and slopes of linear regression models of predictions against observations. The intercept measures the magnitude and direction of bias, with values close to 0 indicating low or no bias. The slope yields information about the consistency in the bias as a function of the mean, with a value of 1 indicating a consistent bias if the intercept is a nonzero value.
## intercept
df3<- data.frame(grp=c("Full model",
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
df3
k3<-ggplot(df3, aes(grp,calibration.intercept,ymin=calibration.intercept-se,ymax=calibration.intercept+se))
k3+geom_pointrange()

## slope
df4<- data.frame(grp=c("Full model",
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
df4
k4<-ggplot(df4, aes(grp,calibration.slope,ymin=calibration.slope-se,ymax=calibration.slope+se))
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

dev_plot(brt4)


# Population size from density predictions ####
#pred_abs_2011<- brick("D:/CHID regional NS BRT/prediction dataset/abs2011_250m.grd")

library(gbm)

rast <-  predict(pred_abs_2011,
          brt5,
          type = "response",
          n.trees = brt5$n.trees)
plot(rast)

preds_brt5<-getValues(rast)[!is.na(getValues(rast))]

density_brt5_250m<-preds_brt5*6.25

abundance_brt5<-matrix(NA,length(density_brt5_250m),100)

for(i in 1:length(density_brt5_250m)){
    abundance_brt5[i,]<- rpois(100,density_brt5_250m[i])
}


mean(colSums(abundance_brt5))

mean(colSums(abundance_brt5))*2

sd(colSums(abundance_brt5))

sum(density_brt5_250m) 

sum(getValues(rast), na.rm= TRUE)


upper<-raster("D://CHID regional NS BRT/BRT_outputs/GFsigma750m/confint/ confint_upper.tif")
upper_brt5<-getValues(upper)[!is.na(getValues(upper))]
upper_brt5_250m<-upper_brt5*6.25
abundance_upper_brt5<-matrix(NA,length(upper_brt5_250m),100)
for(i in 1:length(upper_brt5_250m)){
  abundance_upper_brt5[i,]<- rpois(100,upper_brt5_250m[i])
}
mean(colSums(abundance_upper_brt5))
mean(colSums(abundance_upper_brt5))*2

lower<-raster("D://CHID regional NS BRT/BRT_outputs/GFsigma750m/confint/ confint_lower.tif")
lower_brt5<-getValues(lower)[!is.na(getValues(lower))]
lower_brt5_250m<-lower_brt5*6.25
abundance_lower_brt5<-matrix(NA,length(lower_brt5_250m),100)
for(i in 1:length(lower_brt5_250m)){
  abundance_lower_brt5[i,]<- rpois(100,lower_brt5_250m[i])
}
mean(colSums(abundance_lower_brt5))
mean(colSums(abundance_lower_brt5))*2

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
    
    
    sample<-sample(1:nrow(data),size=nrow(data),replace=T)
    data2<-data[sample,]
    x1<- try(brt <-
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
    
    if(class(x1)=="NULL"){
      sample<-sample(1:nrow(data),size=nrow(data),replace=T)
      x1<- try(brt <-
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
    
    if(class(x1)=="NULL"){
      sample<-sample(1:nrow(data),size=nrow(data),replace=T)
      x1<- try(brt <-
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
    
    if(brt$n.trees>9999){
      sample<-sample(1:nrow(data),size=nrow(data),replace=T)
      x1<- try(brt <-
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
    
    writeRaster(stack, filename=paste0(z, "samples.grd"), format="raster",overwrite=TRUE)
    
    
    
  }
  
  fun0.05 <- function(x) {quantile(x, probs = 0.05, na.rm = TRUE)}
  lower<- calc(stack[[-1]],fun0.05)
  fun0.95 <- function(x) {quantile(x, probs = 0.95, na.rm = TRUE)}
  upper<- calc(stack[[-1]],fun0.95)
  
  writeRaster(lower, filename=paste0(z, " confint_lower.tif"), format="GTiff",overwrite=TRUE)
  writeRaster(upper, filename=paste0(z, " confint_upper.tif"), format="GTiff",overwrite=TRUE)
  
  return(stack)
}

start_time <- Sys.time()
confintbrt5<-boot_brt(datcombo,brt5,sp_block_GF750,pred_abs_2011,nsamples=250,output.folder = "D://CHID regional NS BRT/BRT_outputs/GFsigma750m/confint/")
end_time <- Sys.time()
end_time - start_time

save.image("D:/CHID regional NS BRT/BRT_outputs/outputs.RData")

names(confintbrt5)

#### Residual trend ####
# the first 7290 rows in datcombo are from 2005 and earlier period. the rest is from 2006 and later.
datcombo$period<-"2005 and earlier"
datcombo$period[7291:15059]<-"2006 and later"

fit<-log(brt5$fitted)

m <- glm(datcombo$ABUND ~ datcombo$period, family=poisson, offset=fit+datcombo$logoffset)
trend <-  100*(exp(coef(m)[2])-1)


# Assessing sampling representativeness ####
# The method measures the similarity of any given point to a reference set ofpoints, with respect to the chosen predictor variables. It reports the closeness of the point to the distribution of reference points, gives negative values for dissimilar points and maps these values across the whole prediction region.



mess1<-mess(pred_abs_2011[[brt5$var.names]], extract(pred_abs_2011[[brt5$var.names]],cbind(datcombo$X,datcombo$Y)))

plot(mess1)

writeRaster(mess1, filename="D://CHID regional NS BRT/BRT_outputs/GFsigma250m/mess_brt_predsGF250m.tif", format="GTiff",overwrite=TRUE)




