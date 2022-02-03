#Read in libraries
library(dplyr)
library(tidyr)
library(dismo)
library(raster)
library(data.table)
library(MASS) 
library(magrittr) 
library(maptools) 
require(rJava)
library(sf)
library(stringr)
library(ggplot2)
library(spatialEco)
library(rgdal)
library(sjmisc)
library(rgeos)
library(GISTools)
library(sdmvspecies)
library(Rfast)
library(png)
library(patchwork)
library(Rfast)
library(maps)
library(GISTools)
library(PerformanceAnalytics)
library(usdm)
library(psych)

#Shapefile location 
shapeDir <- "C:\\Users\\jeffr\\Documents\\Academic\\Stanford\\CR_birds_project\\PNAS_revision\\shapefiles"

#Raw tifs
rawTiffs <- "C:\\Users\\jeffr\\Documents\\Academic\\Stanford\\CR_birds_project\\PNAS_revision\\Tiffs_raw"


#Raster location
rasterDir <- "C:\\Users\\jeffr\\Documents\\Academic\\Stanford\\CR_birds_project\\PNAS_revision\\Tiffs"

#SDM direcotry 
sdmOut <- "C:\\Users\\jeffr\\Documents\\Academic\\Stanford\\CR_birds_project\\PNAS_revision\\SDMs"

#Tourism Maxnet out path 
maxentOut <- "C:\\Users\\jeffr\\Documents\\Academic\\Stanford\\CR_birds_project\\PNAS_revision\\tourismMaxent"



#Plotting multiple GGplot graphs together
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#Pre-process environmental layers
abioticPredictors <- function(){
  
  #Projection systems
  wgsProj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  utmProj <- newproj <- "+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs"
  
  
  #CR shapefile
  cr <- rgdal::readOGR(file.path(shapeDir, "CR.shp"))
  
  #Clip and reproject existing rasters 
  cAndR <- function(file, name, baseMap = NULL, band = NULL){
    
    #Read the band you want or the raster if there is only one band
    if(is.null(band)){
      rasterF <- raster(file)
    }else{
      rasterF <- raster(file, band = band)}
    
    #Deal with no data
    rasterF[rasterF==-9999] <- NA
    
    #Specify that is projection is missing it is WGS (for bioclim)
    if(is.na(as.character(rasterF@crs))){
      crs(rasterF) <- wgsProj
    }
    
    #Reproject and resample
    if(is.null(baseMap)){
      newRaster <- projectRaster(rasterF, crs = newproj, method = 'bilinear', res = 915)
    }else{
      baseMapR <- raster(baseMap)
      newRaster <- projectRaster(rasterF, baseMapR, crs = newproj, method = 'bilinear', res = 915)
    }
    
    #Clip to CR
    cr <- spTransform(cr, newproj)
    rasterF <- crop(rasterF, cr)
    rasterF <- mask(rasterF, cr)
    
    #Save raster
    writeRaster(file.path(rasterDir, newRaster), name, overwrite = T)
  }
  
  #Process the raw tiffs 
  cAndR(file.path(rawTiffs, 'bioclim', 'bio_1.bil'),  
        file.path(rasterDir ,'bio1.tif'))
  
  cAndR(file.path(rawTiffs, 'bioclim', 'bio_4.bil'),  
        file.path(rasterDir ,'bio4.tif'), file.path(rasterDir ,'bio1.tif'))
  
  cAndR(file.path(rawTiffs, 'bioclim', 'bio_12.bil'),   
        file.path(rasterDir ,'bio12.tif'), file.path(rasterDir ,'bio1.tif'))
  
  cAndR(file.path(rawTiffs, 'bioclim', 'bio_15.bil'),  
        file.path(rasterDir ,'bio15.tif'), file.path(rasterDir ,'bio1.tif'))
  
  cAndR(file.path(rawTiffs, 'rasters_chris', 'tree-cover.tif'),  
        file.path(rasterDir ,'tree-cover.tif'), file.path(rasterDir ,'bio1.tif'))
  
  cAndR(file.path(rawTiffs, 'rasters_chris', 'fcover.tif'),  
        file.path(rasterDir ,'fcover1.tif'), file.path(rasterDir ,'bio1.tif'), 1)
  
  cAndR(file.path(rawTiffs, 'rasters_chris', 'fcover.tif'),  
        file.path(rasterDir ,'fcover2.tif'), file.path(rasterDir ,'bio1.tif'), 2)
  
  cAndR(file.path(rawTiffs, 'rasters_chris', 'fcover.tif'),  
        file.path(rasterDir ,'fcover3.tif'), file.path(rasterDir ,'bio1.tif'), 3)
  
  
  
  
  #Read in roads
  roads <- rgdal::readOGR(file.path(shapeDir, "groads-v1-americas-shp\\groads-v1-americas-shp\\gROADS-v1-americas.shp"))
  #Interesect
  roads <- intersect(roads, cr)
  #Reproject
  roads <- spTransform(roads, newproj)
  #Rasterive and save
  roadRaster <- rasterize(roads, raster(file.path(rasterDir ,'bio1.tif')), field = 1)
  roadRaster <- distance(roadRaster)
  roadRaster <- roadRaster + (raster(file.path(rasterDir ,'bio1.tif')) * 0)
  writeRaster(roadRaster, file.path(rasterDir ,'distRoads.tif'), overwrite = T)
  
  
  #Read in lakes, clip to CR, transform, and create distance layer 
  sFile <- readShapeSpatial(file.path(shapeDir, "CRI_wat\\CRI_water_areas_dcw.shp"))
  sFile <- intersect(sFile, cr)
  sFile <- spTransform(sFile, wgsProj)
  toSave <- rasterize(sFile, file.path(rasterDir ,'bio1.tif'), field = 1, update = F)
  toSave <- distance(toSave)
  
  #Read in rivers, clip to CR, transform, and create distance layer 
  sFile <-  readShapeSpatial(file.path(shapeDir, "CRI_wat\\CRI_water_lines_dcw.shp"))
  sFile <- intersect(sFile, cr)
  sFile <- spTransform(sFile, newproj)
  toSave2 <- rasterize(sFile, file.path(rasterDir ,'bio1.tif'), field = 1, update = F)
  toSave2 <- distance(toSave2)
  
  #Read in coast,  transform, and create distance layer 
  sFile <- readShapeSpatial(file.path(shapeDir, "crCoast.shp"))
  sFile <- spTransform(sFile, newproj)
  toSave3 <- rasterize(sFile, file.path(rasterDir ,'bio1.tif'), field = 1, update = F)
  toSave3 <- distance(toSave3)
  
  #Take minimum distance and save
  toSave <- min(toSave, toSave2, toSave3)
  toSave <- toSave + (file.path(rasterDir ,'bio1.tif') * 0)
  writeRaster(toSave, file.path(rasterDir ,'distWater.tif'), overwrite = T)
  
  
  #Read in PAs, transform, calculate distance and save 
  sFile <- rgdal::readOGR(file.path(shapeDir, "CR_protected_areas\\WDPA_Dec2018_CRI-shapefile-polygons.shp"))
  sFile <- spTransform(sFile, newproj)
  toSave <- rasterize(sFile, file.path(rasterDir ,'bio1.tif'), field = 1, update = F)
  toSave <- distance(toSave)
  toSave <- toSave + (file.path(rasterDir ,'bio1.tif') * 0)
  writeRaster(toSave, file.path(rasterDir ,'distPAs.tif'), overwrite = T)
  
  #Read in hotels
  hotels <- read.csv(file.path(shapeDir, "Hotels_revised_duplicates_deleted.csv"))
  hotels <- hotels[c('Name', "Longitude", "Latitude")]
  coordinates(hotels)<-~Longitude+Latitude
  proj4string(hotels) <- CRS(wgsProj)
  
  #reproject
  hotels <- spTransform(hotels, newproj)
  
  #Set up raster
  bm <- raster(file.path(rasterDir ,'distRoads.tif'))
  rest <- c(bm@nrows, bm@ncols)
  lims <- c(bm@extent@xmin, bm@extent@xmax, bm@extent@ymin, bm@extent@ymax)
  
  #Do kernal density 
  hotelsR <- kde2d(hotels$Longitude, hotels$Latitude, h=25000, n = rest, lims = lims)
  hotelsR <- raster(hotelsR)
  proj4string(hotelsR) <- newproj
  
  #Resammple and save 
  hotelsR <- projectRaster(hotelsR, bm, crs = newproj, method = 'bilinear', res = 915)
  hotelsR <- hotelsR + bm * 0
  hotelsR@extent <- bm@extent
  writeRaster(hotelsR, file.path(rasterDir ,'hotelDensity.tif'), overwrite = T)
}

#Line to run code to pre-process rasters
#abioticPredictors()

#Run SDMs
doSDMs <- function(){
  #Projection systems
  wgsProj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  utmProj <- "+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs"
  
  
  #Read in CSV of GBIF occurances 
  colToUse <- c('year', 'kingdom', 'phylum', 'class', 'order', 'genus', 'species', 'decimalLatitude', 'decimalLongitude', 'coordinateUncertaintyInMeters', 'coordinatePrecision', 'basisOfRecord', 'issue')
  raw <- fread(file.path(shapeDir, "raw2.csv"), select = colToUse)
  
  #Exclude anything without lat/long
  raw <- raw[!is.na(raw$decimalLongitude),]
  raw <- raw[!is.na(raw$decimalLatitude),]
  
  #Exclude fossils
  raw <- raw[(raw$basisOfRecord != 'FOSSIL_SPECIMEN'),]
  
  #Exclude records with errors
  raw <- raw[!str_detect(raw$issue, 'PRESUMED_NEGATED_LONGITUDE'),]
  raw <- raw[!str_detect(raw$issue, 'PRESUMED_NEGATED_LATITUDE'),]
  raw <- raw[!str_detect(raw$issue, 'COORDINATE_REPROJECTION_SUSPICIOUS'),]
  raw <- raw[!str_detect(raw$issue, 'COORDINATE_PRECISION_INVALID'),]
  raw <- raw[!str_detect(raw$issue, 'COORDINATE_ROUNDED'),]
  raw <- raw[!str_detect(raw$issue, 'GEODETIC_DATUM_INVALID'),]
  
  
  #Excludea nything before 2020
  raw <- raw[(raw$year >= 2000),]
  
  
  #DF formatting
  raw <- as.data.frame(raw)
  coordinates(raw)<-~decimalLongitude+decimalLatitude
  
  #Reproject
  proj4string(raw) <- CRS(wgsProj)
  raw <- spTransform(raw, utmProj)
  
  
  
  #CR shapefile
  cr <- readShapePoly(file.path(shapeDir, "CR.shp"))
  proj4string(cr) <- CRS(wgsProj)
  cr <- spTransform(cr, newproj)
  
  
  
  #Make bias file
  climdat <- brick(file.path(rasterDir, 'bio1.tif'))
  occur.ras <- rasterize(raw@coords, climdat, 1)
  occur.states <- mask(occur.ras, cr) %>% crop(cr)
  
  
  
  
  presences <- which(values(occur.states) == 1)
  pres.locs <- coordinates(occur.states)[presences, ]
  pres.locs <- sample(nrow(pres.locs), 100000)
  
  
  rest <- c(climdat@nrows, climdat@ncols)
  lims <- c(climdat@extent@xmin, climdat@extent@xmax, climdat@extent@ymin, climdat@extent@ymax)
  
  dens <- kde2d(pres.locs[,1], pres.locs[,2], n = rest, lims = lims, h = 5000)
  dens.ras <- raster(dens)
  
  
  
  
  
  dens.ras <- resample(dens.ras, climdat)
  dens.ras <- dens.ras + (raster(file.path(rasterDir, 'bio1.tif')) *0)
  
  
  #Write bias bile
  writeRaster(dens.ras, file.path(rasterDir, 'bias.tif'), overwrite = TRUE)
  
  #Use bias file to generate pseudo absences 
  bckrdPts <- as.data.frame(rasterToPoints(dens.ras))
  pseudos <- bckrdPts[sample(nrow(bckrdPts), 10000, replace = TRUE, prob = bckrdPts$layer), ]
  pseudos <- pseudos[c('x', 'y')]
  
  
  
  #Define directory to save output files
  dir.create(sdmOut, showWarnings = FALSE)
  
  #Predictors
  #first import all files in a single folder as a list 
  rastlist0 <- list.files(path = rasterDir, pattern='.tif$', all.files=TRUE, full.names=TRUE)
  rastlist <- rastlist0[str_detect(rastlist0,pattern="bio")]
  rastlist <- append(rastlist, rastlist0[str_detect(rastlist0,pattern="fcover")])
  rastlist <- append(rastlist, rastlist0[str_detect(rastlist0,pattern="tree")])
  rastlist <- rastlist[!str_detect(rastlist,pattern="fcover3")]
  allrasters <- stack(rastlist)
  
  #Sort by taxa
  #Loop through taxa
  taxaToLoop <- c('Amphibia') #, 'Aves', 'Mammalia', 'Reptilia')
  for (taxa in taxaToLoop){
    
    taxaDF <- raw[(raw@data$class == taxa),]
    
    
    
    #Create an output folder for just the taxon you're workng on
    
    taxaOut <- file.path(sdmOut, taxa)
    dir.create(taxaOut, showWarnings = FALSE)
    
    #Loop through the species in the taxon of interest 
    sppToLoop <- unique(taxaDF@data$species)
    for(k in seq(1, length(sppToLoop), 1)){
      
      try({
        #Select only occs from that species
        speciesName <- sppToLoop[k]
        sppDF <- taxaDF[(taxaDF@data$species == speciesName),]
        
        sppDF <- sppDF[cr,]
        sppDF <- sppDF@coords
        
        #Create folder for each species outputs 
        print(speciesName)
        fname <-file.path(taxaOut, speciesName)
        dir.create(fname, showWarnings = FALSE)
        
        #Execute the maxent run
        myModel <- maxent(allrasters, sppDF, a = pseudos)
        r <- predict(myModel, allrasters) 
        
        #Save prediction raster
        writeRaster(r, file.path(fname, 'meanPrev.tif'), overwrite=TRUE)
        
        #Summarize statistics from model
        summStats<- evaluate(sppDF, pseudos, myModel, allrasters)
        auc <- summStats@auc
        np <- summStats@np
        tv <- as.data.frame(myModel@results)['Equal.training.sensitivity.and.specificity.Cloglog.threshold',]
        
        
        #Save summary stats
        tS <- as.data.frame(transpose(list(c(auc, np, tv))))
        colnames(tS) <-    c('auc', 'points', 'threshold')
        write.csv(tS, file.path(fname, 'summary.csv'), quote = FALSE, row.names = FALSE)
        print(tS)
        print('scuess')
      })
    }    
  }
}

#Line to run biodiveristy SDMs
#doSDMs()

#Make summary biodiversity maps 
summaryBiodiversity <- function(){
  #Species Jim Zook said to exclude as they don't occur in CR
  toExclude <- c('Arremon torquatus',
                 'Basileuterus melanotis', 
                 'Campylorhynchus capistratus',
                 'Cantorchilus elutus',
                 'Colibri thalassinus',
                 'Dendrocolaptes certhia',
                 'Gallinula chloropus',
                 'Trogon aurantiiventris')
  
  
  

  makeMap <- function(bioDivDirectory, biodivAsc){
    #Get asciis for each speices
    medAscs <- list.files(path =bioDivDirectory, pattern = "*meanPrev.tif", full.names = TRUE, recursive = TRUE)
    summStats<- list.files(path =bioDivDirectory, pattern = ".csv", full.names = TRUE, recursive = TRUE)
    
    a = 0
    
    #Add up total specie richness
    for(i in c(1:length(medAscs))){
      
      sumCSV <- read.csv(summStats[i])
      sppName <- dirname(summStats[i])
      sppName <- tail(str_split(sppName, '/')[[1]],1)
      if(sppName %in% toExclude == F){
        if(sumCSV$auc >= 0.75 & sumCSV$points >= 25){
          
          data <- raster(medAscs[i])  
          data[data > sumCSV$threshold] <- 1
          data[data != 1] <- 0
          
          
          if (a == 0){
            data2 = data
            a = 1 
          }
          else{
            data2 = data  + data2
          }
        }
      }
      
    } 
    
    
    
    #Save that 
    writeRaster(data2, biodivAsc, overwrite=TRUE)
  }
  
  makeMap(file.path(sdmOut, 'Amphibia'), file.path(rasterDir, 'amphibians.tif'))
  makeMap(file.path(sdmOut, 'Mammalia'), file.path(rasterDir, 'mammals.tif'))
  makeMap(file.path(sdmOut, 'Aves'), file.path(rasterDir, 'abirdDiversity.tif'))
  makeMap(file.path(sdmOut, 'Reptilia'), file.path(rasterDir, 'reptiles.tif'))

  armDiversity <- raster(file.path(rasterDir, 'amphibians.tif')) + 
                  raster(file.path(rasterDir, 'mammals.tif')) + 
                  raster(file.path(rasterDir, 'reptiles.tif'))
  writeRaster(armDiversity, file.path(rasterDir, 'arm_diversity.tif'), overwrite = T)

  makeMapTE <- function(bioDivDirectory, biodivAsc){
    #Get asciis for each speices
    medAscs <- list.files(path =bioDivDirectory, pattern = "*meanPrev.tif", full.names = TRUE, recursive = TRUE)
    summStats<- list.files(path =bioDivDirectory, pattern = ".csv", full.names = TRUE, recursive = TRUE)
    
    threatenedEndangered <-read.csv('~/Academic/Stanford/CR_Birds_Project/Ale/288_endemic_or_threatened_birds.csv')
    threatenedEndangered <- threatenedEndangered$Endemic_or_threatened
    threatenedEndangered <- str_replace(threatenedEndangered, 'Aves_', '')
    threatenedEndangered <- str_replace(threatenedEndangered, '_binary', '')
    threatenedEndangered <- str_replace(threatenedEndangered, '_', ' ')
    
    
    a = 0
    
    #Add up total specie richness
    for(i in c(1:length(medAscs))){
      
      sumCSV <- read.csv(summStats[i])
      sppName <- dirname(summStats[i])
      sppName <- tail(str_split(sppName, '/')[[1]],1)
      if(sppName %in% threatenedEndangered == T){
        if(sumCSV$auc >= 0.75 & sumCSV$points >= 25){
          
          data <- raster(medAscs[i])  
          data[data > sumCSV$threshold] <- 1
          data[data != 1] <- 0
          
          
          if (a == 0){
            data2 = data
            a = 1 
          }
          else{
            data2 = data  + data2
          }
        }
      }
      
    } 
    
    
    
    #Save that 
    writeRaster(data2, biodivAsc, overwrite=TRUE)
  }  
  makeMapTE(file.path(sdmOut, 'Aves'), file.path(rasterDir, 'abirdDiversity_te.tif'))
}

#Line to run code to make maps that summarize biodiveirsty
#summaryBiodiversity()

#Summarize data
summarizeData <- function(){
  #Projection systems
  wgsProj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  utmProj <- "+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs"
  
  
  #Read in CR shapefile and transform to UTM
  cr <- rgdal::readOGR(file.path(shapeDir, 'cr.shp'))
  cr <- spTransform(cr, utmProj)
  
  
  #Read in flickr data
  flickrP <- as.data.frame(fread(file.path(shapeDir, 'flickr.csv')))
  
  flickrP <- flickrP[c('decimalLongitude', 'decimalLatitude')]
  coordinates(flickrP)<-~decimalLongitude+decimalLatitude
  
  #Transform to utm
  proj4string(flickrP) <- CRS(wgsProj)
  flickrP <- spTransform(flickrP, utmProj)
  
  #Clip to CR
  flickrP <- flickrP[cr,]
  print(paste('Number flickr points = ', as.character(length(flickrP))), sep ='')
  
  #Read in eBird data
  ebirdP <- as.data.frame(fread(file.path(shapeDir, "checklist_uploads.csv")))  
  
  #Clip years 
  ebirdP <- ebirdP[ebirdP$YEAR >= 2005, ]
  ebirdP <- ebirdP[ebirdP$YEAR <= 2017, ]
  
  #Reproject 
  ebirdP <- ebirdP[c('LONGITUDE', 'LATITUDE')]
  coordinates(ebirdP)<-~LONGITUDE+LATITUDE
  proj4string(ebirdP) <- CRS(wgsProj)
  ebirdP <- spTransform(ebirdP, utmProj)
  
  #Clip to CR
  ebirdP <- ebirdP[cr,]
  print(paste('Number ebird points = ', as.character(length(ebirdP))))
  
  
  #Read in and reproject PAs
  sFile <- rgdal::readOGR(file.path(shapeDir, "CR_protected_areas\\WDPA_Dec2018_CRI-shapefile-polygons.shp"))
  sFile <- spTransform(sFile, utmProj)
  sFile <- raster::intersect(sFile, cr)
  
  #Calculate area 
  sFile$area_sqkm <- raster::area(sFile) / 1000000
  
  #Define all the predictors you're going to use
  rastlist <- list.files(path = rasterDir, pattern='tif$', all.files=TRUE, full.names=TRUE)
  rastlist <- rastlist[!str_detect(rastlist,pattern="cover")]
  rastlist <- rastlist[!str_detect(rastlist,pattern="bio")]
  rastlist <- rastlist[!str_detect(rastlist,pattern="amphibian")]
  rastlist <- rastlist[!str_detect(rastlist,pattern="mammal")]
  rastlist <- rastlist[!str_detect(rastlist,pattern="reptile")]
  

  
  #Stack the rasters 
  allrasters <- stack(rastlist)
  
  
  a <- c()
  a <- rbind(a, cellStats(allrasters, 'mean'))
  a <- rbind(a, cellStats(allrasters, 'sd'))
  a <- rbind(a, cellStats(allrasters, 'median'))
  a <- rbind(a, cellStats(allrasters, 'min'))
  a <- rbind(a, cellStats(allrasters, 'max'))
  a <- cbind(c('mean', 'sd', 'median', 'min', 'max'), a)
  print('Cell stats for nation-wide layers')
  print(t(a))
  write.csv(t(a), file.path(maxentOut, 'nationwide_layers_summary.csv'))
  
  #Create masks for protected areas 
  pa_mask <- rasterize(sFile, allrasters$distRoads, field = 1, update = F)
  pa_size <- rasterize(sFile, allrasters$distRoads, field = 'area_sqkm', update = F)
  
  #Add the PA size layer 
  rename <- names(allrasters)
  rename <- append(rename, 'pa_size')
  allrasters <- stack(allrasters,pa_size)
  names(allrasters) <- rename
  

  
  
  #Function to score VIF
  scoreVIF <- function(names){
    vifRasters <- allrasters2[[names]]
    print(vif(vifRasters))
  }
  
  #Score the VIF runs for the various model runs 
  #National 
  allrasters2 <- allrasters
  print('National, birds')
  scoreVIF(c('abirdDiversity', 'distRoads', 'distWater', 'hotelDensity', 'distPAs'))
  print('National, birds_te')
  scoreVIF(c('abirdDiversity_te', 'distRoads', 'distWater', 'hotelDensity', 'distPAs'))
  print('National, ARM')
  scoreVIF(c('arm_diversity', 'distRoads', 'distWater', 'hotelDensity', 'distPAs'))
  
  #Protected areas
  allrasters2 <- allrasters * pa_mask
  names(allrasters2) <- names(allrasters)
  print('PAs, birds')
  scoreVIF(c('abirdDiversity', 'distRoads', 'distWater', 'hotelDensity', 'pa_size'))
  print('PAs, birds_te')
  scoreVIF(c('abirdDiversity_te', 'distRoads', 'distWater', 'hotelDensity', 'pa_size'))
  print('PAs, ARM')
  scoreVIF(c('arm_diversity', 'distRoads', 'distWater', 'hotelDensity', 'pa_size'))
  
  
  a <- c()
  allrasters2 <- allrasters * pa_mask
  names(allrasters2) <- names(allrasters)
  a <- rbind(a, cellStats(allrasters2, 'mean'))
  a <- rbind(a, cellStats(allrasters2, 'sd'))
  a <- rbind(a, cellStats(allrasters2, 'median'))
  a <- rbind(a, cellStats(allrasters2, 'min'))
  a <- rbind(a, cellStats(allrasters2, 'max'))
  a <- cbind(c('mean', 'sd', 'median', 'min', 'max'), a)
  print('Cell stats for PA layers')
  print(t(a))
  write.csv(t(a), file.path(maxentOut, 'PAs_layers_summary.csv'))
  
  
  #Make maps of varaibles 
  outFile <- file.path(maxentOut, "variable_maps.png")
  png(outFile, width = 10, height = 6, units = 'in', res = 600) 
  par(mar=c(2,2,2,5))
  names(allrasters) <- c('Bird_richness', 'Threatened_bird_richness', 'Amphibian_reptile_and_mammal_richness',
                         'Distance_to_PAs', 'Distance_to_roads', 'Distance_to_water',
                         'Hotel_density', 'Protected_area_size')
  
  plot(allrasters, axes=FALSE, box=FALSE)
  dev.off()
  
  #Print summary stats
  print(allrasters)
  cellStats(allrasters, mean)
  cellStats(allrasters, sd)
  cellStats(allrasters, median)
  
  #Create histogram for national variables 
  outFile <- file.path(maxentOut, "national_hist.png")
  png(outFile, width = 10, height = 6, units = 'in', res = 600)  
  allrasters2 <- allrasters[[names(allrasters)[names(allrasters)!='Protected_area_size']]]
  hist(allrasters2)
  dev.off()
  
  #Create correlation plot for national variables
  outFile <-file.path(maxentOut, "national_corr.png")
  png(outFile, width = 10, height = 10, units = 'in', res = 600)  
  flatDF <- as.matrix(allrasters2)
  flatDF <- flatDF[complete.cases(flatDF),]
  flatDF <- as.data.frame(flatDF)
  flatDF <- flatDF[sample(nrow(flatDF), 1000), ]
  chart.Correlation(flatDF, histogram=F, pch=12, method = 'spearman', cex = 4)
  dev.off()
  
  #Create histogram for just PAs
  outFile <- file.path(maxentOut, "PA_hist.png")
  png(outFile, width = 10, height = 6, units = 'in', res = 600)  
  allrasters2 <- allrasters[[names(allrasters)[names(allrasters)!='Distance_to_PAs']]]
  renameX <- names(allrasters2)
  toPlot <- (allrasters2*pa_mask)
  names(toPlot) <-renameX
  hist(toPlot)
  dev.off()
  
  #Create correlation plot for just PAs
  outFile <- file.path(maxentOut, "PA_corr.png")
  png(outFile, width = 10, height = 10, units = 'in', res = 600)  
  flatDF <- as.matrix(toPlot)
  flatDF <- flatDF[complete.cases(flatDF),]
  flatDF <- as.data.frame(flatDF)
  flatDF <- flatDF[sample(nrow(flatDF), 1000), ]
  chart.Correlation(flatDF, histogram=F, pch=12, method = 'spearman', cex = 4)
  dev.off()
  
}

#Execute function to make summary figures 
summarizeData()


#Code to create maxent tourism models 
runModels <- function(taxa, extent, baseOut){
  #Create directory to save results
  dir.create(baseOut)
  
  
  #Define all the predictors you're going to use
  rastlist <- list.files(path = rasterDir, pattern='tif$', all.files=TRUE, full.names=TRUE)
  rastlist <- rastlist[!str_detect(rastlist,pattern="cover")]
  rastlist <- rastlist[!str_detect(rastlist,pattern="bio")]
  rastlist <- rastlist[!str_detect(rastlist,pattern="amphibian")]
  rastlist <- rastlist[!str_detect(rastlist,pattern="mammal")]
  rastlist <- rastlist[!str_detect(rastlist,pattern="reptile")]
  
  
  
  #Drop out predictors that you aren't going to use
  if(taxa == 'birds'){
    rastlist <- rastlist[!str_detect(rastlist,pattern="arm_diversity")]
    rastlist <- rastlist[!str_detect(rastlist,pattern="abirdDiversity_te")]
  }
  if(taxa == 'birds_te'){
    rastlist <- rastlist[!str_detect(rastlist,pattern="arm_diversity")]
    rastlist <- rastlist[!str_detect(rastlist,pattern="abirdDiversity.tif")]
  }
  if(taxa == 'arm'){
    rastlist <- rastlist[!str_detect(rastlist,pattern="abirdDiversity.tif")]
    rastlist <- rastlist[!str_detect(rastlist,pattern="abirdDiversity_te")]
  }
  
  
  
  #If doing just in PAs then don't use distance to PAs
  if(extent == 'PAs'){
    rastlist <- rastlist[!str_detect(rastlist,pattern="distPAs")]
  }
  
  #Make raster stack and do transformations 
  allrasters <- stack(rastlist)
  allrasters$distRoads <- sqrt(allrasters$distRoads )
  allrasters$distWater <- sqrt(allrasters$distWater )
  allrasters$hotelDensity <- log10(allrasters$hotelDensity + 0.1)

  #Only transform distance to PAs if needed
  if(extent != 'PAs'){
    allrasters$distPAs <- sqrt(allrasters$distPAs)
  }
  
  #Get names out 
  rename <- names(allrasters)
  
  
  #Projection systems
  wgsProj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  utmProj <- "+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs"
  
  
  #Read in CR shapefile and transform to UTM
  cr <- rgdal::readOGR(file.path(shapeDir, 'cr.shp'))
  proj4string(cr) <- CRS(wgsProj)
  cr <- spTransform(cr, utmProj)
  
  
  #Read in flickr data
  flickrP <- as.data.frame(fread(file.path(shapeDir, 'flickr.csv')))
  flickrP <- flickrP[c('decimalLongitude', 'decimalLatitude')]
  coordinates(flickrP)<-~decimalLongitude+decimalLatitude
  
  #Transform to utm
  proj4string(flickrP) <- CRS(wgsProj)
  flickrP <- spTransform(flickrP, utmProj)
  
  #Clip to CR
  flickrP <- flickrP[cr,]
  
  
  
  
  #Read in eBird data
  ebirdP <- as.data.frame(fread(file.path(shapeDir, "checklist_uploads.csv")))  
  
  #Clip years 
  ebirdP <- ebirdP[ebirdP$YEAR >= 2005, ]
  ebirdP <- ebirdP[ebirdP$YEAR <= 2017, ]
  
  #Reproject 
  ebirdP <- ebirdP[c('LONGITUDE', 'LATITUDE')]
  coordinates(ebirdP)<-~LONGITUDE+LATITUDE
  proj4string(ebirdP) <- CRS(wgsProj)
  ebirdP <- spTransform(ebirdP, utmProj)
  
  #Clip to CR
  ebirdP <- ebirdP[cr,]
  
  
  
  #Do stuff specific to PAs
  if(extent == 'PAs'){
    #Read in and reproject PAs
    sFile <- rgdal::readOGR(file.path(shapeDir, "CR_protected_areas\\WDPA_Dec2018_CRI-shapefile-polygons.shp"))
    sFile <- spTransform(sFile, utmProj)
    sFile <- raster::intersect(sFile, cr)
    
    #Calculate area 
    sFile$area_sqkm <- raster::area(sFile) / 1000000
    
    #Create masks 
    pa_mask <- rasterize(sFile, allrasters$distRoads, field = 1, update = F)
    pa_size <- rasterize(sFile, allrasters$distRoads, field = 'area_sqkm', update = F)
    #pa_size <- log(pa_size)
    #Used this to try dropping out amistad
    #pa_mask[pa_size > 2000] <- NA
    
    #Clip to masks 
    allrasters <- pa_mask*allrasters
    names(allrasters) <- rename
    
    #update raster stack and names 
    allrasters <- stack(allrasters, pa_size)
    names(allrasters)[5] <- 'pa_size'
    
    #Clip points to just protected areas 
    flickrP <- flickrP[sFile,]
    ebirdP <- ebirdP[sFile,]
  }
  
  #Normalize all varaibles between 0 and 1 (makes maxent run better) 
  allrasters <- (allrasters - minValue(allrasters)) / (maxValue(allrasters) - minValue(allrasters))
  
  #Specify out directory) 
  outDir <- file.path(baseOut, 'ebird')
  dir.create(outDir)
  
  #Run ebird model
  birdModel <- maxent(allrasters, ebirdP,  removeDuplicates=T, path = outDir,
                      args=c("-J", "-P", "replicates=5", 
                             'linear=true', 
                             'quadratic=true',
                             'product=false',
                             'hinge=false',
                             'threshold=false',
                             'replicatetype=Crossvalidate',
                             'writeplotdata=true',
                             'extrapolate=false'))
  
  
  #Create tourism raster using flickr model (mean of all runs)
  ebirdRaster <- mean(predict(birdModel, allrasters))
  
  #Write rasters
  writeRaster(ebirdRaster, file.path(outDir, 'predictedTourism.tif'), overwrite = T)
  writeRaster(allrasters,  file.path(outDir, 'predictors.tif'), overwrite = T)
  write.csv(names(allrasters), file.path(outDir, 'names.csv'))
  
  #Specify out directory 
  outDir <- file.path(baseOut, 'flickr')
  dir.create(outDir)
  
  #Run flickr model
  flickrModel <- maxent(allrasters, flickrP,  removeDuplicates=T, path = outDir,
                        args=c("-J", "-P", "replicates=5", 
                               'linear=true', 
                               'quadratic=true',
                               'product=false',
                               'hinge=false',
                               'threshold=false',
                               'replicatetype=Crossvalidate',
                               'writeplotdata=true',
                               'extrapolate=false'))
  

  #Create tourism raster using flickr model (mean of all runs)
  flickrRaster <- mean(predict(flickrModel, allrasters))
  
  #Write rasters
  writeRaster(flickrRaster, file.path(outDir, 'predictedTourism.tif'), overwrite = T)
  writeRaster(allrasters,  file.path(outDir, 'predictors.tif'), overwrite = T)
  write.csv(names(allrasters), file.path(outDir, 'names.csv'))
}

#Code to loop through models
runAllModels <- function(){
  #Run birds, national model
  runModels('birds', 'national', file.path(maxentOut, 'birds_national'))
  
  #Run birds PAs model
  runModels('birds', 'PAs',  file.path(maxentOut, 'birds_PAs'))
  
  #Run TE birds, national model
  runModels('birds_te', 'national',  file.path(maxentOut, 'birdsTE_national'))
  
  #Run TE birds, PA model 
  runModels('birds_te', 'PAs',  file.path(maxentOut, 'birdsTE_PAs'))
  
  #Run tetrapod, national model
  runModels('arm', 'national',  file.path(maxentOut, 'arm_national'))
  
  #Run tetrapod, PA model
  runModels('arm', 'PAs',  file.path(maxentOut, 'arm_PAs'))
}
runAllModels()
  

#Function to make SI graphs 
make_SI_scatter <- function(){
  #Base data you'll need for this function 
  #Projection systems
  wgsProj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  utmProj <- "+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs"
  
  
  #Read in CR shapefile and transform to UTM
  cr <- rgdal::readOGR(file.path(shapeDir, 'cr.shp'))
  cr <- spTransform(cr, utmProj)
  
  
  #Read in flickr data
  flickrP <- as.data.frame(fread(file.path(shapeDir, 'flickr.csv')))
  flickrP <- flickrP[c('decimalLongitude', 'decimalLatitude')]
  coordinates(flickrP)<-~decimalLongitude+decimalLatitude
  
  #Transform to utm
  proj4string(flickrP) <- CRS(wgsProj)
  flickrP <- spTransform(flickrP, utmProj)
  
  #Clip to CR
  flickrP <- flickrP[cr,]
  
  
  
  
  #Read in eBird data
  ebirdP <- as.data.frame(fread(file.path(shapeDir, "checklist_uploads.csv")))  
  
  #Clip years 
  ebirdP <- ebirdP[ebirdP$YEAR >= 2005, ]
  ebirdP <- ebirdP[ebirdP$YEAR <= 2017, ]
  
  #Reproject 
  ebirdP <- ebirdP[c('LONGITUDE', 'LATITUDE')]
  coordinates(ebirdP)<-~LONGITUDE+LATITUDE
  proj4string(ebirdP) <- CRS(wgsProj)
  ebirdP <- spTransform(ebirdP, utmProj)
  
  #Clip to CR
  ebirdP <- ebirdP[cr,]
  
  #Make the individual GGplot elements
  makeGraph <- function(varName, dirPath, xlabel, pointsToUse){
    
    #Location of the data 
    outDir <- file.path(maxentOut,dirPath, pointsToUse)
    allrasters <- stack(file.path(outDir, 'predictors.tif'))
    names(allrasters) <- read.csv(file.path(outDir, 'names.csv'), row.names = 1)$x
    
    
    #Get values of predicted tourism and predictor values for ebird or flickr 
    if(pointsToUse == 'ebird'){
      
      #Read in out rast
      outRast <- raster(file.path(outDir, 'predictedTourism.tif'))
      
      #Extract the env variable you are going to plot
      x <- extract(allrasters[[varName]] , ebirdP)
      
      #Extract predicted tourism
      y <- extract(outRast, ebirdP)
    }
    
    #if flickr 
    if(pointsToUse == 'flickr'){
      
      #Create tourism raster using ebird model (mean of all runs)
      outRast <- raster(file.path(outDir, 'predictedTourism.tif'))
      
      #Extract the env variable you are going to plot
      x <- extract(allrasters[[varName]] , flickrP)
      
      #Extract predicted tourism
      y <- extract(outRast, flickrP)
    }  
    
    #Make those a dataframe 
    toPlot2 <- as.data.frame(cbind(x, y))
    
    
    #Set up datafame to hold results
    holding <- c()
    
    #Loop through all replicated maxent runs 
    for(i in c(0:4)){
      
      #Read in the results from that replicate
      hold <- read.csv(file.path(outDir, 'plots', 
                                 paste('species_', as.character(i), '_', varName,'.dat', sep = '')))
      
      #Add the predicted tourism from that run to the holding DF 
      holding <- cbind(holding, hold$y)
      
      #Get out the x values 
      xVals <- hold$x
    }
    
    #Get the minimum values from all the runs 
    minValues <- rowMins(holding, value = T)
    
    #Get the maximum values from all the runs 
    maxValues <- rowMaxs(holding, value = T)
    
    #Make dataframe for stuff to plot
    toPlot <- as.data.frame(cbind(xVals, minValues, maxValues))
    names(toPlot) <- c('x', 'minValues', 'maxValues')
    
    #Get mean values 
    toPlot$meanValues <- rowMeans(holding)
    
    #Make plot 
    m1 <- ggplot(toPlot, aes(x = x, y = meanValues)) +
      geom_line() +
      xlim(0,1) +
      ylim(0.01,1) +
      theme_classic() + 
      geom_ribbon(data = toPlot, aes(ymin = minValues, ymax = maxValues), 
                  alpha = 0.1) +
      geom_jitter(data = toPlot2[sample(nrow(toPlot2), 5000), ], 
                  mapping = aes(x = x, y = y),
                  alpha = 0.1) +
      xlab(xlabel) + 
      ylab(paste('Tourism (', pointsToUse, ')', sep = ''))
    
    
    return(m1)
  }
  
  #Make all the graphs for national scale 
  a1 <- makeGraph('abirdDiversity', 'birds_national', 'Bird richness\n', 'ebird')
  a2 <- makeGraph('abirdDiversity_te', 'birdsTE_national', 'Threatened and endemic\nbird richness', 'ebird')
  a3 <- makeGraph('abirdDiversity', 'birds_national', 'Bird richness\n', 'flickr')
  a4 <- makeGraph('abirdDiversity_te', 'birdsTE_national', 'Threatened and endemic\nbird richness', 'flickr')
  a5 <- makeGraph('arm_diversity', 'arm_national', 'Amphibian, reptile, and\nmammal richness', 'flickr')
  
  a6 <- makeGraph('distRoads', 'birds_national', 'Distance to roads\n(meters^0.5)', 'ebird')
  a7 <- makeGraph('distRoads', 'birdsTE_national', 'Distance to roads\n(meters^0.5)', 'ebird')
  a8 <- makeGraph('distRoads', 'birds_national', 'Distance to roads\n(meters^0.5)', 'flickr')
  a9 <- makeGraph('distRoads', 'birdsTE_national', 'Distance to roads\n(meters^0.5)', 'flickr')
  a10 <- makeGraph('distRoads', 'arm_national', 'Distance to roads\n(meters^0.5)', 'flickr')
  
  a11 <- makeGraph('distWater', 'birds_national', 'Distance to water (sqrt)', 'ebird')
  a12 <- makeGraph('distWater', 'birdsTE_national', 'Distance to water\n(meters^0.5)', 'ebird')
  a13 <- makeGraph('distWater', 'birds_national', 'Distance to water\n(meters^0.5)', 'flickr')
  a14 <- makeGraph('distWater', 'birdsTE_national', 'Distance to water\n(meters^0.5)', 'flickr')
  a15 <- makeGraph('distWater', 'arm_national', 'Distance to w0ter\n(meters^0.5)', 'flickr')
  
  a16 <- makeGraph('distPAs', 'birds_national', 'Distance to protected areas\n(meters^0.5)', 'ebird')
  a17 <- makeGraph('distPAs', 'birdsTE_national', 'Distance to protected areas\n(meters^0.5)', 'ebird')
  a18 <- makeGraph('distPAs', 'birds_national', 'Distance to protected areas\n(meters^0.5)', 'flickr')
  a19 <- makeGraph('distPAs', 'birdsTE_national', 'Distance to protected areas\n(meters^0.5)', 'flickr')
  a20 <- makeGraph('distPAs', 'arm_national', 'Distance to protected areas\n(meters^0.5)', 'flickr')
  
  a21 <- makeGraph('hotelDensity', 'birds_national', 'Hotel density\n(log)', 'ebird')
  a22 <- makeGraph('hotelDensity', 'birdsTE_national', 'Hotel density\n(log)', 'ebird')
  a23 <- makeGraph('hotelDensity', 'birds_national', 'Hotel density\n(log)', 'flickr')
  a24 <- makeGraph('hotelDensity', 'birdsTE_national', 'Hotel density\n(log)', 'flickr')
  a25 <- makeGraph('hotelDensity', 'arm_national', 'Hotel density\n(log)', 'flickr')
  
  #Save to a fig
  outFile <- file.path(maxentOut, "national_full.png")
  png(outFile, width = 15, height = 10, units = 'in', res = 600)  
  multiplot(a1, a2, a3, a4, a5,
            a6, a7, a8, a9, a10,
            a11, a12, a13, a14, a15,
            a16, a17, a18, a19, a20,
            a21, a22, a23, a24, a25,
            cols=5)
  dev.off()
  
  #Make all the graphs for protected area
  a1 <- makeGraph('abirdDiversity', 'birds_PAs', 'Bird richness\n', 'ebird')
  a2 <- makeGraph('abirdDiversity_te', 'birdsTE_PAs', 'Threatened and endemic\nbird richness', 'ebird')
  a3 <- makeGraph('abirdDiversity', 'birds_PAs', 'Bird richness\n', 'flickr')
  a4 <- makeGraph('abirdDiversity_te', 'birdsTE_PAs', 'Threatened and endemic\nbird richness', 'flickr')
  a5 <- makeGraph('arm_diversity', 'arm_PAs', 'Amphibian, reptile, and\nmammal richness', 'flickr')
  
  a6 <- makeGraph('distRoads', 'birds_PAs', 'Distance to roads\n(meters^0.5)', 'ebird')
  a7 <- makeGraph('distRoads', 'birdsTE_PAs', 'Distance to roads\n(meters^0.5)', 'ebird')
  a8 <- makeGraph('distRoads', 'birds_PAs', 'Distance to roads\n(meters^0.5)', 'flickr')
  a9 <- makeGraph('distRoads', 'birdsTE_PAs', 'Distance to roads\n(meters^0.5)', 'flickr')
  a10 <- makeGraph('distRoads', 'arm_PAs', 'Distance to roads\n(meters^0.5)', 'flickr')
  
  a11 <- makeGraph('distWater', 'birds_PAs', 'Distance to water\n(meters^0.5)', 'ebird')
  a12 <- makeGraph('distWater', 'birdsTE_PAs', 'Distance to water\n(meters^0.5)', 'ebird')
  a13 <- makeGraph('distWater', 'birds_PAs', 'Distance to water\n(meters^0.5)', 'flickr')
  a14 <- makeGraph('distWater', 'birdsTE_PAs', 'Distance to water\n(meters^0.5)', 'flickr')
  a15 <- makeGraph('distWater', 'arm_PAs', 'Distance to water\n(meters^0.5)', 'flickr')
  
  a16 <- makeGraph('pa_size', 'birds_PAs', 'Protected area size\n(sq. km)', 'ebird')
  a17 <- makeGraph('pa_size', 'birdsTE_PAs', 'Protected area size\n(sq. km)', 'ebird')
  a18 <- makeGraph('pa_size', 'birds_PAs', 'Protected area size\n(sq. km)', 'flickr')
  a19 <- makeGraph('pa_size', 'birdsTE_PAs', 'Protected area size\n(sq. km)', 'flickr')
  a20 <- makeGraph('pa_size', 'arm_PAs', 'Protected area size\n(sq. km)', 'flickr')
  
  a21 <- makeGraph('hotelDensity', 'birds_PAs', 'Hotel density\n(log)', 'ebird')
  a22 <- makeGraph('hotelDensity', 'birdsTE_PAs', 'Hotel density\n(log)', 'ebird')
  a23 <- makeGraph('hotelDensity', 'birds_PAs', 'Hotel density\n(log)', 'flickr')
  a24 <- makeGraph('hotelDensity', 'birdsTE_PAs', 'Hotel density\n(log)', 'flickr')
  a25 <- makeGraph('hotelDensity', 'arm_PAs', 'Hotel density\n(log)', 'flickr')
  
  
  #Save to a fig
  outFile <- file.path(maxentOut, "PA_full.png")
  png(outFile, width = 15, height = 10, units = 'in', res = 600)  
  multiplot(a1, a2, a3, a4, a5,
            a6, a7, a8, a9, a10,
            a11, a12, a13, a14, a15,
            a16, a17, a18, a19, a20,
            a21, a22, a23, a24, a25,
            cols=5)
  dev.off()
}

#Run the code to make the scatter plots 
make_SI_scatter()


#Function to make main graphs 
make_main_lineGraphs <- function(){
  
  #Data you'll need
  #Projection systems
  wgsProj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  utmProj <- "+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs"
  
  
  #Read in CR shapefile and transform to UTM
  cr <- rgdal::readOGR(file.path(shapeDir, 'cr.shp'))
  proj4string(cr) <- CRS(wgsProj)
  cr <- spTransform(cr, utmProj)
  
  
  #Read in flickr data
  flickrP <- as.data.frame(fread(file.path(shapeDir, 'flickr.csv')))
  flickrP <- flickrP[c('decimalLongitude', 'decimalLatitude')]
  coordinates(flickrP)<-~decimalLongitude+decimalLatitude
  
  #Transform to utm
  proj4string(flickrP) <- CRS(wgsProj)
  flickrP <- spTransform(flickrP, utmProj)
  
  #Clip to CR
  flickrP <- flickrP[cr,]
  
  
  
  
  #Read in eBird data
  ebirdP <- as.data.frame(fread(file.path(shapeDir, "checklist_uploads.csv")))  
  
  #Clip years 
  ebirdP <- ebirdP[ebirdP$YEAR >= 2005, ]
  ebirdP <- ebirdP[ebirdP$YEAR <= 2017, ]
  
  #Reproject 
  ebirdP <- ebirdP[c('LONGITUDE', 'LATITUDE')]
  coordinates(ebirdP)<-~LONGITUDE+LATITUDE
  proj4string(ebirdP) <- CRS(wgsProj)
  ebirdP <- spTransform(ebirdP, utmProj)
  
  #Clip to CR
  ebirdP <- ebirdP[cr,]
  
  #Read in and reproject PAs
  sFile <- rgdal::readOGR(file.path(shapeDir, "CR_protected_areas\\WDPA_Dec2018_CRI-shapefile-polygons.shp"))
  sFile <- spTransform(sFile, utmProj)
  sFile <- raster::intersect(sFile, cr)
  
  #Calculate area 
  sFile$area_sqkm <- raster::area(sFile) / 1000000
  
  #Return single panels 
  makeGraph2 <- function(extentZ, varName, xlabel, biodiv){
    
    
    #Function to exract data 
    getTables <- function(pointsToUse, dirPath){
      
      #Location of the data 
      outDir <- file.path(maxentOut,dirPath, pointsToUse)
      allrasters <- stack(file.path(outDir, 'predictors.tif'))
      names(allrasters) <- read.csv(file.path(outDir, 'names.csv'), row.names = 1)$x
      
      
      #Get values of predicted tourism and predictor values for ebird or flickr 
      if(pointsToUse == 'ebird'){
        
        #Read in out rast
        outRast <- raster(file.path(outDir, 'predictedTourism.tif'))
        
        #Extract the env variable you are going to plot
        x <- extract(allrasters[[varName]] , ebirdP)
        
        #Extract predicted tourism
        y <- extract(outRast, ebirdP)
      }
      
      #if flickr 
      if(pointsToUse == 'flickr'){
        
        #Create tourism raster using ebird model (mean of all runs)
        outRast <- raster(file.path(outDir, 'predictedTourism.tif'))
        
        #Extract the env variable you are going to plot
        x <- extract(allrasters[[varName]] , flickrP)
        
        #Extract predicted tourism
        y <- extract(outRast, flickrP)
      } 
      
      toPlot2 <- as.data.frame(cbind(x, y))
      
      
      #Define all the predictors you're going to use
      rastlist <- list.files(path = rasterDir, pattern='tif$', all.files=TRUE, full.names=TRUE)
      rastlist <- rastlist[!str_detect(rastlist,pattern="cover")]
      rastlist <- rastlist[!str_detect(rastlist,pattern="bio")]
      rastlist <- rastlist[!str_detect(rastlist,pattern="amphibian")]
      rastlist <- rastlist[!str_detect(rastlist,pattern="mammal")]
      rastlist <- rastlist[!str_detect(rastlist,pattern="reptile")]
      
      allrasters <- stack(rastlist)
      allrasters$distRoads <- sqrt(allrasters$distRoads )
      allrasters$distWater <- sqrt(allrasters$distWater )
      allrasters$hotelDensity <- log10(allrasters$hotelDensity + 0.1)
      
      
      pa_mask <- rasterize(sFile, allrasters$distRoads, field = 1, update = F)
      pa_size <- rasterize(sFile, allrasters$distRoads, field = 'area_sqkm', update = F)
      allrasters$pa_size <- pa_size
      
      
      #Add replicate runs to a ddat frame 
      holding <- c()
      for(i in c(0:4)){
        hold <- read.csv(file.path(outDir, 'plots', 
                                   paste('species_', as.character(i), '_', varName,'.dat', sep = '')))
        holding <- cbind(holding, hold$y)
        xVals <- hold$x
        
        #back-transform axes (except hotel density)
        if(varName != 'hotelDensity'){
          xVals <- xVals * (cellStats(allrasters[[varName]], max) - cellStats(allrasters[[varName]], min)) + cellStats(allrasters[[varName]], min)
          
        }
      }
      
      #Create plotting DF
      minValues <- rowMins(holding, value = T)
      maxValues <- rowMaxs(holding, value = T)
      toPlot <- as.data.frame(cbind(xVals, minValues, maxValues))
      names(toPlot) <- c('x', 'minValues', 'maxValues')
      toPlot$meanValues <- rowMeans(holding)
      
      #Reutrn the dataframe to make points and lines 
      toReturn <- c(toPlot, toPlot2)
    }
    
    #If it is not a biodiveristy variable get the data
    if(biodiv == 'no'){
      if(extentZ == 'national'){
        t1 <- getTables('ebird', 'birds_national')
        t2 <- getTables('ebird', 'birdsTE_national')
        t3 <- getTables('flickr', 'birds_national')
        t4 <- getTables('flickr', 'birdsTE_national')
        t5 <- getTables('flickr', 'arm_national')
      }
      if(extentZ=='PAs'){
        t1 <- getTables('ebird', 'birds_PAs')
        t2 <- getTables('ebird', 'birdsTE_PAs')
        t3 <- getTables('flickr', 'birds_PAs')
        t4 <- getTables('flickr', 'birdsTE_PAs')
        t5 <- getTables('flickr', 'arm_PAs')
      }
      
      #Split the dataframes into the dataframes to make lines (a) and scatter (b)
      t1a <- as.data.frame(t1[1:4])  
      t1b <- as.data.frame(t1[5:6])
      
      
      t2a <- as.data.frame(t2[1:4])  
      t2b <- as.data.frame(t2[5:6]) 
      
      
      t3a <- as.data.frame(t3[1:4])  
      t3b <- as.data.frame(t3[5:6]) 
      
      
      t4a <- as.data.frame(t4[1:4])  
      t4b <- as.data.frame(t4[5:6]) 
      
      
      t5a <- as.data.frame(t5[1:4])  
      t5b <- as.data.frame(t5[5:6])
      
      
      
      
      #Make graphs 
      m1 <- ggplot(t1a, aes(x = x, y = meanValues)) +
        geom_line(col = c1, linetype='dashed', size = 0.25) +
        ylim(0.01,1) +
        theme_classic() + 
        geom_ribbon(data = t1a, aes(ymin = minValues, ymax = maxValues), 
                    alpha = 0.5, fill = c1) +
        xlab(xlabel) + 
        ylab('Tourism') +
        
        geom_line(data = t2a, (aes(x=x, y =meanValues)), color = c2, linetype='dashed', size = 0.25) +
        geom_ribbon(data = t2a, aes(ymin = minValues, ymax = maxValues), 
                    alpha = 0.5, fill = c2) +
        geom_line(data = t3a, (aes(x=x, y =meanValues)), color = c3, size = 0.25) +
        geom_ribbon(data = t3a, aes(ymin = minValues, ymax = maxValues), 
                    alpha = 0.5, fill = c3) +
        geom_line(data = t4a, (aes(x=x, y =meanValues)), color = c4,  size = 0.25) +
        geom_ribbon(data = t4a, aes(ymin = minValues, ymax = maxValues), 
                    alpha = 0.1, col = c4, fill = c4) +
        geom_line(data = t5a, (aes(x=x, y =meanValues)), color = c5,  size = 0.25) +
        geom_ribbon(data = t5a, aes(ymin = minValues, ymax = maxValues), 
                    alpha = 0.5, fill = c5) 
    }
    
    #Do the same for biodiversity variables 
    if(biodiv == 'birds'){
      if(extentZ == 'national'){
        t1 <- getTables('ebird', 'birds_national')
        t2 <- getTables('flickr', 'birds_national')
        
      }
      if(extentZ=='PAs'){
        t1 <- getTables('ebird', 'birds_PAs')
        t2 <- getTables('flickr', 'birds_PAs')
        
      }
      t1a <- as.data.frame(t1[1:4])  
      t1b <- as.data.frame(t1[5:6])
      
      t2a <- as.data.frame(t2[1:4])  
      t2b <- as.data.frame(t2[5:6]) 
      
      m1 <- ggplot(t1a, aes(x = x, y = meanValues)) +
        geom_line(col = c1, linetype='dashed', size = 0.25) +
        ylim(0.01,1) +
        theme_classic() + 
        geom_ribbon(data = t1a, aes(ymin = minValues, ymax = maxValues), 
                    alpha = 0.5, fill = c1) +
        xlab(xlabel) + 
        ylab('Tourism') +
        
        geom_line(data = t2a, (aes(x=x, y =meanValues)), color = c3, size = 0.25) +
        geom_ribbon(data = t2a, aes(ymin = minValues, ymax = maxValues), 
                    alpha = 0.5, fill = c3) 
      
      
      
    }
    
    if(biodiv == 'birdsTE'){
      if(extentZ == 'national'){
        t1 <- getTables('ebird', 'birdsTE_national')
        t2 <- getTables('flickr', 'birdsTE_national')
        
      }
      if(extentZ=='PAs'){
        t1 <- getTables('ebird', 'birdsTE_PAs')
        t2 <- getTables('flickr', 'birdsTE_PAs')
        
      }
      t1a <- as.data.frame(t1[1:4])  
      t1b <- as.data.frame(t1[5:6])
      
      t2a <- as.data.frame(t2[1:4])  
      t2b <- as.data.frame(t2[5:6]) 
      
      m1 <- ggplot(t1a, aes(x = x, y = meanValues)) +
        geom_line(col = c2, linetype='dashed', size = 0.25) +
        ylim(0.01,1) +
        theme_classic() + 
        geom_ribbon(data = t1a, aes(ymin = minValues, ymax = maxValues), 
                    alpha = 0.5, fill = c2) +
        xlab(xlabel) + 
        ylab('Tourism') +
        
        geom_line(data = t2a, (aes(x=x, y =meanValues)), color = c4, size = 0.25) +
        geom_ribbon(data = t2a, aes(ymin = minValues, ymax = maxValues), 
                    alpha = 0.5, fill = c4) 
      
    }  
    
    if(biodiv == 'arm'){
      if(extentZ == 'national'){
        t1 <- getTables('flickr', 'arm_national')
        
      }
      if(extentZ=='PAs'){
        t1 <- getTables('flickr', 'arm_PAs')
        
      }
      t1a <- as.data.frame(t1[1:4])  
      t1b <- as.data.frame(t1[5:6])
      
      
      
      m1 <- ggplot(t1a, aes(x = x, y = meanValues)) +
        geom_line(col = c5, size = 0.25) +
        ylim(0.01,1) +
        theme_classic() + 
        geom_ribbon(data = t1a, aes(ymin = minValues, ymax = maxValues), 
                    alpha = 0.5, fill = c5) +
        
        xlab(xlabel) + 
        ylab('Tourism') 
      
    }   
    
    
    
    
    
    
    return(m1)
  }
  
  #Make legend 
  c1 <- '#b2df8a'
  c2 <- '#a6cee3'
  c3 <- '#33a02c'
  c4 <- '#1f78b4'
  c5 <- '#ff91a4'
  legendFile <- file.path(maxentOut, "legend.png")
  png(legendFile, width = 1800, height = 2600)
  plot(0,0,pch=19,col="white",xlim=c(0,4), ylim=c(0.5,5.5),
       frame.plot=F, axes=F, xlab="",ylab="", 
       main="")
  lineWidth = 30
  segments(0,5,1,5, col = c1,  lty = 5, lwd = lineWidth)
  segments(0,4,1,4, col = c2,  lty = 5, lwd = lineWidth)
  segments(0,3,1,3, col = c3,  lty = 1, lwd = lineWidth)
  segments(0,2,1,2, col = c4,  lty = 1, lwd = lineWidth)
  segments(0,1,1,1, col = c5,  lty = 1, lwd = lineWidth)
  
  textSize <- 10
  text(1.25,5, '\nBirds\nResponse = eBird\n', adj = 0, cex = textSize)
  text(1.25,4, 'Threatened and\nendemic birds\nResponse = eBird', adj = 0, cex = textSize)
  text(1.25,3, 'Birds\nResponse = Flickr', adj = 0, cex = textSize)
  text(1.25,2, 'Threatened and\nendemic birds\nResponse = Flickr', adj = 0, cex = textSize)
  text(1.25,1, 'Amphibians, reptiles\nand mammals\nResponse = Flickr', adj = 0, cex =textSize)
  dev.off()
  
  
  #Make panels for non biodiveristy variables 
  a1 <- makeGraph2('national', 'distRoads', 'Distance to roads\n(meters^0.5)', 'no')
  a2 <- makeGraph2('national', 'distPAs', 'Distance to protected areas\n(meters^0.5)', 'no')
  a3 <- makeGraph2('national', 'distWater', 'Distance to water\n(meters^0.5)', 'no')
  a4 <- makeGraph2('national', 'hotelDensity', 'Hotel density\n(normalized on log scale)', 'no')
  
  #Make biodiveristy panels 
  a5 <- makeGraph2('national', 'abirdDiversity', 'Bird richness\n', 'birds')
  a6 <- makeGraph2('national', 'abirdDiversity_te', 'Threatened and\nendemic bird richness', 'birdsTE')
  a7 <- makeGraph2('national', 'arm_diversity', 'Amphibian, reptile,\nand mammal richness', 'arm')
  
  #Add text
  a1 <- a1 + annotate('text', x=0, y=1 , label='D.')
  a2 <- a2 + annotate('text', x=0, y=1 , label='E.')
  a3 <- a3 + annotate('text', x=0, y=1 , label='F.')
  a4 <- a4 + annotate('text', x=0, y=1 , label='G.')
  a5 <- a5 + annotate('text', x=25, y=1 , label='A.')
  a6 <- a6 + annotate('text', x=15, y=1 , label='B.')
  a7 <- a7 + annotate('text', x=10, y=1 , label='C.')
  
  #Read in legend 
  img <- readPNG(legendFile, native = TRUE)
  a0 <- ggplot() + 
    xlim(0,1) + 
    ylim(0,1) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()) +
    inset_element(img, 0,0,1,1) 
  
  #Plot and save fig        
  outFile <- file.path(maxentOut, "national.png")
  png(outFile, width = 10, height = 6, units = 'in', res = 600)  
  multiplot(a0, a1, a5, a2, a6, a3, a7, a4, cols = 4)
  dev.off()
  
  #Get non-bio panels
  a1a <- makeGraph2('PAs', 'distRoads', 'Distance to roads\n(meters^0.5)', 'no')
  a2a <- makeGraph2('PAs', 'pa_size', 'Protected area size\n(square km)', 'no')
  a3a <- makeGraph2('PAs', 'distWater', 'Distance to water\n(meters^0.5)', 'no')
  a4a <- makeGraph2('PAs', 'hotelDensity', 'Hotel density\n(normalized on log scale)', 'no')
  
  #Get bio panels
  a5a <- makeGraph2('PAs', 'abirdDiversity', 'Bird richness\n', 'birds')
  a6a <- makeGraph2('PAs', 'abirdDiversity_te', 'Threatened and\nendemic bird richness', 'birdsTE')
  a7a <- makeGraph2('PAs', 'arm_diversity', 'Amphibian, reptile,\nand mammal richness', 'arm')
  
  #Add text
  a1a <- a1a + annotate('text', x=0, y=1 , label='D.')
  a2a <- a2a + annotate('text', x=0, y=1 , label='E.')
  a3a <- a3a + annotate('text', x=0, y=1 , label='F.')
  a4a <- a4a + annotate('text', x=0, y=1 , label='G.')
  a5a <- a5a + annotate('text', x=25, y=1 , label='A.')
  a6a <- a6a + annotate('text', x=15, y=1 , label='B.')
  a7a <- a7a + annotate('text', x=10, y=1 , label='C.')
  
  #Plot and save figure
  outFile <- file.path(maxentOut, "PAs.png")
  png(outFile, width = 10, height = 6, units = 'in', res = 600)  
  multiplot(a0, a1a, a5a, a2a, a6a, a3a, a7a, a4a, cols = 4)
  dev.off()
  
}


#Run code to make main figures 
make_main_lineGraphs()




#Function to make result maps
makeResultMaps <- function(){
  #Projection systems
  wgsProj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  utmProj <- "+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs"
  
  
  #Read in CR shapefile and transform to UTM
  cr <- rgdal::readOGR(file.path(shapeDir, 'cr.shp'))
  cr <- spTransform(cr, utmProj)
  
  
  #Read in flickr data
  flickrP <- as.data.frame(fread(file.path(shapeDir, 'flickr.csv')))
  flickrP <- flickrP[c('decimalLongitude', 'decimalLatitude')]
  coordinates(flickrP)<-~decimalLongitude+decimalLatitude
  
  #Transform to utm
  proj4string(flickrP) <- CRS(wgsProj)
  flickrP <- spTransform(flickrP, utmProj)
  
  #Clip to CR
  flickrP <- flickrP[cr,]
  
  
  
  
  #Read in eBird data
  ebirdP <- as.data.frame(fread(file.path(shapeDir, "checklist_uploads.csv")))
  
  #Clip years 
  ebirdP <- ebirdP[ebirdP$YEAR >= 2005, ]
  ebirdP <- ebirdP[ebirdP$YEAR <= 2017, ]
  
  #Reproject 
  ebirdP <- ebirdP[c('LONGITUDE', 'LATITUDE')]
  coordinates(ebirdP)<-~LONGITUDE+LATITUDE
  proj4string(ebirdP) <- CRS(wgsProj)
  ebirdP <- spTransform(ebirdP, utmProj)
  
  #Clip to CR
  ebirdP <- ebirdP[cr,]
  
  
  #Read in and reproject PAs
  sFile <- rgdal::readOGR(file.path(shapeDir, "CR_protected_areas\\WDPA_Dec2018_CRI-shapefile-polygons.shp"))
  sFile <- spTransform(sFile, utmProj)
  sFile <- raster::intersect(sFile, cr)
  
  
  #Define all the predictors you're going to use
  rastlist <- list.files(path = rasterDir, pattern='tif$', all.files=TRUE, full.names=TRUE)
  rastlist <- rastlist[!str_detect(rastlist,pattern="cover")]
  rastlist <- rastlist[!str_detect(rastlist,pattern="bio")]
  rastlist <- rastlist[!str_detect(rastlist,pattern="amphibian")]
  rastlist <- rastlist[!str_detect(rastlist,pattern="mammal")]
  rastlist <- rastlist[!str_detect(rastlist,pattern="reptile")]
  rastlist <- rastlist[!str_detect(rastlist,pattern="abirdDiversity.tif")]
  rastlist <- rastlist[!str_detect(rastlist,pattern="arm_diversity")]
  
  allrasters <- stack(rastlist)
  allrasters$distRoads <- sqrt(allrasters$distRoads )
  allrasters$distWater <- sqrt(allrasters$distWater )
  allrasters$hotelDensity <- log10(allrasters$hotelDensity + 0.1)
  
  
  allrasters <- (allrasters - minValue(allrasters)) / (maxValue(allrasters) - minValue(allrasters))
  
  #Make basemap for ebird
  birdModel <- maxent(allrasters, ebirdP,  removeDuplicates=T,
                      args=c("-J", "-P", "replicates=5", 
                             'linear=true', 
                             'quadratic=true',
                             'product=false',
                             'hinge=false',
                             'threshold=false',
                             'replicatetype=Crossvalidate',
                             'writeplotdata=true',
                             'extrapolate=false'))
  
  baseMap <- mean(predict(birdModel, allrasters))
  
  #Make bird map with no biodiversity
  birdModel_noBio <- maxent(allrasters, ebirdP,  removeDuplicates=T, 
                            args=c("-J", "-P", "replicates=5", 
                                   'linear=true', 
                                   'quadratic=true',
                                   'product=false',
                                   'hinge=false',
                                   'threshold=false',
                                   'replicatetype=Crossvalidate',
                                   'writeplotdata=true',
                                   'extrapolate=false',
                                   'togglelayerselected=abirdDiversity_te'))
  noBio <- mean(predict(birdModel_noBio, allrasters))
  
  #Make ebird map with no infrastructure
  birdModel_noInfra <- maxent(allrasters, ebirdP,  removeDuplicates=T, 
                              args=c("-J", "-P", "replicates=5", 
                                     'linear=true', 
                                     'quadratic=true',
                                     'product=false',
                                     'hinge=false',
                                     'threshold=false',
                                     'replicatetype=Crossvalidate',
                                     'writeplotdata=true',
                                     'extrapolate=false',
                                     'togglelayerselected=distRoads',
                                     'togglelayerselected=hotelDensity'))
  noInfra <- mean(predict(birdModel_noInfra, allrasters))
  
  #Make flickr basemap
  birdModel_flickr <- maxent(allrasters, flickrP,  removeDuplicates=T,
                             args=c("-J", "-P", "replicates=5", 
                                    'linear=true', 
                                    'quadratic=true',
                                    'product=false',
                                    'hinge=false',
                                    'threshold=false',
                                    'replicatetype=Crossvalidate',
                                    'writeplotdata=true',
                                    'extrapolate=false'))
  
  baseMap_flickr <- mean(predict(birdModel_flickr, allrasters))
  
  #Make flickr map with no biodiveristy 
  birdModel_noBio_flickr <- maxent(allrasters, flickrP,  removeDuplicates=T, 
                                   args=c("-J", "-P", "replicates=5", 
                                          'linear=true', 
                                          'quadratic=true',
                                          'product=false',
                                          'hinge=false',
                                          'threshold=false',
                                          'replicatetype=Crossvalidate',
                                          'writeplotdata=true',
                                          'extrapolate=false',
                                          'togglelayerselected=abirdDiversity_te'))
  noBio_flickr <- mean(predict(birdModel_noBio_flickr, allrasters))
  
  #make flickr map with no infra
  birdModel_noInfra_flickr <- maxent(allrasters, flickrP,  removeDuplicates=T, 
                                     args=c("-J", "-P", "replicates=5", 
                                            'linear=true', 
                                            'quadratic=true',
                                            'product=false',
                                            'hinge=false',
                                            'threshold=false',
                                            'replicatetype=Crossvalidate',
                                            'writeplotdata=true',
                                            'extrapolate=false',
                                            'togglelayerselected=distRoads',
                                            'togglelayerselected=hotelDensity'))
  noInfra_flickr <- mean(predict(birdModel_noInfra_flickr, allrasters))
  
  mapDir <- file.path(maxentOut, 'mapTiffs')
  dir.create(mapDir)
  writeRaster(baseMap, file.path(mapDir, 'ebird_base.tif'), overwrite=TRUE)
  writeRaster(noBio, file.path(mapDir, 'ebird_noBio.tif'), overwrite=TRUE)
  writeRaster(noInfra, file.path(mapDir, 'ebird_noInfra.tif'), overwrite=TRUE)
  writeRaster(baseMap_flickr, file.path(mapDir, 'flickr_base.tif'), overwrite=TRUE)
  writeRaster(noBio_flickr, file.path(mapDir, 'flickr_noBio.tif'), overwrite=TRUE)
  writeRaster(noInfra_flickr, file.path(mapDir, 'flickr_noInfra.tif'), overwrite=TRUE)
}

#Make rasters that you will use to make graph
makeResultMaps()


#Function to make map figure in main paper
makeMapFigure <- function(){
  
  #Projection systems
  wgsProj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  utmProj <- "+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs"
  
  
  #Read in CR shapefile and transform to UTM
  cr <- rgdal::readOGR(file.path(shapeDir, 'cr.shp'))
  cr <- spTransform(cr, utmProj)
  
  #Read in and reproject PAs
  sFile <- rgdal::readOGR(file.path(shapeDir, "CR_protected_areas\\WDPA_Dec2018_CRI-shapefile-polygons.shp"))
  sFile <- spTransform(sFile, utmProj)
  sFile <- raster::intersect(sFile, cr)
  
  
  #Map directory 
  mapDir <- file.path(maxentOut, 'mapTiffs')
  
  #Read in maps 
  baseMap <- raster(file.path(mapDir, 'ebird_base.tif'))
  noBio <- raster(file.path(mapDir, 'ebird_noBio.tif'))
  noInfra <- raster(file.path(mapDir, 'ebird_noInfra.tif'))
  baseMap_flickr <- raster(file.path(mapDir, 'flickr_base.tif'))
  noBio_flickr <- raster(file.path(mapDir, 'flickr_noBio.tif'))
  noInfra_flickr <- raster(file.path(mapDir, 'flickr_noInfra.tif'))
  
  #Make figure
  outFile <- file.path(maxentOut, "maps.png")
  png(outFile, width = 10, height = 6, units = 'in', res = 600)  
  
  
  par(mfrow=c(2,3),
      mar=c(1,1,1,8))
  par(xpd = TRUE)
  corners = par("usr")
  cuts1 <- seq(0,1,0.05)
  cuts2 <- seq(-0.75,0.75,0.05)
  
  pal1 <- colorRampPalette(c('#f7f7f7', "#006d2c"))
  pal2 <- colorRampPalette(c("#ca0020","#f4a582", "#f7f7f7", "#92c5de", "#0571b0"))
  parkColor <- 'darkgoldenrod4'
  
  sFile2 <- spTransform(sFile, utmProj)
  
  plot(baseMap, breaks = cuts1, col = pal1(length(cuts1)), 
       axes = F,  box = F,  main = 'A. Predicted tourism (eBird)',
       legend.args=list(text='Higher tourism', side=3,  line=1, cex =  0.6))
  plot(baseMap, breaks = cuts1, col = pal1(length(cuts1)), 
       axes = F,  box = F, main = 'A. Predicted tourism (eBird)',
       legend.args=list(text='Lower tourism', side=1,  line=1, cex =  0.6),
       add = T)
  plot(sFile2, add = T, border = parkColor, lty = 3)
  
  plot(noBio-baseMap, breaks = cuts2, col = pal2(length(cuts1)), 
       axes = F,  box = F,  main = 'B. Change if bird richness is excluded',
       legend.args=list(text='Increased tourism', side=3,  line=1, cex =  0.6))
  plot(noBio-baseMap, breaks = cuts2, col = pal2(length(cuts2)), 
       axes = F,  box = F, main = 'B. Change if bird richness is excluded',
       legend.args=list(text='Decreased tourism', side=1,  line=1, cex =  0.6),
       add = T)
  plot(sFile2, add = T, border = parkColor, lty = 3)
  
  plot(noInfra-baseMap, breaks = cuts2, col = pal2(length(cuts2)), 
       axes = F,  box = F,main = 'C. Change if infrastructure is excluded',
       legend.args=list(text='Increased tourism', side=3,  line=1, cex =  0.6))
  plot(noInfra-baseMap, breaks = cuts2, col = pal2(length(cuts2)), 
       axes = F,  box = F, main = 'C. Change if infrastructure is excluded',
       legend.args=list(text='Decreased tourism', side=1,  line=1, cex =  0.6),
       add = T)
  plot(sFile2, add = T, border = parkColor, lty = 3)
  
  
  plot(baseMap_flickr, breaks = cuts1, col = pal1(length(cuts1)), 
       axes = F,  box = F, main = 'D. Predicted tourism (Flickr)',
       legend.args=list(text='Higher tourism', side=3,  line=1, cex =  0.6))
  plot(baseMap_flickr, breaks = cuts1, col = pal1(length(cuts1)), 
       axes = F,  box = F, main = 'D. Predicted tourism (Flickr)',
       legend.args=list(text='Lower tourism', side=1,  line=1, cex =  0.6),
       add = T)
  plot(sFile2, add = T, border = parkColor, lty = 3)
  north.arrow(xb=bbox(baseMap)['s1', 'min'] *1.1, yb=bbox(baseMap)['s2', 'min'] *1.1, len=15000, lab="N")  
  map.scale(x=bbox(baseMap)['s1', 'min'] *1.05, y=bbox(baseMap)['s2', 'min'] *1.05, ratio=FALSE, relwidth=0.4, cex=0.8)  
  rect(bbox(baseMap)['s1', 'min'] *1.18,
       bbox(baseMap)['s2', 'min'] *1.1,
       bbox(baseMap)['s1', 'min'] *1.29,
       bbox(baseMap)['s2', 'min'] *1.13,
       border = parkColor, lty = 3)
  text(bbox(baseMap)['s1', 'min'] *1.24,
       bbox(baseMap)['s2', 'min'] *1.085,
       'Protected areas')
  text(x = corners[2], y = mean(corners[3:4])*1.1, "Strength", srt = 0)
  
  
  plot(noBio_flickr-baseMap_flickr, breaks = cuts2, col = pal2(length(cuts2)), 
       axes = F,  box = F, main = 'E. Change if bird richness is excluded',
       legend.args=list(text='Increased tourism', side=3,  line=1, cex =  0.6))
  plot(noBio_flickr-baseMap_flickr, breaks = cuts2, col = pal2(length(cuts2)), 
       axes = F,  box = F, main = 'E. Change if bird richness is excluded',
       legend.args=list(text='Decreased tourism', side=1,  line=1, cex =  0.6),
       add = T)
  plot(sFile2, add = T, border = parkColor, lty = 3)
  
  plot(noInfra_flickr-baseMap_flickr, breaks = cuts2, col = pal2(length(cuts2)), 
       axes = F,  box = F, main = 'F. Change if infrastucture is excluded',
       legend.args=list(text='Increased tourism', side=3,  line=1, cex =  0.6))
  plot(noInfra_flickr-baseMap_flickr, breaks = cuts2, col = pal2(length(cuts2)), 
       axes = F,  box = F, main = 'F. Change if infrastucture is excluded',
       legend.args=list(text='Decreased tourism', side=1,  line=1, cex =  0.6),
       add = T)
  plot(sFile2, add = T, border = parkColor, lty = 3)
  
  
  dev.off()
}

#Execute function 
makeMapFigure()

#Pull perm importance

pullPC <- function(points, metric){

  fName <- paste(points, '_', metric, '.csv', sep = '')
  outName <- file.path(maxentOut, fName)
  
  
  outDF  <- data.frame(row.names = c('arm_national', 'arm_PAs',
                                                    'birds_national', 'birds_PAs',
                                                    'birdsTE_national', 'birdsTE_PAs'))
  
  for(i in rownames(outDF)){
    values <- read.csv(file.path(maxentOut, i, points, 'maxentResults.csv'))
    AUC <- values[6, 'Test.AUC']
    values <- values[6,grep(metric, colnames(values), value = T)]
    
    a <- 1 
    for(j in values){
      k <- colnames(values)[a]
      outDF[i, k] <- j 
      a <- a + 1
    }
    outDF[i,'AUC'] <- AUC
    
    
  }
  write.csv(outDF, outName)
}

pullPC('eBird', 'permutation.importance')
pullPC('flickr', 'permutation.importance')
pullPC('eBird', 'contribution')
pullPC('flickr', 'contribution')

