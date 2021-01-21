#Libraries
library(raster)
library(rgbif)
library(vegan)
library(ggplot2)
library(dismo)
library(rjava)
library(sf)
library(tidyr)


#download costa rica shapefile 
cr  <- getData("GADM",country="Costa Rica",level=0)

#download climate data 
climStack <- getData("worldclim", var='bio', res=0.5, lon=bbox(cr)[1,], lat=bbox(cr)[2,])

#crop climate data to CR
climStack <- crop(climStack, cr)
climStack <- mask(climStack, cr)


#makeSDM takes a scientific name and downloads the GBIF points you want to use 
makeSDM <- function(species){

  #download GBIF points
  gbifOccs <- occ_search(scientificName = species, geometry = bbox(cr), limit = 1000000)
  gbifOccs <- as.data.frame(gbifOccs$data)
  points <- gbifOccs[,c('decimalLongitude', 'decimalLatitude')]

  #Run maxent model and make predictions 
  myModel <- maxent(climStack, points)
  myRaster <- predict(myModel, climStack)
  
  #Set up return elements
  toReturn <- list("mxaentModel" = myModel, "maxentRaster" = myRaster)
  
  }

#Example 
modelResults <- makeSDM('Apis mellifera')
