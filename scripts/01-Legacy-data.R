#
# Digital Soil Mapping
# Soil Sampling Design
# Evaluation of Legacy Data
#
# GSP-Secretariat
# Contact: Luis.RodriguezLado@fao.org

#________________________________________________________________

# Empty environment and cache 
  rm(list = ls())
  gc()

# Content of this script ========================================

# Script for evaluation the degree of representativeness of a soil legacy dataset
# relative to the diversity of the environmental conditions described in a set 
# of rasterb covariates.
# 
# 0 - Set working directory and load packages
# 1 - User-defined variables 
# 2 - Extract environmental data from rasters at soil locations
# 3 - Extract environmental data from rasters at soil locations
# 4 - Compute variability matrix in covariates
# 5 - Calculate hypercube of "covariates" distribution (P)
# 6 - Calculate hypercube of "sample" distribution (Q)
# 7 - Calculate Representativeness of the Legacy Dataset
#________________________________________________________________

start_time <- Sys.time()

## 0 - Set working directory and load packages =================================

  #remotes::install_github("lemuscanovas/synoptReg")

  # Set working directory to source file location
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  setwd("../") # Move wd down to main folder

  # List of packages
  packages <- c("sp","terra","raster","sf","clhs", "sgsR","entropy", "tripack",
                "manipulate","dplyr","synoptReg")
  # Load packages
  lapply(packages, require, character.only = TRUE)
  # Remove object to save memory space
  rm(packages) 
  
  
## 1 - User-defined variables ==================================================
  # Path to rasters
  raster.path <- "data/rasters"
  # Path to shapes
  shp.path <- "data/shapes"
  # Path to results
  results.path <- "data/results/"
  # Aggregation factor for up-scaling raster covariates (optional)
  agg.factor = 10
  
  
## 2 - Prepare data ============================================================
  ## Load raster covariate data
  # Read Spatial data covariates as rasters with terra
  cov.dat <-  list.files(raster.path, pattern = "tif$",  recursive = TRUE, full.names = TRUE)
  cov.dat <- terra::rast(cov.dat) # SpatRaster from terra
  # Aggregate stack to simplify data rasters for calculations 
  cov.dat <- aggregate(cov.dat, fact=agg.factor, fun="mean")
  
  # Load shape of district
  nghe <- sf::st_read(file.path(paste0(shp.path,"/Nghe_An.shp")),quiet=TRUE)
  
  # Crop covariates on administrative boundary
  cov.dat <- crop(cov.dat, nghe, mask=TRUE)
  
  # Transform raster information with PCA
  pca <- raster_pca(cov.dat)
  
  # Get SpatRaster layers
  cov.dat <- pca$PCA
  # Create a raster stack to be used as input in the clhs::clhs function 
  cov.dat.ras <- raster::stack(cov.dat) 
  # Subset rasters
  cov.dat <- pca$PCA[[1:first(which(pca$summaryPCA[3,]>0.99))]]
  cov.dat.ras <-  cov.dat.ras[[1:first(which(pca$summaryPCA[3,]>0.99))]]

## 3 - Extract environmental data from rasters at soil locations ===============
  # Load legacy soil data
  p.dat <- terra::vect(file.path(paste0(shp.path,"/legacy_soils.shp")))
  # Extract data
  p.dat_I <- terra::extract(cov.dat, p.dat)
  p.dat_I <- na.omit(p.dat_I) # Remove soil points outside study area
  p.dat_I.df <- p.dat_I[,-1]
  str(p.dat_I.df)

## 4 - Compute variability matrix in the covariates ====================================
  # Define Number of bins
  nb <- 25
  # Quantile matrix of the covariate data
  q.mat <- matrix(NA, nrow = (nb + 1), ncol = nlyr(cov.dat))
  j = 1
  for (i in 1:nlyr(cov.dat)) {
    ran1 <- minmax(cov.dat[[i]])[2] - minmax(cov.dat[[i]])[1]
    step1 <- ran1 / nb
    q.mat[, j] <-
      seq(minmax(cov.dat[[i]])[1],
          to = minmax(cov.dat[[i]])[2],
          by = step1)
    j <- j + 1
  }
  q.mat

## 5 - Calculate hypercube of "covariates" distribution (P)  ===================  
  # Convert SpatRaster to dataframe for calculations
  cov.dat.df <- as.data.frame(cov.dat) 
  cov.mat <- matrix(1, nrow = nb, ncol = ncol(q.mat))
  cov.dat.mx <- as.matrix(cov.dat.df)
  for (i in 1:nrow(cov.dat.mx)) {
    for (j in 1:ncol(cov.dat.mx)) {
      dd <- cov.dat.mx[[i, j]]
      
      if (!is.na(dd)) {
        for (k in 1:nb) {
          kl <- q.mat[k, j]
          ku <- q.mat[k + 1, j]
          
          if (dd >= kl && dd <= ku) {
            cov.mat[k, j] <- cov.mat[k, j] + 1
          }
        }
      }
    }
  }
  cov.mat
  
## 6 - Calculate hypercube of "sample" distribution (Q) ========================  

  h.mat <- matrix(1, nrow = nb, ncol = ncol(q.mat))
  for (i in 1:nrow(p.dat_I.df)) {
    for (j in 1:ncol(p.dat_I.df)) {
      dd <- p.dat_I.df[i, j]
      
      if (!is.na(dd)) {
        for (k in 1:nb) {
          kl <- q.mat[k, j]
          ku <- q.mat[k + 1, j]
          
          if (dd >= kl && dd <= ku) {
            h.mat[k, j] <- h.mat[k, j] + 1
          }
        }
      }
    }
  }
  h.mat 

  
## 7 - Calculate Representativeness of the Legacy Dataset ==================    

  ## Calculate the proportion of "variables" in the covariate spectra that fall within the convex hull of variables in the "environmental sample space"
  # Principal component of the legacy data sample
  pca.s = prcomp(p.dat_I[,2:(ncol(cov.dat.df)+1)],scale=TRUE, center=TRUE)
  scores_pca1 = as.data.frame(pca.s$x)
  # Plot the first 2 principal components and convex hull
  rand.tr <- tri.mesh(scores_pca1[,1],scores_pca1[,2],"remove") # Delaunay triangulation 
  rand.ch <- convex.hull(rand.tr, plot.it=F) # convex hull
  pr_poly = cbind(x=c(rand.ch$x),y=c(rand.ch$y)) # save the convex hull vertices
  plot(scores_pca1[,1], scores_pca1[,2], xlab="PCA 1", ylab="PCA 2", xlim=c(min(scores_pca1[,1:2]), max(scores_pca1[,1:2])),ylim=c(min(scores_pca1[,1:2]), max(scores_pca1[,1:2])), main='Convex hull of soil legacy data')
  lines(c(rand.ch$x,rand.ch$x[1]), c(rand.ch$y,rand.ch$y[1]),col="red",lwd=1) # draw the convex hull (domain of legacy data)
  
  # PCA projection of study area population onto the principal components
  PCA_projection <- predict(pca.s, cov.dat.df) # Project study area population onto sample PC
  newScores = cbind(x=PCA_projection[,1],y=PCA_projection[,2]) # PC scores of projected population
  
  # Plot the polygon and all points to be checked
  plot(newScores, xlab="PCA 1", ylab="PCA 2", xlim=c(min(newScores[,1:2]), max(newScores[,1:2])), ylim=c(min(newScores[,1:2]), max(newScores[,1:2])), col='black', main='Environmental space plots over the convex hull of soil legacy data')
  polygon(pr_poly,col='#99999990')
  # Check which points fall within the polygon
  pip <- point.in.polygon(newScores[,2], newScores[,1], pr_poly[,2],pr_poly[,1],mode.checked=FALSE)
  newScores <- data.frame(cbind(newScores, pip))
  # Plot points outside convex hull  
  points(newScores[which(newScores$pip==0),1:2],pch='X', col='red')
  
  #  Proportion of the conditions in the study area that fall within the convex hull of sample conditions
  sum(nrow(newScores[newScores$pip>0,]))/nrow(newScores)*100 
  
