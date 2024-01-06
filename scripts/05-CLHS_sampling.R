#
# Digital Soil Mapping
# Soil Sampling Design
# conditional Latin Hypercube Sampling
#
# GSP-Secretariat
# Contact: Luis.RodriguezLado@fao.org

#________________________________________________________________

  # Empty environment and cache 
  rm(list = ls())
  gc()

# Content of this script ========================================

# Script for creating a sampling design based on conditioned Latin Hypercube Sampling.
# Given a suite of covariates this algorithm will assess the optimal location of
#  samples based on the amount of information in the set of covariates.
# 
# 0 - Set working directory and load packages
# 1 - User-defined variables 
# 2 - Import national data 
# 3 - Compute clhs
# 4 - Including existing legacy data in a cLHS sampling design 
# 5 - Working with large raster data
# 6 - Cost–constrained cLHS sampling
# 7 - Replacement areas in cLHS design
# 8 - Polygonize replacement areas by similarity   
# 9 - Constrained cLHS sampling accounting for accessibility and legacy data  
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
  raster.path <- "data/rasters/"
  # Path to shapes
  shp.path <- "data/shapes/"
  # Path to results
  results.path <- "data/results/"
  # Path to additional data
  other.path <- "data/other/"
  # Buffer distance for replacement areas (clhs)
  D <- 1000 # Buffer distance to calculate replacement areas 
  # Define the minimum sample size. By default it uses the value calculated previously
  minimum_n <- 180
  # Aggregation factor for up-scaling raster covariates (optional)
  agg.factor = 10

## 2 - Import national data ====================================================

  # Read Spatial data covariates as rasters with terra
  cov.dat <-  list.files(raster.path, pattern = "tif$",  recursive = TRUE, full.names = TRUE)
  cov.dat <- terra::rast(cov.dat) # SpatRaster from terra
  # Load shape of district
  nghe <- sf::st_read(file.path(paste0(shp.path,"/Nghe_An.shp")),quiet=TRUE)

  # Crop covariates on administrative boundary
  cov.dat <- crop(cov.dat, nghe, mask=TRUE)
  # Store elevation and slope separately
  elevation <- cov.dat$dtm_elevation_250m
  slope <- cov.dat$dtm_slope_250m
  
  # Load roads
  roads <-  sf::st_read(file.path(paste0(shp.path,"/roads.shp")),quiet=TRUE)
  roads <-  st_intersection(roads, nghe)

  # Simplify raster information with PCA
  pca <- raster_pca(cov.dat)
  
  # Get SpatRaster layers
  cov.dat <- pca$PCA
  # Subset rasters
  cov.dat <- pca$PCA[[1:first(which(pca$summaryPCA[3,]>0.99))]]
  # Remove pca stack
  rm(pca)

  # Aggregate stack to simplify data rasters for calculations 
    # cov.dat <- aggregate(cov.dat, fact=10, fun="mean")
  
  # Plot of covariates
  plot(cov.dat)

  
## 3 - Compute clhs ============================================================
 
  # Distribute sampling points with clhs
  pts <-  sample_clhs(cov.dat, nSamp = 100, iter = 10000, progress = FALSE, simple = FALSE)
  # Plot cLHS samples on map
  plot(cov.dat[[1]], main="cLHS samples")
  points(pts, col="red", pch = 1)
  
  
## 4 - Including existing legacy data in a cLHS sampling design ================
  
  # Create an artificial random legacy dataset of 50 samples over the study area as an example
  legacy.data <- sample_srs( raster = cov.dat, nSamp = 50)
 
  # Calculate clhs 100 points plus locations of legacy data
  res <-  sample_clhs(cov.dat, nSamp = 150, existing=legacy.data, iter = 10000, progress = FALSE)

  # Plot points
  plot(cov.dat[[1]], main="cLHS samples (blue circles) and legacy samples (red diamonds)")
  points(res, col="navy", pch = 1)
  points(res[res$type=="existing",], col="red", pch = 5, cex=2)
  
  
## 5 - Working with large raster data ==========================================
  
  # Scaling covariates
    # Aggregation of covariates by a factor of 10. 
    # The original grid resolution is up-scaled using the mean value of the pixels in the grid
    cov.dat2 <- aggregate(cov.dat, fact=10, fun="mean")
    # Create clhs samples upon the resamples rasters  
    resampled.clhs <- sample_clhs(cov.dat2, nSamp = 150, iter = 10000, progress = FALSE)
    # Plot the points over the 1st raster
    plot(cov.dat2[[1]], main="Regular resampled data")
    points(resampled.clhs , col="red", pch = 1)
    rm(cov.dat2)

  # Sampling to regular points
    # Create a regular grid of 1000 points on the covariate space
    regular.sample <- spatSample(cov.dat, size = 1000, xy=TRUE, method="regular", na.rm=TRUE)
    # plot the points over the 1st raster
    plot(cov.dat[[1]], main="Regular resampled data")
    points(regular.sample, col="red", pch = 1)
    
    # Create clhs samples upon the regular grid  
    regular.sample.clhs <- clhs(regular.sample, size = 100, progress = FALSE, iter = 10000, simple = FALSE)
    # Plot points of clhs samples
    points <- regular.sample.clhs$sampled_data # Get point coordinates of clhs sampling
    plot(cov.dat[[1]], main="cLHS samples (red) and covariate resampled points (blue)")
    points(regular.sample, col="navy", pch = 1)
    points(points, col="red", cex=1)
   
  
## 6 - Cost–constrained cLHS sampling
      # Use 'Distance to roads' as cost surface
      # Calculate distance to roads with te same spatial definition than the covariates
      # dist2access <- terra::distance(cov.dat[[1]], roads, progress=TRUE)
      # names(dist2access) <- "dist2access"
      # Save cost surface to disk
      # writeRaster(dist2access, paste0(other.path,"nghe_d2roads.tif"), overwrite=TRUE)
      
    # Load pre-calculated distance–to–roads surface
    dist2access <- terra::rast(paste0(other.path,"nghe_d2roads.tif"))
    # Aggregate to the same soatial definition
     # dist2access <- aggregate(dist2access, fact=10, fun="mean")
    plot(dist2access)
    plot(nghe, col="transparent", add=TRUE)
    
    # Add cost surface (distance to roads) as layer
    cov.dat <- c(cov.dat,dist2access)
    names(cov.dat)
    
    # Harmonize NAs
    cov.dat$dist2access <- cov.dat$dist2access * cov.dat[[1]]/cov.dat[[1]]
    plot(cov.dat$dist2access)
    plot(nghe, col="transparent",add=TRUE)
    
    # Compute sampling points
    cost.clhs <- sample_clhs(cov.dat, nSamp = minimum_n, iter = 10000, progress = FALSE, cost = 'dist2access',  use.cpp = TRUE)

    # Plot samples
    plot(cov.dat[['dist2access']], main="cLHS samples with 'cost' constraints")
    plot(cost.clhs, col="red", cex=1,add=T)
    

## 7 - Replacement areas in cLHS design (requires raster package)
    
    # Determine the similarity to points in a buffer of distance 'D'
    # Compute the buffers around points # cov25??
    gw <- similarity_buffer(raster::stack(cov.dat) ,  as(cost.clhs, "Spatial"), buffer = D)

    # Plot results
    plot(elevation, legend=TRUE,main=paste("Similarity probability over elevation"))
    ## Overlay points
    points(cost.clhs, col = "dodgerblue", pch = 3)
    ## Overlay probability stack for point 1
    colors <- c((RColorBrewer::brewer.pal(9, "YlOrRd")))
    terra::plot(gw[[1]], add=TRUE ,  legend=FALSE, col=colors)
    ## Overlay 1st cLHS point
    points(cost.clhs[1,], col = "red", pch = 12,cex=1)
    
# 8 - Polygonize replacement areas by similarity    

    # Determine a threshold break to delineate replacement areas
    similarity_threshold <- 0.90
    # Reclassify buffer raster data according to the threshold break of probability
    # 1 = similarity >= similarity_break; NA =  similarity <  similarity_break
    # Define a vector with the break intervals and the output values (NA,1) 
    breaks <- c(0, similarity_threshold, NA, similarity_threshold, 1, 1)
    # Convert to a matrix
    breaks <- matrix(breaks, ncol=3, byrow=TRUE)
    # Reclassify the data in the layers from probabilities to (NA,)
    s = stack(lapply(1:raster::nlayers(gw), function(i){raster::reclassify(gw[[i]], breaks, right=FALSE)}))
    
    # Polygonize replacement areas 
    s = lapply(as.list(s), rasterToPolygons, dissolve=TRUE)
    s <- bind(s,keepnames=TRUE)
    # Add the identifier of the corresponding target point
    for(i in 1: length(s)){
      s@data$ID[i] <- as.integer(stringr::str_replace(s@polygons[[i]]@ID,"1.",""))
    }
    # Clean the data by storing target ID data only
    s@data <- s@data["ID"]
    
    # Plot results
    plot(cov.dat[[1]], main=paste("cLHS samples and replacement areas for threshold = ", similarity_threshold))
    plot(s,add=TRUE, col=NA, border="gray40")
    points(cost.clhs, col="red", pch = 3)
    
    # Export replacement areas to shapefile 
    s <- st_as_sf(s)
    st_write(s, file.path(paste0(results.path,'replacement_areas_', D, '.shp')), delete_dsn = TRUE)
    
    # Write cLHS sampling points to shapefile
    cost.clhs$ID <- row(cost.clhs)[,1] # Add identifier
    out.pts <- st_as_sf(cost.clhs)
    st_write(out.pts, paste0(results.path,'target_clhs.shp'), delete_dsn = TRUE)

         
## 9 - Constrained cLHS sampling taking into account accessibility and legacy data  
    # Load legacy data 
    legacy <- sf::st_read(file.path(paste0(shp.path,"/legacy_soils.shp")),quiet=TRUE)
    
    # Add slope to the stack
    cov.dat$slope <- slope
    # Calculate clhs points with legacy, cost and buffer to roads
    buff_inner=20;
    buff_outer=3000
    # Convert roads to sf object and cast to multilinestring
    roads <- st_as_sf(roads) %>%
      st_cast("MULTILINESTRING") 
    
    # Calculate clhs samples using slope as cost surface, distance to roads as
    # access limitations, and including existing legacy data
    aa <- sgsR::sample_clhs(mraster = cov.dat, nSamp = minimum_n, existing = legacy,
                            iter = 10000, details = TRUE, cost="slope", access=roads,
                            buff_inner=buff_inner, buff_outer=buff_outer)
    
    ## Plot distances, roads, clhs points and legacy data 
    plot(cov.dat$dist2access, main="New (red) and existing (blue) samples")
    plot(roads,add=TRUE, col="gray60")
    plot(aa$samples[aa$samples$type=="new",], col= "tomato",add=TRUE)
    plot(aa$samples[aa$samples$type=="existing",], col= "dodgerblue2",add=TRUE)
    
    # Write samples as shapefile
    aa$samples[c("type","dist2access")] %>%
      st_write(paste0(results.path,'const_clhs.shp'), delete_dsn = TRUE)
    
    
## 10 - Calculate COOBS ========================================================
    calculate_coobs(
      mraster = cov.dat,
      existing = cost.clhs,
      plot = TRUE,
      cores = 4
    )
    
    
