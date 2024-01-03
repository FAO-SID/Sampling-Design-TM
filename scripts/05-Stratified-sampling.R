#
# Digital Soil Mapping
# Soil Sampling Design
# Stratified Sampling on 'soil-landuse' strata
#
# Contact: Luis.RodriguezLado@fao.org
#          Marcos.Angelini@fao.org
#________________________________________________________________

  # Empty environment and cache 
    rm(list = ls())
    gc()

# Content of this script =======================================================
# The goal of this script is to perform random and regular stratified
# soil sampling designs using soil/land-use classes as 'strata'
# The distribution of samples on the strata is based on a
# weighted area distribution with the requirement of at least 2 sampling points
# by strata. Small polygons (< 100 has) are eliminated from calculations.
# Only 'strata' classes prone to be samples are included in the analyses.
# The final sample data includes 'target' points - i.e. primary sampling points,
# and 'replacement' points - i.e. points to be used as replacements in case
# that the corresponding 'target' point cannot be sampled for any reason.
# 
# 0 - Set working directory and load packages
# 1 - User-defined variables 
# 2 - Import data 
# 3 - Delineate soil strata
# 4 - Delineate land-use strata 
# 5 - Create sampling strata
# 6 - Accommodate strata to requirements
# 7 - Plot strata
# 8 - Stratified random sampling
# 9 - Plot points over strata
# 10 - Calculate replacement areas for points
# 11 - Stratified regular sampling
#________________________________________________________________

  start_time <- Sys.time()

## 0 - Set working directory and load packages =================================

  # Load packages as a vector objects
    packages <- c("sf", "terra", "tidyverse", "rmapshaper",
                  "units","plyr", "mapview", "leaflet")
    lapply(packages, require, character.only = TRUE) # Load packages
    rm(packages) # Remove object to save memory space 
  
  # Set working directory to source file location
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


## 1 - User-defined variables ==================================================
  # Path to rasters
    raster.path <- "vnm/rasters"
  # Path to shapes
    shp.path <- "vnm/shapes"
  # Path to results
    results.path <- "vnm/results/"
  # Define the number of samples to allocate
    n <- 242

    
## 2 - Import data ====================================================
  # List existing shapefiles
    list.files(shp.path, pattern = "shp$",recursive = TRUE, full.names = TRUE)
  
  # Load soil map data
    soil <- st_read(paste0(shp.path,"/soils.shp"))
  # Load land-use data
    lc <- st_read(paste0(shp.path,"/landuses.shp"))
  
  # Define a bounding area (optional).
  # If BOUNDING = TRUE, then a bounding–area shapefile (bb) must be specified
  # and the line in the code below must be uncommented
    BOUNDING = FALSE
    #bb <- st_read(paste0(shp.path,"/your_bounding_area.shp"))
  
  # Compute replacement areas (optional).
    # If REPLACEMENT = TRUE, then replacement areas are calculated around a 
    # specified buffer distance from target points within the same stratum class
    REPLACEMENT = FALSE
    distance.buffer = 500 # Distance must be adjusted to the coordinate system
    #bb <- st_read(paste0(shp.path,"/your_bounding_area.shp"))
    

# 3 - Delineate soil strata =====================================================
  
    # Clip soil data  by bounding area if REPLACEMENT = TRUE
    if (BOUNDING) {
      soil <- st_intersection(soil, bb)
    } else {
      soil <- soil
    }
  
  # Identify existing soil groups
    unique(soil$USDA_CLASS)
  # Make corrections on data
    soil$USDA_CLASS[soil$USDA_CLASS=="molisol"] <- "mollisol"
    soil$USDA_CLASS[soil$USDA_CLASS=="mollisols"] <- "mollisol"
    soil$USDA_CLASS[soil$USDA_CLASS=="ultisols"] <- "ultisol"
    soil$USDA_CLASS[soil$USDA_CLASS=="histisol"] <- "histosol"
    soil$USDA_CLASS[soil$SOIL_PH=="River Wash"] <- "riverwash"
    soil$USDA_CLASS[is.na(soil$USDA_CLASS)] <- "undetermined"
    unique(soil$USDA_CLASS)
    
  # Select classes to use (NA and water are not used)
    selected.classes <- c("inceptisol","vertisol","entisol","oxisol",
                          "mollisol","alfisol", "ultisol","histosol",
                          "undetermined")
    soil <- soil[which(soil$USDA_CLASS %in% selected.classes),]
    
  # Check polygon geometry
    soil <- st_make_valid(soil)
  # Explode polygons
    soil <- ms_explode(soil )
    # Cast geometry to polygon
    soil <- st_cast(soil, 'POLYGON')
  # Write soil strata to disk
    st_write(soil, paste0(results.path,"soil_classes.shp"), delete_dsn = TRUE)
    
  # Plot aggregated soil classes
    map = leaflet(options = leafletOptions(minZoom = 11.4)) %>%
      addTiles()
    mv <- mapview(soil["USDA_CLASS"], alpha=0, homebutton=T,
                  layer.name = "Soils", map=map)
    mv@map
    
    
# 4 - Delineate land-use strata =================================================
    
  # Clip land-use data  by bounding area if REPLACEMENT = TRUE
    if (BOUNDING) {
      lc <- st_intersection(lc, bb)
    } else {
      lc <- lc
    }
  # Identify existing land-use groups 
    unique(lc$DESCRIPTIO)
    # Load look-up table to reclassify land-use classes
    lu <- read_csv("../soil_sampling/JAM/lu_classes.csv",show_col_types = FALSE)
    # Join shapefile and look-up table
    lc <- left_join(lc,lu)
    unique(lc$LU)
    
  # Combine classes and delete classes # DELETE MANGROVES
    selected.classes <- c("Agriculture","Savanna","Forest","Grassland")
    lc <- lc[which(lc$LU %in% selected.classes),]
    unique(lc$LU)
    
  # Check and simplify geometry
    lc <- st_make_valid(lc)
  # Cast to polygons
    lc <- st_cast(lc, 'POLYGON')
  # Omit empty polygons
    lc <- na.omit(lc)
    
  # Write land use strata as shapefile
    st_write(lc, paste0(results.path,"lc_classes.shp"), delete_dsn = TRUE)
  
  # Plot aggregated land-use classes   
    map = leaflet(options = leafletOptions(minZoom = 11.4)) %>%
      addTiles()
    mv <- mapview(lc["LU"], alpha=0, homebutton=T,
                  layer.name = "Landuse", map=map)
    mv@map
    

# 5 - Create sampling strata ===================================================
    
  # Combine soil and land use layers
    soil_lc <- st_intersection(soil, lc)  
    soil_lc$soil_lc <- paste0(soil_lc$USDA_CLASS, "_", soil_lc$LU)
    soil_lc <- soil_lc %>% dplyr::select(soil_lc, geometry)
    

# 6 - Accommodate strata to requirements =======================================

  # Select by Area. Convert to area to ha and select polygons > 100 has
    soil_lc$area <- st_area(soil_lc)/10000 
    soil_lc$area <- as.vector(soil_lc$area)
    soil_lc <- soil_lc %>% 
      group_by(soil_lc) %>% 
      mutate(area = sum(area))
    soil_lc <- soil_lc[soil_lc$area > 100,]
    plot(soil_lc[1])
    
  # Replace blank spaces with underscore symbol to keep names uniform
    soil_lc$soil_lc <- str_replace_all(soil_lc$soil_lc, " ", "_")
    
  # Create a column of strata numeric codes
    soil_lc$code <- as.character(as.numeric(as.factor(soil_lc$soil_lc)))
    
  # List final strata
    unique(soil_lc$soil_lc)

  # Write final sampling strata map
    st_write(soil_lc, paste0(results.path,"strata.shp"), delete_dsn = TRUE)

# 7 - Plot strata ==============================================================
  
    # Plot strata map
    map = leaflet(options = leafletOptions(minZoom = 11.4)) %>%
      addTiles()
    mv <- mapview(soil_lc["soil_lc"], alpha=0, homebutton=T,
                  layer.name = "Strata", map=map)
    mv@map
    

# 8 - Stratified random sampling ===============================================

  # Read already created strata shapefile
    polygons <- st_read(paste0(results.path,"strata.shp"))
    if(REPLACEMENT){
      polygons = st_intersection(polygons,distance.buffer)
    }
    
  # Calculate the area of each polygon  
    polygons$area <- st_area(polygons) 
    
  # Create a new column to group polygons by a common attribute
    polygons$group <- polygons$soil_lc
  # Drop units to allow computations
    polygons <- drop_units(polygons)
    
  # Calculate the total area of all polygons in each group
    group_areas <- polygons %>%
      dplyr::group_by(group)  %>% 
      dplyr::summarize(total_area = sum(area))
  # Add a code to each group
    group_codes <- polygons %>% group_by(group) %>%
      dplyr::summarize(code = first(code)) 
  # Join polygon strata and codes   
    group_areas <- left_join(group_areas,st_drop_geometry(group_codes), by = "group")
    
  # Ensure minimum of 2 samples at each polygon in each group
    group_areas$sample_count <- 2
    
  # Calculate the number of samples per group based on relative area
    group_areas$sample_count <- 
      group_areas$sample_count + 
      round(group_areas$total_area/sum(group_areas$total_area) *
              (n-sum(group_areas$sample_count)))
  # Adjust sample size on classes
    while (sum(group_areas$sample_count) != n) {
      if (sum(group_areas$sample_count) > n) {
      # Reduce sample count for the largest polygon until total count is n
        max_index <- which.max(group_areas$sample_count)
        group_areas$sample_count[max_index] <- group_areas$sample_count[max_index] - 1
      } else {
      # Increase sample count for the smallest polygon until total count is n
        min_index <- which.min(group_areas$sample_count)
        group_areas$sample_count[min_index] <- group_areas$sample_count[min_index] + 1
      }
    }
  # Count the total samples. It must be equal to the sampling size
    sum(group_areas$sample_count) 
    
    polygons <- left_join(polygons, st_drop_geometry(group_areas),
                          by = c("soil_lc"="group"))
    polygons <- dplyr::select(polygons, soil_lc, code.x, sample_count, geometry)
    
  # Generate random points within each strata of size 3 times
  #the required samples for each strata
    x <- spatSample(x = vect(group_areas),
                    size = group_areas$sample_count * 3, method = "random")
    
  # Compute sampling points for strata
    z <- x %>% 
      st_as_sf() %>% 
      dplyr::group_by(code) %>% 
      dplyr::mutate(sample_count = as.numeric(sample_count),
                    order = seq_along(code),
                    ID = paste0(code, ".", order),
                    type = ifelse(sample_count >= order, "Target", "Replacement")) %>% 
      vect()
    
  # Find classes with missing samples
    missing.data <- left_join(group_areas,data.frame(z) %>%
                                dplyr::filter(type=="Target") %>%
                                dplyr::group_by(code) %>%
                                tally()) %>%
      dplyr::mutate(diff=sample_count-n)
    
  # Determine missing sampled strata
    missing.strata <- which(is.na(missing.data$diff))
    
  # Determine missing sampling point in strata (undersampled strata)
    missing.sample = which(missing.data$diff != 0)
    missing.number <- as.numeric(unlist(st_drop_geometry(missing.data[(missing.sample <- which(missing.data$diff != 0)),7])))
    
  # Compute sampling points for missing sampled strata
    x.missing.strata <- x[1]
    x.missing.strata$sample_count<- 0
    
    for(i in missing.strata){
      xx.missing.strata <- x[1]
      xx.missing.strata$sample_count<- 0
      nn=0
      while (sum(xx.missing.strata$sample_count) < 
             group_areas[group_areas$code==i,][["sample_count"]]*5) {
        
        while(nn < group_areas[group_areas$code==i,][["sample_count"]]*3){
          my.missing.strata <- spatSample(x = vect(group_areas[group_areas$code %in% i,]),
                                          size =  group_areas[group_areas$code==i,][["sample_count"]]*5,
                                          method = "random")
          nn <- nn + nrow(data.frame(my.missing.strata))
        }
        xx.missing.strata <- rbind(xx.missing.strata,my.missing.strata)
        print(sum(xx.missing.strata$sample_count))
      }
      print(i)
      print(xx.missing.strata)
      x.missing.strata <- rbind(x.missing.strata,xx.missing.strata)
    }
    
  # Join initial sampling with missing sampling strata data
    x <- rbind(x, x.missing.strata)
    
  # Compute sampling points for missing samples (random sampling)
    x.missing.sample <- x[1]
    
    for(i in missing.sample){
      xx.missing.sample <- x[1]
      xx.missing.sample$sample_count<- 0
      while (sum(xx.missing.sample$sample_count) < (group_areas[group_areas$code==i,][["sample_count"]]*3)) {
        my.missing.sample <- spatSample(x = vect(group_areas[group_areas$code %in% i,]),
                                        size = as.numeric(vect(group_areas[group_areas$code %in% i,])[[4]])+
                                          (group_areas[group_areas$code==i,][["sample_count"]]*3), method = "random")
        
        xx.missing.sample <- rbind(xx.missing.sample,my.missing.sample)
        print(sum(xx.missing.sample$sample_count))
      }
      print(i)
      print(xx.missing.sample)
      x.missing.sample <- rbind(x.missing.sample,xx.missing.sample)
    }
    
  # Join initial sampling with missing sampling strata data and with missing samples 
    x <- rbind(x, x.missing.sample)
    
  # Remove extra artificial replacements 
    x <- x[x$sample_count > 0,]
    
  # Convert and export to shapefile
    z <- x %>% 
      st_as_sf() %>% 
      dplyr::group_by(code) %>% 
      dplyr::mutate(sample_count = as.numeric(sample_count),
                    order = seq_along(code),
                    ID = paste0(code, ".", order),
                    type = ifelse(sample_count >= order, "Target", "Replacement")) %>% 
      vect()
    writeVector(z,
                paste0(results.path,"strat_randm_samples.shp", overwrite=TRUE))

  # Check if the number of initial target points equals the final target points 
    n==nrow(z[z$type=="Target",])
    n;nrow(z[z$type=="Target",])
    
# 9 - Plot points over strata ==================================================

    map = leaflet(options = leafletOptions(minZoom = 11.4)) %>%
      addTiles()
    mv <- mapview(soil_lc["soil_lc"], alpha=0, homebutton=T, layer.name = "Strata") + 
      mapview(sf::st_as_sf(z), zcol = 'type', color = "white", col.regions = c('royalblue', 'tomato'), cex=3, legend = TRUE,layer.name = "Samples")
    mv@map
    
    
# 10 - Calculate replacement areas for points ==================================    
    
    if(REPLACEMENT){
      # Load strata
      soil_lc <- st_read(paste0(results.path,"strata.shp"))
      
      # Read sampling. points from previous step
      z <- st_read(paste0(results.path,"strat_randm_samples.shp"))
      
      # Apply buffer of 500 meters
      buf.samples <- st_buffer(z, dist=distance.buffer)
      
      # Intersect buffers
      samples_buffer = st_intersection(soil_lc, buf.samples)
      samples_buffer <- samples_buffer[samples_buffer$type=="Target",]
      samples_buffer <- samples_buffer[samples_buffer$soil_lc==samples_buffer$group,]
      
      # Save Sampling areas
      st_write(samples_buffer, paste0(results.path,"replacement_areas.shp"),
                                      delete_dsn = TRUE)
      # Write target points only
      targets <- z[z$type=="Target",]
      st_write(targets, paste0(results.path,"target_points.shp"), delete_dsn = TRUE)
    }
    
    
# 11 - Stratified regular sampling ==================================    
    # Read already created strata shapefile
    polygons <- st_read(paste0(results.path,"strata.shp"))
    if(REPLACEMENT){
      polygons = st_intersection(polygons,distance.buffer)
    }
    
    # Calculate the area of each polygon  
    polygons$area <- st_area(polygons) 
    
    # Create a new column to group polygons by a common attribute
    polygons$group <- polygons$soil_lc
    # Drop units to allow computations
    polygons <- drop_units(polygons)
    
    # Calculate the total area of all polygons in each group
    group_areas <- polygons %>%
      dplyr::group_by(group)  %>% 
      dplyr::summarize(total_area = sum(area))
    # Add a code to each group
    group_codes <- polygons %>% group_by(group) %>%
      dplyr::summarize(code = first(code)) 
    # Join polygon strata and codes   
    group_areas <- left_join(group_areas,st_drop_geometry(group_codes), by = "group")
    
    # Ensure minimum of 2 samples at each polygon in each group
    group_areas$sample_count <- 2
    
    # Calculate the number of samples per group based on relative area
    group_areas$sample_count <- 
      group_areas$sample_count + 
      round(group_areas$total_area/sum(group_areas$total_area) *
              (n-sum(group_areas$sample_count)))
    # Adjust sample size on classes
    while (sum(group_areas$sample_count) != n) {
      if (sum(group_areas$sample_count) > n) {
        # Reduce sample count for the largest polygon until total count is n
        max_index <- which.max(group_areas$sample_count)
        group_areas$sample_count[max_index] <- group_areas$sample_count[max_index] - 1
      } else {
        # Increase sample count for the smallest polygon until total count is n
        min_index <- which.min(group_areas$sample_count)
        group_areas$sample_count[min_index] <- group_areas$sample_count[min_index] + 1
      }
    }
    # Count the total samples. It must be equal to the sampling size
    sum(group_areas$sample_count) 
    
    polygons <- left_join(polygons, st_drop_geometry(group_areas),
                          by = c("soil_lc"="group"))
    polygons <- dplyr::select(polygons, soil_lc, code.x, sample_count, geometry)
    
    # Generate regular points within each strata of size 3 times
    #the required samples for each strata
    x <- spatSample(x = vect(group_areas),
                    size = group_areas$sample_count * 3, method = "regular")
    
    # Compute sampling points for strata
    z <- x %>% 
      st_as_sf() %>% 
      dplyr::group_by(code) %>% 
      dplyr::mutate(sample_count = as.numeric(sample_count),
                    order = seq_along(code),
                    ID = paste0(code, ".", order),
                    type = ifelse(sample_count >= order, "Target", "Replacement")) %>% 
      vect()
    
    # Find classes with missing samples
    missing.data <- left_join(group_areas,data.frame(z) %>%
                                dplyr::filter(type=="Target") %>%
                                dplyr::group_by(code) %>%
                                tally()) %>%
      dplyr::mutate(diff=sample_count-n)
    
    # Determine missing sampled strata
    missing.strata <- which(is.na(missing.data$diff))
    
    # Determine missing sampling point in strata (undersampled strata)
    missing.sample = which(missing.data$diff != 0)
    missing.number <- as.numeric(unlist(st_drop_geometry(missing.data[(missing.sample <- which(missing.data$diff != 0)),7])))
    
    # Compute sampling points for missing sampled strata
    x.missing.strata <- x[1]
    x.missing.strata$sample_count<- 0
    
    for(i in missing.strata){
      xx.missing.strata <- x[1]
      xx.missing.strata$sample_count<- 0
      nn=0
      while (sum(xx.missing.strata$sample_count) < 
             group_areas[group_areas$code==i,][["sample_count"]]*5) {
        
        while(nn < group_areas[group_areas$code==i,][["sample_count"]]*3){
          my.missing.strata <- spatSample(x = vect(group_areas[group_areas$code %in% i,]),
                                          size =  group_areas[group_areas$code==i,][["sample_count"]]*5,
                                          method = "regular")
          nn <- nn + nrow(data.frame(my.missing.strata))
        }
        xx.missing.strata <- rbind(xx.missing.strata,my.missing.strata)
        print(sum(xx.missing.strata$sample_count))
      }
      print(i)
      print(xx.missing.strata)
      x.missing.strata <- rbind(x.missing.strata,xx.missing.strata)
    }
    
    # Join initial sampling with missing sampling strata data
    x <- rbind(x, x.missing.strata)
    
    # Compute sampling points for missing samples (regular sampling)
    x.missing.sample <- x[1]
    
    for(i in missing.sample){
      xx.missing.sample <- x[1]
      xx.missing.sample$sample_count<- 0
      while (sum(xx.missing.sample$sample_count) < (group_areas[group_areas$code==i,][["sample_count"]]*3)) {
        my.missing.sample <- spatSample(x = vect(group_areas[group_areas$code %in% i,]),
                                        size = as.numeric(vect(group_areas[group_areas$code %in% i,])[[4]])+
                                          (group_areas[group_areas$code==i,][["sample_count"]]*3), method = "regular")
        
        xx.missing.sample <- rbind(xx.missing.sample,my.missing.sample)
        print(sum(xx.missing.sample$sample_count))
      }
      print(i)
      print(xx.missing.sample)
      x.missing.sample <- rbind(x.missing.sample,xx.missing.sample)
    }
    
    # Join initial sampling with missing sampling strata data and with missing samples 
    x <- rbind(x, x.missing.sample)
    
    # Remove extra artificial replacements 
    x <- x[x$sample_count > 0,]
    
    # Convert and export to shapefile
    z <- x %>% 
      st_as_sf() %>% 
      dplyr::group_by(code) %>% 
      dplyr::mutate(sample_count = as.numeric(sample_count),
                    order = seq_along(code),
                    ID = paste0(code, ".", order),
                    type = ifelse(sample_count >= order, "Target", "Replacement")) %>% 
      vect()
    writeVector(z,
                paste0(results.path,"strat_randm_samples.shp", overwrite=TRUE))
    
    # Check if the number of initial target points equals the final target points 
    n==nrow(z[z$type=="Target",])
    n;nrow(z[z$type=="Target",])
    
    # Plot results
    map = leaflet(options = leafletOptions(minZoom = 11.4)) %>%
      addTiles()
    mv <- mapview(soil_lc["soil_lc"], alpha=0, homebutton=T, layer.name = "Strata") + 
      mapview(sf::st_as_sf(z), zcol = 'type', color = "white", col.regions = c('royalblue', 'tomato'), cex=3, legend = TRUE,layer.name = "Samples")
    mv@map
    
    
  ## Random Sampling based on a stratified raster
    strata <- st_read("../soil_sampling/JAM/strata.shp", quiet = TRUE)
    strata$code <- as.integer(strata$code)
    
    # Create stratification raster 
    strata <- rast(st_rasterize(strata["code"],st_as_stars(st_bbox(strata), nx = 250, ny = 250)))
    names(strata) <- "strata"
    
    # Create stratified random sampling
    srs <- sample_strat(
      sraster = strata,
      nSamp = 200
    )
    
    # Plot samples over strata
    plot(strata, main="Strata and random samples")
    points(srs,col="red")
    
    # Histogram of frequencies
    calculate_representation(
      sraster = strata,
      existing = srs,
      plot = TRUE 
    )
    
####  Time required to compute the script
Sys.time() - start_time 

