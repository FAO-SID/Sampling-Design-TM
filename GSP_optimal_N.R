## Sample number optimisation
## Adaptation of the Script for Sample number optimisation by Malone et al.
## The script has been modified to work with raster (tif) and shape data
## The original data set in the publication has been transformed to tif and shape for as a reproducible example.
## Code chunks for raster aggregation and cropping study areas have been included but are commented in the file. 
## The aim of this script is to compare the collected data to the actual population in terms of covariate coverage
#https://www.youtube.com/watch?v=RdSmYvbQkhs&t=304s

# Load required libraries
  packages <- c("sp", "terra", "clhs", "entropy", "tripack","manipulate")
  lapply(packages, require, character.only = TRUE)
  rm(packages)

# Set working directory to source file location
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path));

## Load soil legacy point data
  p.dat <- terra::vect("data/shapes/legacy_soils.shp")

## Load raster covariate data----
  # Read Spatial data covariates as rasters with terra
  rasters <- "data/rasters"
  cov.dat <-  list.files(rasters, pattern = "tif$",  recursive = TRUE, full.names = TRUE)
  cov.dat <- terra::rast(cov.dat)
  plot(cov.dat)
  
  # Aggregate the information at pixel level for higher speed
  # cov.aggregated <- terra::aggregate(cov.dat, fact = 10, na.rm=TRUE) # Aggregate data for simplicity
  # Recalculate the aggregation values for a categorical dataset (example with landuse)
  # cov.aggregated[["landuse"]] <- resample(cov.dat[["landuse"]], cov.aggregated[["landuse"]], "near") 
  # cov.dat <- cov.aggregated 
  # plot(cov.dat)
  # rm(cov.aggregated) # remove intermediary object to save memory space

  # Crop rasters by study area
  # crop.area <- terra::vect("data/shapes/crop_area.shp") # Load shape with cropping boundaries
  # cov.dat <- crop(cov.dat, crop.area)
  # plot(cov.dat)
  
# Extract environmental data from rasters at soil locations ----
  p.dat_I <- terra::extract(cov.dat, p.dat)
  p.dat_I <- na.omit(p.dat_I) # Remove soil points outside study area
  str(p.dat_I)
  #p.dat_I <- p.dat_I[,2:(ncol(cov.dat.df)+1)]
  
## Comparison of population and sample distributions ----

  # Kullback-Leibler (KL) divergence
    
  # Quantiles of the population
    # Number of bins
      nb<- 25
    #quantile matrix (of the covariate data)
      q.mat<- matrix(NA, nrow=(nb+1), ncol= nlyr(cov.dat))
      j=1
      for (i in 1:nlyr(cov.dat)){ #note the index start here
      #get a quantile matrix together of the covariates
        ran1 <- minmax(cov.dat[[i]])[2] - minmax(cov.dat[[i]])[1]
        step1<- ran1/nb 
        q.mat[,j]<- seq(minmax(cov.dat[[i]])[1], to = minmax(cov.dat[[i]])[2], by =step1)
        j<- j+1}
      q.mat
    
    # Hypercube of population
      cov.dat.df <- as.data.frame(cov.dat) # convert SpatRaster to dataframe
      cov.mat<- matrix(1, nrow=nb, ncol=ncol(q.mat))
      for (i in 1:nrow(cov.dat.df)){ # the number of pixels
        cntj<- 1 
        for (j in 1:ncol(cov.dat.df)){ #for each column
          dd<- cov.dat.df[i,j]  
          for (k in 1:nb){  #for each quantile
            kl<- q.mat[k, cntj] 
            ku<- q.mat[k+1, cntj] 
            if (is.na(dd)) {
              print('Missing')
            }
            else if (dd >= kl & dd <= ku){cov.mat[k, cntj]<- cov.mat[k, cntj] + 1} 
          }
          cntj<- cntj+1
        }
      }
      cov.mat
    
    
    # Compare whole study area covariate space with the selected sample
      # Sample data hypercube (essentially the same script as for the grid data but just doing it on the sample data)
      h.mat<- matrix(1, nrow=nb, ncol=ncol(q.mat))
    
      for (ii in 1:nrow(p.dat_I)){ # the number of observations
        cntj<- 1 
        for (jj in 2:ncol(p.dat_I)){ #for each column
          dd<- p.dat_I[ii,jj]  
          for (kk in 1:nb){  #for each quantile
            kl<- q.mat[kk, cntj] 
            ku<- q.mat[kk+1, cntj] 
            if (dd >= kl & dd <= ku){h.mat[kk, cntj]<- h.mat[kk, cntj] + 1}
          }
          cntj<- cntj+1
        }
      }
      h.mat 
    
      
   ## Compute Kullback-Leibler (KL) divergence
      kl.index <-c()
      for(i in 1:ncol(cov.dat.df)){
      kl <-    KL.empirical(c(cov.mat[,i]), c(h.mat[,i]))
      kl.index <- c(kl.index,kl)
      klo <-  mean(kl.index)
      }
      kl.index # KL divergences of each covariate
      klo # KL divergence in the existing soil sample ( N= number of soil legacy point data) )


#### Representativeness of the Legacy Dataset: ----
  ## Calculate the proportion of "env. variables" in the covariate spectra that fall within the convex hull of variables in the "environmental sample space"
      
  # Principal component of the legacy data sample
    pca.s = prcomp(p.dat_I[,2:(ncol(cov.dat.df)+1)],scale=TRUE, center=TRUE)
    scores_pca1 = as.data.frame(pca.s$x)
  # Plot the first 2 principal components and convex hull
    rand.tr <- tri.mesh(scores_pca1[,1],scores_pca1[,2],"remove") # Delaunay triangulation 
    rand.ch<-convex.hull(rand.tr, plot.it=F) # convex hull
    pr_poly = cbind(x=c(rand.ch$x),y=c(rand.ch$y)) # save the convex hull vertices
    plot(scores_pca1[,1], scores_pca1[,2], xlab="PCA 1", ylab="PCA 2", xlim=c(min(scores_pca1[,1:2]), max(scores_pca1[,1:2])),ylim=c(min(scores_pca1[,1:2]), max(scores_pca1[,1:2])), main='Convex hull of soil legacy data')
    lines(c(rand.ch$x,rand.ch$x[1]), c(rand.ch$y,rand.ch$y[1]),col="red",lwd=1) # draw the convex hull (domain of legacy data)
    
  # PCA projection of study area population onto the principal components
    PCA_projection<- predict(pca.s, cov.dat.df) # Project study area population onto sample PC
    newScores = cbind(x=PCA_projection[,1],y=PCA_projection[,2]) # PC scores of projected population
    
  # Plot the polygon and all points to be checked
    plot(newScores, xlab="PCA 1", ylab="PCA 2", xlim=c(min(newScores[,1:2]), max(newScores[,1:2])), ylim=c(min(newScores[,1:2]), max(newScores[,1:2])), col='black', main='Environmental space plots over the convex hull of soil legacy data')
    polygon(pr_poly,col='#99999990')
  # Check which points fall within the polygon
    pip <- point.in.polygon(newScores[,2], newScores[,1], pr_poly[,2],pr_poly[,1],mode.checked=FALSE)
    newScores <- data.frame(cbind(newScores, pip))
  # Plot points outside convex hull  
    points(newScores[which(newScores$pip==0),1:2],pch='X', col='red')
  # Proportion of the conditions in the study area that fall within the convex hull
    sum(newScores$pip)/nrow(newScores)*100 

    

    plot(cov.dat[[1]])
    points(p.dat, col="red")
    points(p.dat_I_sp,col="blue")

    
    
