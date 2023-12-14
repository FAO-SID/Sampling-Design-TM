packages <- c("sp", "terra", "raster","clhs", "entropy", "tripack","manipulate","dplyr","plotly")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Path to rasters
raster.path <- "clhs_data/rasters"
# Path to shapes
shp.path <- "clhs_data/shapes"
# Aggregation factor
agg.factor = 2
# Buffer distance for replacement areas in clhs
D <- 200
optimal_n <- 500

# Read Spatial data covariates as rasters with terra
cov.dat <-  list.files(raster.path, pattern = "tif$",  recursive = TRUE, full.names = TRUE)
cov.dat <- terra::rast(cov.dat) # SpatRaster from terra

## Load raster covariate data
# Read Spatial data covariates as rasters with terra
cov.dat <-  list.files(raster.path, pattern = "tif$",  recursive = TRUE, full.names = TRUE)
cov.dat <- terra::rast(cov.dat)

# Aggregate to simplify data rasters for calculations 
cov.dat <- aggregate(cov.dat, fact=agg.factor, fun="mean") # Now you can see the statistics
# Create a raster stack
cov.dat.ras <- raster::stack(cov.dat) 

# Plot of covariates
plot(cov.dat)

# Distribute sampling points with clhs
phychem.pts <- clhs(cov.dat.ras, size = optimal_n * 0.9, iter = 5000, progress = FALSE, simple = FALSE)
# Plot of objective function
plot(phychem.pts, c('obj'))

## Create a cLHS sampling point set----
plot(cov.dat[[1]], main="cLHS samples")
points(phychem.pts$sampled_data, col="red", pch = 1)

# Extract environmental data from rasters at soil locations
phychem.data <- phychem.pts$sampled_data@data


#### Diversity of covariates- -----


## Variability matrix in the covariates
nb<- 25 # Define Number of bins
q.mat<- matrix(NA, nrow=(nb+1), ncol= nlyr(cov.dat)) # quantile matrix of covariate data
j=1
for (i in 1:nlyr(cov.dat)){
  ran1 <- minmax(cov.dat[[i]])[2] - minmax(cov.dat[[i]])[1]
  step1<- ran1/nb 
  q.mat[,j]<- seq(minmax(cov.dat[[i]])[1], to = minmax(cov.dat[[i]])[2], by =step1)
  j<- j+1
}

q.mat


## Hypercube of "objective" distribution (P) - covariates
cov.dat.df <- as.data.frame(cov.dat) # Convert SpatRaster to dataframe for calculations
cov.mat<- matrix(1, nrow=nb, ncol=ncol(q.mat))
for (i in 1:nrow(cov.dat.df)){ # the number of pixels
  cntj<- 1 
  for (j in 1:ncol(cov.dat.df)){ #for each column
    dd<- cov.dat.df[i,j]  
    for (k in 1:nb){  #for each quantile
      kl<- q.mat[k, cntj] 
      ku<- q.mat[k+1, cntj] 
      if (dd >= kl & dd <= ku){cov.mat[k, cntj]<- cov.mat[k, cntj] + 1} 
    }
    cntj<- cntj+1
  }
}
cov.mat

## Sample data hypercube
h.mat<- matrix(1, nrow=nb, ncol=ncol(q.mat))

for (ii in 1:nrow(phychem.data)){ # the number of observations
  cntj<- 1 
  for (jj in 1:ncol(phychem.data)){ #for each column
    dd <- phychem.data[ii,jj]  
    for (kk in 1:nb){  #for each bin
      kl <- q.mat[kk, cntj] 
      ku <- q.mat[kk+1, cntj] 
      if (dd >= kl & dd <= ku){h.mat[kk, cntj] <- h.mat[kk, cntj] + 1}
    }
    cntj<- cntj+1
  }
}

## Compare covariate distributions in P and Q with Kullback-Leibler (KL) divergence
kl.index <-c()
for(i in 1:ncol(cov.dat.df)){
  kl <-    KL.empirical(c(cov.mat[,i]), c(h.mat[,i]))
  kl.index <- c(kl.index,kl)
  klo <-  mean(kl.index)
}
print(kl.index) # KL divergences of each covariate
print(klo) # KL divergence in the existing soil samples


## Calculate the proportion of "env. variables" in the covariate spectra that fall within the convex hull of variables in the "environmental sample space"

# Principal component of the legacy data sample
pca.s = prcomp(phychem.data,scale=TRUE, center=TRUE)
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
polygon(pr_poly,col=adjustcolor("dodgerblue1", alpha.f=0.5), border="gray" )

# Check which points fall within the polygon
pip <- point.in.polygon(newScores[,2], newScores[,1], pr_poly[,2],pr_poly[,1],mode.checked=FALSE)
newScores <- data.frame(cbind(newScores, pip))
# Plot points outside convex hull  
points(newScores[which(newScores$pip==0),1:2],pch=3, col='red')
# Proportion of the conditions in the study area that fall within the convex hull
explained.1 <- sum(nrow(newScores[newScores$pip>0,]))/nrow(newScores)*100 
explained.1

# Second cLHS- Samples on samples ------

# Distribute sampling points with clhs
spect.pts <- clhs(phychem.pts$sampled_data, size = optimal_n *0.1, iter = 5000, progress = FALSE, simple = FALSE)
# Plot of objective function
plot(spect.pts, c('obj'))

spect.data <- spect.pts$sampled_data@data

#### Diversity of first sampling- -----


## Variability matrix in the  first sampling
nb<- 25 # Define Number of bins
q.mat<- matrix(NA, nrow=(nb+1), ncol= ncol(phychem.data)) # quantile matrix of covariate data
j=1
for (i in 1:ncol(phychem.data)){
  ran1 <- max(phychem.data[,i]) - min(phychem.data[,i])
  step1<- ran1/nb 
  q.mat[,j]<- seq(min(phychem.data[,i]), to = max(phychem.data[,i]), by =step1)
  j<- j+1
}

q.mat


## Hypercube of "objective" distribution (P) -  first sampling
cov.dat.df <- phychem.data # Convert SpatRaster to dataframe for calculations
cov.mat<- matrix(1, nrow=nb, ncol=ncol(q.mat))
for (i in 1:nrow(cov.dat.df)){ # the number of pixels
  cntj<- 1 
  for (j in 1:ncol(cov.dat.df)){ #for each column
    dd<- cov.dat.df[i,j]  
    for (k in 1:nb){  #for each quantile
      kl<- q.mat[k, cntj] 
      ku<- q.mat[k+1, cntj] 
      if (dd >= kl & dd <= ku){cov.mat[k, cntj]<- cov.mat[k, cntj] + 1} 
    }
    cntj<- cntj+1
  }
}
cov.mat

## Sample data hypercube
h.mat<- matrix(1, nrow=nb, ncol=ncol(q.mat))

for (ii in 1:nrow(spect.data)){ # the number of observations
  cntj<- 1 
  for (jj in 1:ncol(spect.data)){ #for each column
    dd <- spect.data[ii,jj]  
    for (kk in 1:nb){  #for each bin
      kl <- q.mat[kk, cntj] 
      ku <- q.mat[kk+1, cntj] 
      if (dd >= kl & dd <= ku){h.mat[kk, cntj] <- h.mat[kk, cntj] + 1}
    }
    cntj<- cntj+1
  }
}

## Compare covariate distributions in P and Q with Kullback-Leibler (KL) divergence
kl.index <-c()
for(i in 1:ncol(cov.dat.df)){
  kl <-    KL.empirical(c(cov.mat[,i]), c(h.mat[,i]))
  kl.index <- c(kl.index,kl)
  klo <-  mean(kl.index)
}
print(kl.index) # KL divergences of each covariate
print(klo) # KL divergence in the existing soil samples


## Calculate the proportion of "env. variables" in the covariate spectra that fall within the convex hull of variables in the "environmental sample space"

# Principal component of the second sampling data
pca.s = prcomp(spect.data,scale=TRUE, center=TRUE)
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
polygon(pr_poly,col=adjustcolor("dodgerblue1", alpha.f=0.5), border="gray" )

# Check which points fall within the polygon
pip <- point.in.polygon(newScores[,2], newScores[,1], pr_poly[,2],pr_poly[,1],mode.checked=FALSE)
newScores <- data.frame(cbind(newScores, pip))
# Plot points outside convex hull  
points(newScores[which(newScores$pip==0),1:2],pch=3, col='red')
# Proportion of the conditions in the study area that fall within the convex hull
explained.2 <- sum(nrow(newScores[newScores$pip>0,]))/nrow(newScores)*100 
explained.2

# % representativenes
explained.1
explained.2
explained.2*100/explained.1

plot(cov.dat[[1]])
points(phychem.pts$sampled_data)
points(spect.pts$sampled_data,col="red")


#### Diversity of second samples on covariates- -----


## Variability matrix in the covariates
nb<- 25 # Define Number of bins
q.mat<- matrix(NA, nrow=(nb+1), ncol= nlyr(cov.dat)) # quantile matrix of covariate data
j=1
for (i in 1:nlyr(cov.dat)){
  ran1 <- minmax(cov.dat[[i]])[2] - minmax(cov.dat[[i]])[1]
  step1<- ran1/nb 
  q.mat[,j]<- seq(minmax(cov.dat[[i]])[1], to = minmax(cov.dat[[i]])[2], by =step1)
  j<- j+1
}

q.mat


## Hypercube of "objective" distribution (P) - covariates
cov.dat.df <- as.data.frame(cov.dat) # Convert SpatRaster to dataframe for calculations
cov.mat<- matrix(1, nrow=nb, ncol=ncol(q.mat))
for (i in 1:nrow(cov.dat.df)){ # the number of pixels
  cntj<- 1 
  for (j in 1:ncol(cov.dat.df)){ #for each column
    dd<- cov.dat.df[i,j]  
    for (k in 1:nb){  #for each quantile
      kl<- q.mat[k, cntj] 
      ku<- q.mat[k+1, cntj] 
      if (dd >= kl & dd <= ku){cov.mat[k, cntj]<- cov.mat[k, cntj] + 1} 
    }
    cntj<- cntj+1
  }
}
cov.mat

## Sample data hypercube
h.mat<- matrix(1, nrow=nb, ncol=ncol(q.mat))

for (ii in 1:nrow(spect.data)){ # the number of observations
  cntj<- 1 
  for (jj in 1:ncol(spect.data)){ #for each column
    dd <- spect.data[ii,jj]  
    for (kk in 1:nb){  #for each bin
      kl <- q.mat[kk, cntj] 
      ku <- q.mat[kk+1, cntj] 
      if (dd >= kl & dd <= ku){h.mat[kk, cntj] <- h.mat[kk, cntj] + 1}
    }
    cntj<- cntj+1
  }
}

## Compare covariate distributions in P and Q with Kullback-Leibler (KL) divergence
kl.index <-c()
for(i in 1:ncol(cov.dat.df)){
  kl <-    KL.empirical(c(cov.mat[,i]), c(h.mat[,i]))
  kl.index <- c(kl.index,kl)
  klo <-  mean(kl.index)
}
print(kl.index) # KL divergences of each covariate
print(klo) # KL divergence in the existing soil samples


## Calculate the proportion of "env. variables" in the covariate spectra that fall within the convex hull of variables in the "environmental sample space"

# Principal component of the legacy data sample
pca.s = prcomp(spect.data,scale=TRUE, center=TRUE)
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
polygon(pr_poly,col=adjustcolor("dodgerblue1", alpha.f=0.5), border="gray" )

# Check which points fall within the polygon
pip <- point.in.polygon(newScores[,2], newScores[,1], pr_poly[,2],pr_poly[,1],mode.checked=FALSE)
newScores <- data.frame(cbind(newScores, pip))
# Plot points outside convex hull  
points(newScores[which(newScores$pip==0),1:2],pch=3, col='red')
# Proportion of the conditions in the study area that fall within the convex hull
explained.3 <- sum(nrow(newScores[newScores$pip>0,]))/nrow(newScores)*100 


# % representativenes
explained.1 # phychem.pts over covariates
explained.2 # spect.pts over phychem.pts
explained.3 # spect.pts over covariates

explained.2*100/explained.1

plot(cov.dat[[1]])
points(phychem.pts$sampled_data)
points(spect.pts$sampled_data,col="red")


### ----


# Distribute sampling points with clhs
spec_cov.pts <- clhs(cov.dat.ras, size = optimal_n * 0.1, iter = 5000, progress = FALSE, simple = FALSE)
# Plot of objective function
plot(spec_cov.pts, c('obj'))

## Create a cLHS sampling point set----
plot(cov.dat[[1]], main="cLHS samples")
points(spec_cov.pts$sampled_data, col="red", pch = 1)

# Extract environmental data from rasters at soil locations
spec_cov.data <- spec_cov.pts$sampled_data@data


#### Diversity of covariates- -----


## Variability matrix in the covariates
nb<- 25 # Define Number of bins
q.mat<- matrix(NA, nrow=(nb+1), ncol= nlyr(cov.dat)) # quantile matrix of covariate data
j=1
for (i in 1:nlyr(cov.dat)){
  ran1 <- minmax(cov.dat[[i]])[2] - minmax(cov.dat[[i]])[1]
  step1<- ran1/nb 
  q.mat[,j]<- seq(minmax(cov.dat[[i]])[1], to = minmax(cov.dat[[i]])[2], by =step1)
  j<- j+1
}

q.mat


## Hypercube of "objective" distribution (P) - covariates
cov.dat.df <- as.data.frame(cov.dat) # Convert SpatRaster to dataframe for calculations
cov.mat<- matrix(1, nrow=nb, ncol=ncol(q.mat))
for (i in 1:nrow(cov.dat.df)){ # the number of pixels
  cntj<- 1 
  for (j in 1:ncol(cov.dat.df)){ #for each column
    dd<- cov.dat.df[i,j]  
    for (k in 1:nb){  #for each quantile
      kl<- q.mat[k, cntj] 
      ku<- q.mat[k+1, cntj] 
      if (dd >= kl & dd <= ku){cov.mat[k, cntj]<- cov.mat[k, cntj] + 1} 
    }
    cntj<- cntj+1
  }
}
cov.mat

## Sample data hypercube
h.mat<- matrix(1, nrow=nb, ncol=ncol(q.mat))

for (ii in 1:nrow(spec_cov.data)){ # the number of observations
  cntj<- 1 
  for (jj in 1:ncol(spec_cov.data)){ #for each column
    dd <- spec_cov.data[ii,jj]  
    for (kk in 1:nb){  #for each bin
      kl <- q.mat[kk, cntj] 
      ku <- q.mat[kk+1, cntj] 
      if (dd >= kl & dd <= ku){h.mat[kk, cntj] <- h.mat[kk, cntj] + 1}
    }
    cntj<- cntj+1
  }
}

## Compare covariate distributions in P and Q with Kullback-Leibler (KL) divergence
kl.index <-c()
for(i in 1:ncol(cov.dat.df)){
  kl <-    KL.empirical(c(cov.mat[,i]), c(h.mat[,i]))
  kl.index <- c(kl.index,kl)
  klo <-  mean(kl.index)
}
print(kl.index) # KL divergences of each covariate
print(klo) # KL divergence in the existing soil samples


## Calculate the proportion of "env. variables" in the covariate spectra that fall within the convex hull of variables in the "environmental sample space"

# Principal component of the legacy data sample
pca.s = prcomp(spec_cov.data,scale=TRUE, center=TRUE)
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
polygon(pr_poly,col=adjustcolor("dodgerblue1", alpha.f=0.5), border="gray" )

# Check which points fall within the polygon
pip <- point.in.polygon(newScores[,2], newScores[,1], pr_poly[,2],pr_poly[,1],mode.checked=FALSE)
newScores <- data.frame(cbind(newScores, pip))
# Plot points outside convex hull  
points(newScores[which(newScores$pip==0),1:2],pch=3, col='red')
# Proportion of the conditions in the study area that fall within the convex hull
explained.4 <- sum(nrow(newScores[newScores$pip>0,]))/nrow(newScores)*100 

# % representativenes
explained.1 # phychem.pts over covariates
explained.2 # spect.pts over phychem.pts
explained.3 # spect.pts vs covariates
explained.4 # spect.pts on covariates



