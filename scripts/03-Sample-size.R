#
# Digital Soil Mapping
# Soil Sampling Design
# Optimizing Sample Size
#
# GSP-Secretariat
# Contact: Luis.RodriguezLado@fao.org

#________________________________________________________________

#Empty environment and cache 
rm(list = ls())
gc()

# Content of this script ========================================
# The goal of this script is to determine the minimum sample size required to describe an area
# while retaining for a 95% of coincidence in the environmental variability of covariates
# in the area 
# 
# 0 - Set working directory and load necessary packages
# 1 - User-defined variables 
# 2 - Import national data 
# 3 - Calculate the minimum sample size to describe the area
# 4 - Plot covariate diversity as PCA scores 
# 5 - KL divergence and % similarity results for growing N samples
# 6 - Model KL divergence
# 7 - Determine the minimum sample size for 95% coincidence
# 8 - Determine the optimal iteration according to the minimum N size 
# 9 - Plot minimum points from best iteration
# 10 - Calculate COOBS (sgsR)
# 11 - Calculate minimum and and optimal sample size with opendsm
#________________________________________________________________


## 0 - Set working directory and load necessary packages =======================

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../") # Move wd down to main folder

# List of packages
packages <- c("sp","terra","raster","sf","clhs",
              "sgsR","entropy", "tripack",
              "manipulate","dplyr","plotly","synoptReg")

# Load packages
lapply(packages, require, character.only = TRUE) 
rm(packages) # Remove object to save memory space


## 1 - User-defined variables ==================================================
# Path to rasters
raster.path <- "data/rasters/"
# Path to shapes
shp.path <- "data/shapes/"
# Path to results
results.path <- "data/results/"
# Path to additional data
other.path <- "data/other/"
# Aggregation factor for up-scaling raster covariates (optional)
agg.factor = 10

# Define parameters to determine minimum sampling size
initial.n <- 60 # Initial sampling size to test
final.n <- 360 # Final sampling size to test
by.n <- 20 # Increment size
iters <- 5 # Number of trials on each size


## 2 - Import national data ====================================================
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
# Store elevation and slope separately
elevation <- cov.dat$dtm_elevation_250m
slope <- cov.dat$dtm_slope_250m

# Load roads
roads <-  vect(file.path(paste0(shp.path,"/roads.shp")))
roads <- crop(roads, nghe)

# Simplify raster information with PCA
pca <- raster_pca(cov.dat)

# Get SpatRaster layers
cov.dat <- pca$PCA
# Subset rasters
cov.dat <- pca$PCA[[1:first(which(pca$summaryPCA[3,]>0.99))]]
# Plot covariates
plot(cov.dat)


# 3 - Calculate the minimum sample size to describe the area ===================
# Start computations ----
# Initialize empty vectors to store results
number_of_samples <- c()
prop_explained <- c()
klo_samples <- c()
samples_storage <- list()

# Convert SpatRaster to dataframe
cov.dat.df <-  as.data.frame(cov.dat)

# Start evaluation with growing sample sizes
for (trial in seq(initial.n, final.n, by = by.n)) {
  for (iteration in 1:iters) {
    # Generate stratified clhs samples
    p.dat_I <-  sample_clhs(cov.dat,
                     nSamp = trial, iter = 10000,
                     progress = FALSE, simple = FALSE)
    
    # Get covariate values as dataframe and delete NAs, avoid geometry
    p.dat_I.df <- as.data.frame(p.dat_I) %>%
      dplyr::select(!(geometry)) %>%
      na.omit()
    
    # Store samples as list
    samples_storage[[paste0("N", trial, "_", iteration)]] <- p.dat_I
    
    ## Comparison of population and sample distributions - Kullback-Leibler (KL) divergence
    # Define quantiles of the study area (number of bins)
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
    
    # Hypercube of covariates in study area
    # Initialize the covariate matrix
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
    
    # Compare whole study area covariate space with the selected sample
    # Sample data hypercube (the same as for the raster data but on the sample data)
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
    
    ## Compute Kullback-Leibler (KL) divergence
    kl.index <- c()
    for (i in 1:ncol(cov.dat.df)) {
      kl <-    KL.empirical(c(cov.mat[, i]), c(h.mat[, i]))
      kl.index <- c(kl.index, kl)
      klo <-  mean(kl.index)
    }
    
    ## Calculate the proportion of "env. variables" in the covariate spectra that fall within the convex hull of variables in the "environmental sample space"
    # Principal component of the data sample
    pca.s = prcomp(p.dat_I.df, scale = TRUE, center = TRUE)
    scores_pca1 = as.data.frame(pca.s$x)
    # Plot the first 2 principal components and convex hull
    rand.tr <-
      tri.mesh(scores_pca1[, 1], scores_pca1[, 2], "remove") # Delaunay triangulation
    rand.ch <- convex.hull(rand.tr, plot.it = F) # convex hull
    pr_poly <-
      cbind(x = c(rand.ch$x), y = c(rand.ch$y)) # save the convex hull vertices
    # PCA projection of study area population onto the principal components
    PCA_projection <-
      predict(pca.s, cov.dat.df) # Project study area population onto sample PC
    newScores = cbind(x = PCA_projection[, 1], y = PCA_projection[, 2]) # PC scores of projected population
    # Check which points fall within the polygon
    pip <-
      point.in.polygon(newScores[, 2], newScores[, 1], pr_poly[, 2], pr_poly[, 1], mode.checked =
                         FALSE)
    newScores <- data.frame(cbind(newScores, pip))
    klo_samples <- c(klo_samples, klo)
    prop_explained <-
      c(prop_explained, sum(newScores$pip) / nrow(newScores) * 100)
    number_of_samples <- c(number_of_samples, trial)
    print(
      paste(
        "N samples = ",
        trial,
        " out of ",
        final.n,
        "; iteration = ",
        iteration,
        "; KL = ",
        klo,
        "; Proportion = ",
        sum(newScores$pip) / nrow(newScores) * 100
      )
    )
  }
}


# 4 - Plot covariate diversity as PCA scores ===================================
plot(newScores[,1:2],
     xlab = "PCA 1",
     ylab = "PCA 2",
     xlim = c(min(newScores[,1:2], na.rm = T), max(newScores[,1:2], na.rm = T)),
     ylim = c(min(newScores[,1:2], na.rm = T), max(newScores[,1:2], na.rm = T)),
     col = 'black',
     main = 'Environmental space plots on convex hull of soil samples')

polygon(pr_poly, col = '#99999990')
# # Plot points outside convex hull
points(newScores[which(newScores$pip == 0), 1:2],
       col = 'red',
       pch = 12,
       cex = 1)

# 5 - KL divergence and % coincidence for growing N samples =============
# Merge data from number of samples, KL divergence and % coincidence
results <- data.frame(number_of_samples, klo_samples, prop_explained)
names(results) <- c("N", "KL", "Perc")

# Calculate mean results by N size
mean_result <- results %>%
  group_by(N) %>%
  summarize_all(mean)
mean_result

## Plot dispersion on KL and % by N
par(mar = c(5, 4, 5, 5))
boxplot(
  Perc ~ N,
  data = results,
  col = rgb(1, 0.1, 0, alpha = 0.5),
  ylab = "% coincidence"
)
mtext("KL divergence", side = 4, line = 3)
# Add new plot
par(new = TRUE, mar = c(5, 4, 5, 5))
# Box plot
boxplot(
  KL ~ N,
  data = results,
  axes = FALSE,
  outline = FALSE,
  col = rgb(0, 0.8, 1, alpha = 0.5),
  ylab = ""
)
axis(
  4,
  at = seq(0.02, 0.36, by = .06),
  label = seq(0.02, 0.36, by = .06),
  las = 3
)
# Draw legend
par(xpd=TRUE)
legend("top", inset=c(1,-.15) ,c("% coincidence", "KL divergence"), horiz=T,cex=.9,
       box.lty=0,fill = c(rgb(1, 0.1, 0, alpha = 0.5), rgb(0, 0.8, 1, alpha = 0.5)))
par(xpd=FALSE)

# 6 - Model KL divergence ======================================================
# Create an exponential decay function of the KL divergence
x <- mean_result$N
y = (mean_result$KL - min(mean_result$KL)) / (max(mean_result$KL) - min(mean_result$KL)) #KL

# Parametrize Exponential decay function
start <-  list()     # Initialize an empty list for the starting values

#fit function
k = 2
b0 = 0.01
b1 = 0.01

fit1 <-
  nls(
    y ~ k * exp(-b1 * x) + b0,
    start = list(k = k, b0 = b0, b1 = b1),
    control = list(maxiter = 500),
    trace = T
  )
summary(fit1)

# Plot fit
xx <- seq(1, final.n, 1)
plot(x, y)
lines(xx, predict(fit1, list(x = xx)))
# Predict with vfit function
jj <- predict(fit1, list(x = xx))
normalized = 1 - (jj - min(jj)) / (max(jj) - min(jj))


# 7 - Determine the minimum sample size to account for 95% of cumulative probability of the covariate diversity ====
minimum_n <- length(which(normalized < 0.95)) + 1

# Plot cdf and minimum sampling point
x <- xx
y <- normalized

mydata <- data.frame(x, y)
opti <- mydata[mydata$x == minimum_n, ]

plot_ly(
  mydata,
  x = ~ x,
  y = ~ normalized,
  mode = "lines+markers",
  type = "scatter",
  name = "CDF (1–KL divergence)"
) %>%
  add_trace(
    x = ~ x,
    y = ~ jj,
    mode = "lines+markers",
    type = "scatter",
    yaxis = "y2",
    name = "KL divergence"
  )  %>%
  add_trace(
    x = ~ opti$x,
    y = ~ opti$y,
    yaxis = "y",
    mode = "markers",
    name = "Minimum N",
    marker = list(
      size = 8,
      color = '#d62728',
      line = list(color = 'black', width = 1)
    )
  ) %>%
  layout(
    xaxis = list(
      title = "N",
      showgrid = T,
      dtick = 50,
      tickfont = list(size = 11)
    ),
    yaxis = list(title = "1–KL divergence (% CDF)", showgrid = F),
    yaxis2 = list(
      title = "KL divergence",
      overlaying = "y",
      side = "right"
    ),
    legend = list(
      orientation = "h",
      y = 1.2,
      x = 0.1,
      traceorder = "normal"
    ),
    margin = list(
      t = 50,
      b = 50,
      r = 100,
      l = 80
    ),
    hovermode = 'x'
  )  %>%
  config(displayModeBar = FALSE)

# 8 - Determine the optimal iteration according to the minimum N size ==========
optimal_iteration <- results[which(abs(results$N - minimum_n) == min(abs(results$N - minimum_n))), ] %>%
  mutate(IDX = 1:n()) %>%
  filter(Perc == max(Perc))

# 9 - Plot minimum points from best iteration ==================================
N_final <- samples_storage[paste0("N", optimal_iteration$N, "_", optimal_iteration$IDX)][[1]]
plot(cov.dat[[1]])
plot(N_final,add=T,col="red")

## 10 - Calculate COOBS (sgsR) =================================================
#' COOBS (count of observations) is an algorithm to map relatively
#' adequate and under-sampled areas on a sampling pattern.
#' COOBS algorithms allow one to understand which areas in a spatial domain are adequately and under-sampled. 
#' COOBS ca be used to design an additional survey by limiting the cLHS algorithm to areas where the COOBS
#' value is below some specified threshold. 
coobs <- calculate_coobs(
  mraster = cov.dat,
  existing = N_final,
  plot = TRUE,
  cores = 4,
  filename = paste0(results.path,"/coobs_SampleSize.tif"),
  overwrite = TRUE
)

# Plot COOBS

plot(coobs[[2]], main="Number of samples in similar environment")
plot(nghe[1], col="transparent", add=TRUE)

# coobs.results <-  rast(paste0(results.path,"coobs_SampleSize.tif"))
# plot(coobs.results[[2]], main="Number of samples in similar environment (COOBS)")
# plot(nghe[1], col="transparent", add=TRUE)

## 11 - Calculate minimum and and optimal sample size with opendsm =============

#install.packages('~/Downloads/DescTools_0.99.45.tar', repos=NULL, type='source')
#install.packages('DescTools')
library("DescTools")

#' Source scripts from  https://github.com/newdale/opendsm/tree/main
#'
#' Determine minimum sample size for the clhs algorithm
#' This script is based on the publication:
#' Determining minimum sample size for the conditioned Latin hypercube sampling algorithm
#' DOI: https://doi.org/10.1016/j.pedsph.2022.09.001
source("scripts/clhs_min.R")  
#' Optimal sample size based on the publication:
#' Divergence metrics for determining optimal training sample size in digital soil mapping
#' DOI: https://doi.org/10.1016/j.geoderma.2023.116553   
source("scripts/opt_sample.R")
#' Perform Sequential Variable Inflation Factor Analysis
source("scripts/seqVIF.R")
#' Function to calculate the KL-divergence
#' score between two probability distributions, P and Q.
source("scripts/kldiv.R")
#' Function to calculate the Jensen-Shannon Distance
#' score between two probability distributions, P and Q.
source("scripts/jsdist.R")
#' Function to calculate the Jensen-Shannon-divergence
#' score between two probability distributions, P and Q.
source("scripts/jsdiv.R")

# Convert raster covariates to dataframe
cov.dat.df <- data.frame(cov.dat)
# Calculate minimum sample size
min_N <-clhs_min(cov.dat.df)
# Calculate optimal sample size using normalized KL-div, JS-div and JS distance
opt_N <- opt_sample(alg="clhs",s_min=150,s_max=500, s_step=50, s_reps=4, covs = cov.dat.df,clhs_iter=100, cpus=NULL, conf=0.95)


## END


