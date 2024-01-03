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
  # while retaining for a 95% of representativeness in the environmental variability of covariates
  # in the area 
# 
# 0 - Set working directory and load packages
# 1 - User-defined variables 
# 2 - Import national data 
# 3 - Calculate the minimum sample size to describe the area
# 4 - Plot covariate diversity as PCA scores 
# 5 - KL divergence and similarity results for growing N samples
# 6 - Model KL divergence
# 7 - Determine the minimum sample size for 95% representativeness
# 8 - Determine the optimal iteration according to the minimum N size 
# 9 - Plot minimum points from best iteration
#________________________________________________________________

 start_time <- Sys.time()

## 0 - Set working directory and load packages =================================

  # Set working directory to source file location
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  setwd("../") # Move wd down to main folder
  
  # List of packages
  packages <- c("sp","terra","raster","sf","clhs", "sgsR","entropy", "tripack",
      "manipulate","dplyr","plotly","synoptReg")
  
  # Load packages
  lapply(packages, require, character.only = TRUE) 
  rm(packages) # Remove object to save memory space


## 1 - User-defined variables ==================================================
# Path to rasters
  raster.path <- "data/vnm/rasters"
# Path to shapes
  shp.path <- "data/vnm/shapes"
# Path to results
  results.path <- "data/vnm/results/"
# Aggregation factor for up-scaling raster covariates (optional)
  agg.factor = 10

# Define parameters to determine minimum sampling size
  initial.n <- 50 # Initial sampling size to test
  final.n <- 300 # Final sampling size to test
  by.n <- 10 # Increment size
  iters <- 5 # Number of trials on each size
  

## 2 - Import national data ====================================================
  # Load raster data
  cov.dat <- list.files(
      raster.path, pattern = "tif$",
      recursive = TRUE, full.names = TRUE)
  cov.dat <- terra::rast(cov.dat)
  # Drop lights at night and Population density
  cov.dat <- cov.dat[[-c(70:71)]]

  # Aggregate raster pixels to simplify calculations (optional)
  cov.dat <- aggregate(cov.dat, fact = agg.factor, fun = "mean")

  # Load district
  nghe <- sf::st_read(file.path(paste0(shp.path,"/Nghe_An.shp")),quiet=TRUE)
  
  # Crop covariates on administrative boundary
  cov.dat <- crop(cov.dat, nghe, mask=TRUE)
  
  # Simplify raster information with PCA
  pca <- raster_pca(cov.dat)
  
  # Get SpatRaster layers
  cov.dat <- pca$PCA
  # Create a raster stack to be used as input in the clhs::clhs function 
  cov.dat.ras <- raster::stack(cov.dat) 
  # Subset rasters to first components
  break.pca <- first(which(pca$summaryPCA[3,]>0.99))
  cov.dat <- pca$PCA[[1:break.pca]]
  cov.dat.ras <-  cov.dat.ras[[1:break.pca]]
  
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
      p.dat_I <-  clhs(cov.dat.ras,
          size = trial, iter = 5000,
          progress = FALSE, simple = FALSE)
      
      # Get covariate values for each point
      p.dat_I <- p.dat_I$sampled_data
      # Get the covariate values at points as dataframe and delete NAs
      p.dat_I.df <- as.data.frame(p.dat_I@data) %>%
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
      # q.mat
      
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
  plot(newScores[, 1:2],
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

# 5 - KL divergence and % representativeness for growing N samples =============
  # Merge data from number of samples, KL divergence and % representativeness
  results <- data.frame(number_of_samples, klo_samples, prop_explained)
  names(results) <- c("N", "KL", "Perc")

  # Calculate mean results by N size
    mean_result <- results %>%
      group_by(N) %>%
      summarize_all(mean)
    mean_result

  ## Plot dispersion on KL and % by N
    par(mar = c(5, 4, 1, 6))
    boxplot(
      Perc ~ N,
      data = results,
      col = rgb(1, 0.1, 0, alpha = 0.5),
      ylab = "%"
    )
    mtext("KL divergence", side = 4, line = 3)
  # Add new plot
    par(new = TRUE, mar = c(5, 4, 1, 6))
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
  points(N_final)
  
####  Time required to compute the script
  Sys.time() - start_time 
  
  
