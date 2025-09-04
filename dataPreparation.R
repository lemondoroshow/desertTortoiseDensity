library(dplyr)

############################ SPLINE INTERPOLATION ##############################

# Import data
raw_densities <- read.csv('./data/density/density_data_raw.csv')
strata <- colnames(raw_densities)[colnames(raw_densities) != 'Year']

# Set up adjusted density data frame
adj_densities <- data.frame(matrix(
  nrow = dim(raw_densities)[1],
  ncol = dim(raw_densities)[2]
))
colnames(adj_densities) <- colnames(raw_densities)
adj_densities$Year <- raw_densities$Year

# Prepare to count missing years
num_na <- vector(mode = "list", length = length(strata))
names(num_na) <- strata

# Iterate through strata
for (i in 2:(dim(raw_densities)[2])) {
  
  # Filter data
  ts <- raw_densities[,c(1, i)]
  stratum <- colnames(ts)[2]
  years <- as.matrix(ts)[,1]
  densities <- as.matrix(ts)[,2]
  
  # Count how many missing years there are
  num_na[[stratum]] <- sum(is.na(densities))
  
  # Fit spline function
  fit <- spline(x = years, y = densities, xout = years)$y
  
  # Remove negative values from interpolation
  fit[fit < 0] <- 0
  
  # Add results to dataframe
  adj_densities[stratum] <- lapply(fit, round, 1) |>
    unlist()
}

# Clean up
rm(i, ts, stratum, years, densities, fit)