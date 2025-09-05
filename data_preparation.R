library(dplyr)
library(ggplot2)
library(gridExtra)
library(prism)
library(sf)
library(terra)
library(tidyterra)
library(tmap)

############################ SPLINE INTERPOLATION ##############################

# Import data
raw_densities <- read.csv('./data/density/densities_raw.csv')
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
  where_na <- is.na(densities)
  num_na[[stratum]] <- sum(where_na)
  
  # Fit spline function
  fit <- spline(x = years, y = densities, xout = years)$y
  
  # Remove negative values from interpolation and remove outliers from spline
  fit[fit < 0] <- 0
  lower <- quantile(densities, 0.25, na.rm = TRUE) - 1.5 * IQR(densities, na.rm = TRUE)
  upper <- quantile(densities, 0.75, na.rm = TRUE) + 1.5 * IQR(densities, na.rm = TRUE)
  fit[where_na & fit < lower] <- lower
  fit[where_na & fit > upper] <- upper
  
  # Add results to dataframe
  adj_densities[stratum] <- lapply(fit, round, 1) |>
    unlist()
  
  # Plot raw densities
  raw_df <- data.frame(years = years, densities = densities) |>
    dplyr::filter(!is.na(densities))
  colors <- c('Raw density line' = 'darkgreen', 'Raw density points' = 'black')
  raw_plot <- ggplot(data = raw_df, aes(x = years)) +
    geom_line(aes(y = densities, color = 'Raw density line'), linewidth = 0.7) + 
    geom_point(aes(y = densities, color = 'Raw density points')) +
    labs(x = 'Year', y = 'Density (torts / km²)', color = 'Color',
         title = paste0('Raw densities across ', stratum, ' from range-wide monitoring, 2001 - 2024')) +
    theme(panel.grid = element_line(color = 'grey'),
          panel.background = element_rect(fill = 'white', colour = 'black')) +
    scale_color_manual(values = colors) +
    scale_x_continuous(limits = c(2000, 2025))
  
  # Plot adjusted densities
  adj_df <- data.frame(years = years, densities = fit, points = densities)
  colors <- c('Spline density line' = 'purple', 'Raw density points' = 'black')
  adj_plot <- ggplot(data = adj_df, aes(x = years)) +
    geom_line(aes(y = densities, color = 'Spline density line'), linewidth = 0.7) + 
    geom_point(aes(y = points, color = 'Raw density points')) +
    labs(x = 'Year', y = 'Density (torts / km²)', color = 'Color',
         title = paste0('Adjusted densities across ', stratum, ' from range-wide monitoring, 2001 - 2024')) +
    theme(panel.grid = element_line(color = 'grey'),
          panel.background = element_rect(fill = 'white', colour = 'black')) +
    scale_color_manual(values = colors) +
    scale_x_continuous(limits = c(2000, 2025))
  
  # Put both plots together
  combined_plot <- grid.arrange(raw_plot, adj_plot, ncol = 1)
  
  # Export
  ggsave(paste0('./figures/spline/spline_comparison_', stratum, '.png'), 
         combined_plot, width = 8, height = 5, units = 'in')
}

# Clean up
rm(i, ts, stratum, years, densities, fit, adj_df, adj_plot, combined_plot,
   raw_df, raw_plot, colors, lower, upper, where_na)
num_na <- unlist(num_na)
mean(num_na) # 8.5 missing years on average
mean(num_na[names(num_na) != 'AG']) # 8.93 missing years on average w/o AG

# Export data
write.csv(adj_densities, './data/density/densities_spline_adj.csv', row.names = FALSE, quote = FALSE)

############################# PRISM ACQUISITION ################################

# Set download directory
prism_set_dl_dir('./data/prism')

# Get annual PRISM data
get_prism_annual(
  type = 'ppt',
  years = 2001:2024,
  keepZip = FALSE
)

# List all shapefiles
shapes <- c(
  list.files(path = './tcaShapes/coloradoDesert/', 
             pattern = '.shp', recursive = TRUE, full.names = TRUE),
  list.files(path = './tcaShapes/easternMojave/', 
             pattern = '.shp', recursive = TRUE, full.names = TRUE),
  list.files(path = './tcaShapes/northeasternMojave/', 
             pattern = '.shp', recursive = TRUE, full.names = TRUE),
  list.files(path = './tcaShapes/westernMojave/', 
             pattern = '.shp', recursive = TRUE, full.names = TRUE)
)

# Set up data frame
total_data <- data.frame(matrix(ncol = 0, nrow = 24))
total_data$Year <- 2001:2024

# Iterate through strata
for (shp in shapes) {
  stratum <- strata[which(shapes == shp)]
  
  # Open shapefile
  shp <- sf::read_sf(shp) |>
    terra::vect()
  
  # Iterate through years of interest
  precips <- list()
  for (year in 2001:2024) {
    print(paste0("Beginning ", stratum, " for ", year)) # Debugging
    
    # Open PRISM data
    prism_rast <- prism_archive_subset(
      type = 'ppt',
      temp_period = 'annual',
      years = year
    ) |>
      pd_to_file() |>
      terra::rast()
    
    # Clip raster to stratum, calculate mean
    mean_precip <- terra::project(prism_rast, shp) |>
      terra::mask(shp) |>
      terra::crop(shp) |>
      terra::global(mean, na.rm = TRUE)
    
    # Add to list
    precips <- c(precips, mean_precip$mean)
  }
  
  # Add to data frame
  stratum_data = data.frame(
    Year = 2001:2024,
    placeholder = unlist(precips)
  ) |> rename(!!stratum := placeholder)
  total_data <- left_join(total_data, stratum_data, by = 'Year')
}

# Clean up
rm(mean_precip, precips, prism_rast, shp, stratum_data, stratum, year)

# Export total data
write.csv(total_data, './data/compiledCovariates/annualPrecipitation.csv', row.names = FALSE, quote = FALSE)
