library(ggplot2)
library(gridExtra)

############################ SPLINE INTERPOLATION ##############################

# Import data
raw_densities <- read.csv('./data/density/densities_raw.csv')
adj_densities <- read.csv('./data/density/densities_spline_adj.csv')

# Iterate through data
for (i in 2:(dim(raw_densities)[2])) {
  
  # Filter data
  ts_raw <- raw_densities[,c(1, i)]
  ts_adj <- adj_densities[,c(1, i)]
  stratum <- colnames(ts_raw)[2]
  years <- as.matrix(ts_raw)[,1]
  densities <- as.matrix(ts_raw)[,2]
  fit <- as.matrix(ts_adj)[,2]
  
  print(stratum) # Debugging
  
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

############################ PRISM PRECIPITATION ###############################

# Import data
densities <- read.csv('./data/density/densities_raw.csv')
precips <- read.csv('./data/compiledCovariates/annual_precipitation.csv')

# Iterate through data
for (i in 2:dim(densities)[2]) {
  
  # Open data
  dense <- as.matrix(densities[,c(1, i)])
  precip <- as.matrix(precips[, c(1, i)])
  stratum <- colnames(dense)[2]
  years <- as.matrix(dense)[,1]
  print(stratum) # Debugging
  
  # Scale densities and precipitation for axes' sake
  scl <- mean(precip[, stratum]) / 
    mean(dense[, stratum], na.rm = TRUE)
  
  # Plot densities and precipitation together
  df <- data.frame(year = years, density = dense[, stratum], precip = precip[, stratum]) |>
    dplyr::filter(!is.na(density))
  colors <- c('Density' = 'darkgreen', 'Precipitation' = 'royalblue')
  plot <- ggplot(data = df, aes(x = year)) +
    geom_line(aes(y = density, color = 'Density'), linewidth = 0.7) +
    geom_line(aes(y = precip / scl, color = 'Precipitation'), linewidth = 0.7) +
    scale_y_continuous(name = "Density (torts / km²)",
                       sec.axis = sec_axis(trans = ~ . * scl, name = "Precipitation (mm / yr)")) +
    scale_color_manual(values = colors) +
    labs(x = "Year", title = paste0("Comparison of tortoise densities and annual precipitation for ", stratum)) +
    theme(panel.grid = element_line(color = 'grey'),
          panel.background = element_rect(fill = 'white', colour = 'black')) +
    scale_x_continuous(limits = c(2000, 2025))
  
  # Export plots
  ggsave(paste0('./figures/densityPrecip/density_precip_', stratum, '.png'), plot)
}
