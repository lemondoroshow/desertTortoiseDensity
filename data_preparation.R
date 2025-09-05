library(dplyr)
library(ggplot2)
library(gridExtra)

############################ SPLINE INTERPOLATION ##############################

# Import data
raw_densities <- read.csv('./data/density/densitiesRaw.csv')
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
  ggsave(paste0('./figures/spline/splineComparison', stratum, '.png'), 
         combined_plot, width = 8, height = 5, units = 'in')
}

# Clean up
rm(i, ts, stratum, years, densities, fit, adj_df, adj_plot, combined_plot,
   raw_df, raw_plot, colors, lower, upper, where_na)
num_na <- unlist(num_na)
mean(num_na) # 8.5 missing years on average
mean(num_na[names(num_na) != 'AG']) # 8.93 missing years on average w/o AG

# Export data
write.csv(adj_densities, './data/density/densitiesAdj.csv', row.names = FALSE, quote = FALSE)
