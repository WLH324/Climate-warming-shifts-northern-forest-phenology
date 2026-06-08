# Hovmoller diagrams: latitude x year composites of phenology ratios.

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  source("R/config.R")
})

lat_mean_matrix <- function(r_stack, lat_min = 30, lat_max = 90) {
  lat_r <- init(r_stack[[1]], "y")
  lat_seq <- seq(lat_min, lat_max, by = 1)
  mat <- matrix(NA_real_, nrow = length(lat_seq), ncol = nlyr(r_stack))
  for (i in seq_along(lat_seq)) {
    band <- lat_r >= lat_seq[i] - 0.5 & lat_r < lat_seq[i] + 0.5
    for (t in seq_len(nlyr(r_stack))) {
      v <- values(r_stack[[t]])[values(band)]
      mat[i, t] <- mean(v, na.rm = TRUE)
    }
  }
  list(lat = lat_seq, mat = mat)
}

plot_hovmoller <- function(mat, years, title, ylab = "Latitude (deg N)") {
  df <- as.data.frame(mat$mat)
  colnames(df) <- years
  df$lat <- mat$lat
  long <- tidyr::pivot_longer(df, -lat, names_to = "Year", values_to = "Value")
  long$Year <- as.integer(long$Year)
  ggplot(long, aes(Year, lat, fill = Value)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma") +
    labs(title = title, x = "Year", y = ylab, fill = NULL) +
    theme_minimal()
}

if (sys.nframe() == 0L) {
  root <- get_data_root()
  pheno <- file.path(root, "phenology")
  years <- YEARS_HIST
  vgu <- rast(file.path(pheno, "Greenup_Duration", paste0("Greenup_Duration_days_", years, "_PlantedForest.tif")))
  vss <- rast(file.path(pheno, "Senescence_Duration", paste0("Senescence_Duration_days_", years, "_PlantedForest.tif")))
  ratio <- vgu / vss
  p <- plot_hovmoller(lat_mean_matrix(ratio), years, "VGU/VSS ratio (planted forest)")
  ggsave(file.path(root, "results", "hovmoller_vgu_vss.png"), p, width = 8, height = 5, dpi = 300)
}
