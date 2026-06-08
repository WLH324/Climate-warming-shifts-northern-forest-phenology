# Arctic polar-stereographic map of regression coefficients with hard zero break.

suppressPackageStartupMessages({
  library(terra)
  library(tidyterra)
  library(sf)
  library(ggplot2)
  library(rnaturalearth)
  library(scales)
  source("R/config.R")
})

plot_polar_coefficient <- function(coef_raster, output_png = NULL,
                                   lat_min = 30, fill_limits = c(-1, 1),
                                   title = NULL) {
  geo_crs <- "EPSG:4326"
  polar_crs <- "EPSG:3995"

  lat_extent <- as.polygons(ext(-180, 180, lat_min, 90), crs = geo_crs)
  ratio_geo <- crop(coef_raster, lat_extent)
  ratio_proj <- project(ratio_geo, polar_crs)

  pole <- st_sfc(st_point(c(0, 0)), crs = polar_crs)
  polar_circle <- st_buffer(pole, dist = 7e6)
  ratio_proj <- mask(ratio_proj, vect(polar_circle))

  world <- ne_countries(scale = "medium", returnclass = "sf")
  world_clip <- st_intersection(st_transform(world, polar_crs), polar_circle)
  grat_clip <- st_intersection(
    st_transform(st_graticule(lat = c(40, 50, 60, 70, 80), lon = seq(-180, 180, 30)), polar_crs),
    polar_circle
  )
  outer_circle_sf <- st_sf(geometry = polar_circle)
  lim <- 7.05e6

  trend_colors <- c("#008080", "#66c2a5", "#b3e2e8", "#e0f3f8",
                    "#fc8d59", "#e34a33", "#b30000", "#7f0000")

  p <- ggplot() +
    geom_spatraster(data = ratio_proj, aes(fill = after_stat(value))) +
    scale_fill_gradientn(
      colours = trend_colors,
      values = c(0, 0.25, 0.375, 0.499, 0.501, 0.625, 0.75, 1),
      limits = fill_limits,
      oob = oob_squish,
      breaks = c(-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1),
      labels = c("< -1", "-0.8", "-0.6", "-0.4", "-0.2", "0",
                 "0.2", "0.4", "0.6", "0.8", "> 1"),
      na.value = "transparent"
    ) +
    geom_sf(data = grat_clip, colour = "grey65", linewidth = 0.25, linetype = "dashed") +
    geom_sf(data = outer_circle_sf, fill = NA, colour = "black", linewidth = 1.1) +
    geom_sf(data = world_clip, fill = NA, colour = "black", linewidth = 0.4) +
    coord_sf(crs = st_crs(polar_crs), xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.key.height = unit(3.5, "cm"),
      legend.key.width = unit(0.4, "cm"),
      plot.title = element_text(size = 14, hjust = 0.5)
    )
  if (!is.null(title)) p <- p + labs(title = title)

  if (!is.null(output_png)) {
    dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)
    ggsave(output_png, p, width = 8, height = 8, dpi = 500, bg = "white")
  }
  invisible(p)
}

if (sys.nframe() == 0L) {
  root <- get_data_root()
  coef_path <- file.path(root, "results", "diffreg_l1", "SSP1_2.6",
                         "PlantedForest", "Coef_L3prev_to_L1.tif")
  out_png <- file.path(root, "results", "figures", "Coef_L3prev_to_L1_polar.png")
  plot_polar_coefficient(rast(coef_path), out_png,
                         title = "L3 carry-over coefficient on L1 (SSP1-2.6)")
}
