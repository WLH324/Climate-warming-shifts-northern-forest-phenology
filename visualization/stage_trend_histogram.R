# Histogram of Sen's slope distributions for phenological stages across scenarios.

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(dplyr)
  source("R/config.R")
})

read_slope_values <- function(tif_path) {
  if (!file.exists(tif_path)) return(numeric(0))
  as.numeric(values(rast(tif_path)))
}

plot_stage_histograms <- function(slope_specs, out_png, bins = 80) {
  df <- bind_rows(lapply(seq_len(nrow(slope_specs)), function(i) {
    data.frame(
      value = read_slope_values(slope_specs$path[i]),
      stage = slope_specs$stage[i],
      forest = slope_specs$forest[i],
      scenario = slope_specs$scenario[i]
    )
  })) %>% filter(is.finite(value))

  df$forest <- factor(df$forest, levels = c("NaturalForest", "PlantedForest"))
  df$stage <- factor(df$stage, levels = c("L_total", "L1", "L2", "L3"))

  p <- ggplot(df, aes(value, fill = forest)) +
    geom_histogram(bins = bins, position = "identity", alpha = 0.55) +
    facet_grid(forest ~ stage, scales = "free_y") +
    labs(x = "Sen's slope (days yr-1)", y = "Pixel count", title = "Stage duration trends") +
    theme_minimal()

  ggsave(out_png, p, width = 10, height = 8, dpi = 300)
  invisible(p)
}

if (sys.nframe() == 0L) {
  root <- get_data_root()
  specs <- expand.grid(
    stage = c("L_total", "L1", "L2", "L3"),
    forest = FOREST_TYPES,
    scenario = c("historical", "SSP2_4.5"),
    stringsAsFactors = FALSE
  )
  specs$path <- file.path(root, "results", "mk_sen", specs$scenario, specs$forest,
                          paste0(specs$stage, "_slope.tif"))
  plot_stage_histograms(specs, file.path(root, "results", "stage_trend_histogram.png"))
}
