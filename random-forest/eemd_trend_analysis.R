# EEMD nonlinear trend analysis on future RF-predicted phenology by latitude band.

suppressPackageStartupMessages({
  library(Rlibeemd)
  library(terra)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(zoo)
  library(patchwork)
  source("R/config.R")
})

process_eemd <- function(ts, years, ensemble = 100, noise = 0.2) {
  ts <- na.approx(ts, na.rm = FALSE)
  imfs <- eemd(ts, ensemble_size = ensemble, noise = noise)
  residual <- imfs[, ncol(imfs)]
  data.frame(Year = years, rel_trend = residual - residual[1])
}

run_eemd_by_latband <- function(pred_root, ssps, forest_types = FOREST_TYPES,
                                phenos = c("L_total", "L1", "L2", "L3"),
                                years = YEARS_FUTURE) {
  lat_bands <- c("30-40N", "40-50N", "50-60N", ">60N")
  lat_breaks <- c(30, 40, 50, 60, 90)
  results <- list()

  sample <- list.files(file.path(pred_root, ssps[1], "future_phenology", forest_types[1]),
                       pattern = paste0(phenos[1], "_", years[1], ".tif"), full.names = TRUE)[1]
  zones <- classify(init(rast(sample), "y"),
                    rcl = cbind(lat_breaks[-length(lat_breaks)], lat_breaks[-1], seq_along(lat_bands)))

  for (ssp in ssps) {
    for (ft in forest_types) {
      for (ph in phenos) {
        zone_ts <- setNames(vector("list", length(lat_bands)), lat_bands)
        for (i in seq_along(years)) {
          f <- file.path(pred_root, ssp, "future_phenology", ft, paste0(ph, "_", years[i], ".tif"))
          if (!file.exists(f)) next
          zm <- zonal(rast(f), zones, "mean", na.rm = TRUE)
          for (z in seq_along(lat_bands)) {
            if (is.null(zone_ts[[z]])) zone_ts[[z]] <- rep(NA, length(years))
            zone_ts[[z]][i] <- zm[z, 2]
          }
        }
        for (z in seq_along(lat_bands)) {
          ts <- zone_ts[[z]]
          if (all(is.na(ts))) next
          res <- process_eemd(ts, years)
          res$SSP <- ssp; res$Forest <- ft; res$Pheno <- ph; res$LatBand <- lat_bands[z]
          results[[length(results) + 1]] <- res
        }
      }
    }
  }
  bind_rows(results)
}

if (sys.nframe() == 0L) {
  root <- get_data_root()
  df <- run_eemd_by_latband(file.path(root, "cmip6"), c("SSP1_2.6", "SSP2_4.5", "SSP5_8.5"))
  out <- file.path(root, "results", "eemd_trends.csv")
  write.csv(df, out, row.names = FALSE)
  message("Saved: ", out)
}
