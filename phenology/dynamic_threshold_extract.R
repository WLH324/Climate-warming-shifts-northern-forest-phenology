# Dynamic-threshold phenology from daily GOSIF (20% amplitude SOS/EOS; 90% peak PPOS/APOS).

suppressPackageStartupMessages({
  library(lubridate)
  library(terra)
  source("R/config.R")
  source("R/phenology_io.R")
})

extract_dynamic_threshold_year <- function(year_matrix, year_days_idx, daily_dates,
                                         ratio_green = 0.20, ratio_plat = 0.90) {
  n_pixels <- nrow(year_matrix)
  sos <- ppos <- apos <- eos <- rep(NA_integer_, n_pixels)

  annual_max <- apply(year_matrix, 1, max, na.rm = TRUE)
  annual_min <- apply(year_matrix, 1, min, na.rm = TRUE)
  invalid <- is.na(annual_max) | is.na(annual_min) | annual_max <= 0.1 | (annual_max - annual_min) < 0.05
  range_val <- annual_max - annual_min
  range_val[invalid | range_val == 0] <- NA

  ratio_amp <- (year_matrix - annual_min) / range_val
  ratio_max <- year_matrix / annual_max
  ratio_max[invalid, ] <- NA

  sos_idx <- apply(ratio_amp >= ratio_green, 1, function(x) { h <- which(x); if (length(h)) h[1] else NA })
  eos_idx <- apply(ratio_amp >= ratio_green, 1, function(x) { h <- which(x); if (length(h)) rev(h)[1] else NA })
  ppos_idx <- apply(ratio_max >= ratio_plat, 1, function(x) { h <- which(x); if (length(h)) h[1] else NA })
  apos_idx <- apply(ratio_max >= ratio_plat, 1, function(x) { h <- which(x); if (length(h)) rev(h)[1] else NA })

  SOS_global <- year_days_idx[sos_idx]
  EOS_global <- year_days_idx[eos_idx]
  PPOS_global <- year_days_idx[ppos_idx]
  APOS_global <- year_days_idx[apos_idx]

  doy_seq <- yday(daily_dates[year_days_idx])
  SOS_doy <- doy_seq[sos_idx]
  EOS_doy <- doy_seq[eos_idx]
  PPOS_doy <- doy_seq[ppos_idx]
  APOS_doy <- doy_seq[apos_idx]

  valid <- !is.na(SOS_doy) & SOS_doy >= 70 & SOS_doy <= 200 &
    !is.na(EOS_doy) & EOS_doy >= 220 & EOS_doy <= 310 &
    !is.na(SOS_global) & !is.na(PPOS_global) & !is.na(APOS_global) & !is.na(EOS_global) &
    SOS_global < PPOS_global & PPOS_global <= APOS_global & APOS_global < EOS_global

  list(
    SOS_DOY = ifelse(valid, as.integer(SOS_doy), NA),
    PPOS_DOY = ifelse(valid, as.integer(PPOS_doy), NA),
    APOS_DOY = ifelse(valid, as.integer(APOS_doy), NA),
    EOS_DOY = ifelse(valid, as.integer(EOS_doy), NA),
    Greenup_Duration_days = ifelse(valid, PPOS_global - SOS_global + 1, NA),
    Plateau_Duration_days = ifelse(valid, APOS_global - PPOS_global + 1, NA),
    Senescence_Duration_days = ifelse(valid, EOS_global - APOS_global + 1, NA)
  )
}

run_dynamic_threshold <- function(input_dir, output_dir, years = YEARS_HIST,
                                  daily_dates = seq(ymd("2000-02-26"), ymd("2024-12-31"), by = "day")) {
  dat <- load_daily_stack(input_dir)
  template <- dat$stack[[1]]
  n_pixels <- nrow(dat$matrix)

  init_mat <- function() matrix(NA_integer_, nrow = n_pixels, ncol = length(years))
  layers <- list(
    SOS_DOY = init_mat(), PPOS_DOY = init_mat(), APOS_DOY = init_mat(), EOS_DOY = init_mat(),
    Greenup_Duration_days = init_mat(), Plateau_Duration_days = init_mat(),
    Senescence_Duration_days = init_mat()
  )

  for (i in seq_along(years)) {
    yr <- years[i]
    idx <- which(year(daily_dates) == yr)
    if (!length(idx)) next
    yr_out <- extract_dynamic_threshold_year(dat$matrix[, idx, drop = FALSE], idx, daily_dates)
    for (nm in names(layers)) layers[[nm]][, i] <- yr_out[[nm]]
  }

  save_phenology_layers(template, output_dir, years, layers)
}

plot_region_diagnostic <- function(input_dir, layers, years, daily_dates, region_ext) {
  dat <- load_daily_stack(input_dir)
  template <- rast(dat$files[1])
  region_cells <- cells(crop(template, region_ext), 1)[[1]]
  if (!length(region_cells)) stop("No pixels in region")

  yearly_mean <- matrix(NA_real_, 366, length(years))
  colnames(yearly_mean) <- years
  region_daily <- dat$matrix[region_cells, , drop = FALSE]

  for (i in seq_along(years)) {
    idx <- which(year(daily_dates) == years[i])
    if (!length(idx)) next
    temp <- colMeans(region_daily[, idx, drop = FALSE], na.rm = TRUE)
    yearly_mean[yday(daily_dates[idx]), i] <- temp
  }
  annual_cycle <- rowMeans(yearly_mean, na.rm = TRUE)

  med_doy <- function(mat) {
    median(apply(mat[region_cells, , drop = FALSE], 2,
                 function(x) median(yday(daily_dates[x]), na.rm = TRUE)), na.rm = TRUE)
  }

  dev.new(width = 14, height = 6)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  plot(1:366, annual_cycle, type = "n", xlab = "DoY", ylab = "GOSIF", main = "25-year curves")
  for (i in seq_along(years)) lines(1:366, yearly_mean[, i], col = "#33333344", lwd = 1.5)
  lines(1:366, annual_cycle, col = "red", lwd = 3)
  plot(1:366, annual_cycle, type = "l", lwd = 3, xlab = "DoY", ylab = "GOSIF", main = "Mean cycle + phenology")
  for (lab in c("SOS_DOY", "PPOS_DOY", "APOS_DOY", "EOS_DOY")) {
    d <- med_doy(layers[[lab]])
    if (!is.na(d)) abline(v = d, lty = 2, col = "steelblue")
  }
  invisible(NULL)
}

if (sys.nframe() == 0L) {
  root <- get_data_root()
  run_dynamic_threshold(
    input_dir = file.path(root, "gosif_sg", "region_01"),
    output_dir = file.path(root, "results", "phenology_dynamic")
  )
}
