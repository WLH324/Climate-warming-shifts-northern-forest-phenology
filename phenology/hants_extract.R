# HANTS-based phenology extraction from daily GOSIF time series.

suppressPackageStartupMessages({
  library(lubridate)
  library(terra)
  source("R/config.R")
  source("R/phenology_io.R")
})

hants_fit <- function(y, period = 365, n_harmonics = 4) {
  n <- length(y)
  t <- 0:(n - 1)
  A <- matrix(0, nrow = n, ncol = 2 * n_harmonics + 1)
  A[, 1] <- 1
  for (i in seq_len(n_harmonics)) {
    A[, 2 * i] <- cos(2 * pi * i * t / period)
    A[, 2 * i + 1] <- sin(2 * pi * i * t / period)
  }
  as.numeric(A %*% solve(t(A) %*% A) %*% t(A) %*% y)
}

extract_phenology_hants <- function(daily_matrix, year_days_idx, n_harmonics = 4, period = 365) {
  n_pixels <- nrow(daily_matrix)
  sos <- ppos <- apos <- eos <- rep(NA_integer_, n_pixels)
  year_mat <- daily_matrix[, year_days_idx, drop = FALSE]
  annual_max <- apply(year_mat, 1, max, na.rm = TRUE)
  annual_min <- apply(year_mat, 1, min, na.rm = TRUE)
  invalid <- is.na(annual_max) | annual_max <= 0.1 | (annual_max - annual_min) < 0.05
  range_val <- annual_max - annual_min
  range_val[invalid | range_val == 0] <- NA
  ratio <- (year_mat - annual_min) / range_val

  for (p in seq_len(n_pixels)) {
    if (invalid[p]) next
    y <- ratio[p, ]
    if (sum(!is.na(y)) < 30) next
    fit <- hants_fit(y, period, n_harmonics)
    peak <- which.max(fit)
    thr_lo <- 0.2 * max(fit, na.rm = TRUE)
    thr_hi <- 0.9 * max(fit, na.rm = TRUE)
    sos[p] <- which(fit >= thr_lo)[1]
    ppos[p] <- peak
    apos[p] <- which(fit >= thr_hi)[1]
    eos[p] <- tail(which(fit >= thr_lo), 1)
  }
  list(sos = sos, ppos = ppos, apos = apos, eos = eos)
}

run_hants_batch <- function(input_dir, output_dir, years = YEARS_HIST,
                            dates = seq(ymd("2000-02-26"), ymd("2024-12-31"), by = "day")) {
  dat <- load_daily_stack(input_dir)
  n_pixels <- nrow(dat$matrix)
  init_mat <- function() matrix(NA_integer_, nrow = n_pixels, ncol = length(years))
  layers <- list(
    SOS_DOY = init_mat(), PPOS_DOY = init_mat(), APOS_DOY = init_mat(), EOS_DOY = init_mat(),
    Greenup_Duration_days = init_mat(), Plateau_Duration_days = init_mat(),
    Senescence_Duration_days = init_mat()
  )
  for (i in seq_along(years)) {
    idx <- which(year(dates) == years[i])
    if (!length(idx)) next
    doy_seq <- yday(dates[idx])
    ph <- extract_phenology_hants(dat$matrix, idx)
    ok <- !is.na(ph$sos) & !is.na(ph$eos) & !is.na(ph$ppos) & !is.na(ph$apos)
    layers$SOS_DOY[ok, i] <- doy_seq[ph$sos[ok]]
    layers$PPOS_DOY[ok, i] <- doy_seq[ph$ppos[ok]]
    layers$APOS_DOY[ok, i] <- doy_seq[ph$apos[ok]]
    layers$EOS_DOY[ok, i] <- doy_seq[ph$eos[ok]]
    layers$Greenup_Duration_days[ok, i] <- ph$ppos[ok] - ph$sos[ok] + 1
    layers$Plateau_Duration_days[ok, i] <- ph$apos[ok] - ph$ppos[ok] + 1
    layers$Senescence_Duration_days[ok, i] <- ph$eos[ok] - ph$apos[ok] + 1
  }
  save_phenology_layers(dat$stack[[1]], output_dir, years, layers)
}

if (sys.nframe() == 0L) {
  root <- get_data_root()
  run_hants_batch(file.path(root, "gosif_sg", "region_01"), file.path(root, "results", "phenology_hants"))
}
