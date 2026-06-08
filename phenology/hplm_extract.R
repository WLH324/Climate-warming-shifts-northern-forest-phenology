# HPLM phenology: Logistic green-up + Weibull senescence on normalized daily GOSIF.

suppressPackageStartupMessages({
  library(lubridate)
  library(terra)
  library(minpack.lm)
  source("R/config.R")
  source("R/phenology_io.R")
})

hplm_fit <- function(t, y, doy_seq) {
  if (sum(!is.na(y)) < 50) return(list(fitted = rep(NA, length(y)), peak_idx = NA))

  peak_candidates <- which(diff(sign(diff(y, na.rm = FALSE)), na.rm = FALSE) < 0) + 1
  valid_peaks <- peak_candidates[doy_seq[peak_candidates] >= 75 & doy_seq[peak_candidates] <= 315]
  if (length(valid_peaks)) {
    peak_idx <- valid_peaks[which.max(y[valid_peaks])]
  } else {
    g <- which.max(y)
    peak_idx <- if (doy_seq[g] >= 75 && doy_seq[g] <= 315) g else NA
  }
  if (is.na(peak_idx) || peak_idx < 10 || peak_idx > length(y) - 10) {
    return(list(fitted = rep(NA, length(y)), peak_idx = NA))
  }

  left <- y[1:peak_idx]; right <- y[peak_idx:length(y)]
  t_left <- seq_along(left); t_right <- seq_along(right)
  base_left <- min(left, na.rm = TRUE); amp_left <- max(left, na.rm = TRUE) - base_left
  growth_model <- function(t, a1, b1) amp_left * 0.95 / (1 + exp(a1 + b1 * t)) + base_left * 1.05

  popt_growth <- tryCatch(
    nlsLM(left ~ growth_model(t_left, a1, b1), start = list(a1 = 5, b1 = -0.1),
          lower = c(-10, -0.5), upper = c(20, -0.01), control = nls.lm.control(maxiter = 100)),
    error = function(e) NULL
  )

  base_right <- min(right, na.rm = TRUE); amp_right <- max(right, na.rm = TRUE) - base_right
  weibull_model <- function(t, k, lambda, offset = 0.5) {
    amp_right * 0.95 * exp(-((t + offset) / lambda)^k) + base_right * 1.05
  }
  popt_stress <- tryCatch(
    nlsLM(right ~ weibull_model(t_right, k, lambda),
          start = list(k = 3, lambda = length(right) / 3),
          lower = c(0.5, 5), upper = c(15, length(right) * 1.5),
          control = nls.lm.control(maxiter = 200)),
    error = function(e) NULL
  )

  pred <- rep(NA, length(y))
  if (!is.null(popt_growth)) {
    pg <- coef(popt_growth)
    pred[1:peak_idx] <- growth_model(t_left, pg["a1"], pg["b1"])
  } else pred[1:peak_idx] <- left
  if (!is.null(popt_stress)) {
    ps <- coef(popt_stress)
    pred[peak_idx:length(y)] <- weibull_model(t_right, ps["k"], ps["lambda"])
  } else pred[peak_idx:length(y)] <- right

  list(fitted = pmax(pmin(pred, 1.5), -0.2), peak_idx = peak_idx)
}

extract_hplm_year <- function(year_matrix, year_days_idx, daily_dates, ratio_plat = 0.90) {
  n_pixels <- nrow(year_matrix)
  current_doy <- yday(daily_dates[year_days_idx])
  t_seq <- seq_along(year_days_idx)

  annual_max <- apply(year_matrix, 1, max, na.rm = TRUE)
  annual_min <- apply(year_matrix, 1, min, na.rm = TRUE)
  invalid <- is.na(annual_max) | is.na(annual_min) | (annual_max - annual_min) < 0.05
  range_val <- annual_max - annual_min
  range_val[range_val <= 0] <- NA
  year_ratio <- (year_matrix - annual_min) / range_val

  fitted_curve <- matrix(NA, nrow = n_pixels, ncol = length(year_days_idx))
  peak_indices <- rep(NA_integer_, n_pixels)
  for (p in seq_len(n_pixels)) {
    if (invalid[p]) next
    res <- hplm_fit(t_seq, year_ratio[p, ], current_doy)
    if (!all(is.na(res$fitted))) {
      fitted_curve[p, ] <- res$fitted
      peak_indices[p] <- res$peak_idx
    }
  }

  sos_idx <- eos_idx <- ppos_idx <- apos_idx <- rep(NA_integer_, n_pixels)
  for (p in seq_len(n_pixels)) {
    fit <- fitted_curve[p, ]
    if (all(is.na(fit))) next
    peak <- peak_indices[p]
    if (is.na(peak) || peak <= 10 || peak >= length(fit) - 10) next
    left <- fit[1:peak]; right <- fit[peak:length(fit)]
    sos_idx[p] <- which.min(abs(left - (min(left) + max(left)) / 2))
    eos_idx[p] <- peak + which.min(abs(right - (min(right) + max(right)) / 2)) - 1
    fmax <- max(fit, na.rm = TRUE)
    if (fmax > 0) {
      above <- which(fit >= ratio_plat * fmax)
      if (length(above)) {
        ppos_idx[p] <- above[1]
        apos_idx[p] <- above[length(above)]
      }
    }
  }

  SOS_doy <- yday(daily_dates[year_days_idx[sos_idx]])
  EOS_doy <- yday(daily_dates[year_days_idx[eos_idx]])
  PPOS_doy <- yday(daily_dates[year_days_idx[ppos_idx]])
  APOS_doy <- yday(daily_dates[year_days_idx[apos_idx]])

  valid <- !is.na(SOS_doy) & SOS_doy >= 60 & SOS_doy <= 200 &
    !is.na(EOS_doy) & EOS_doy >= 220 & EOS_doy <= 330 &
    SOS_doy < PPOS_doy & PPOS_doy <= APOS_doy & APOS_doy < EOS_doy

  list(
    SOS_DOY = ifelse(valid, SOS_doy, NA),
    PPOS_DOY = ifelse(valid, PPOS_doy, NA),
    APOS_DOY = ifelse(valid, APOS_doy, NA),
    EOS_DOY = ifelse(valid, EOS_doy, NA),
    Greenup_Duration_days = ifelse(valid, PPOS_doy - SOS_doy + 1, NA),
    Plateau_Duration_days = ifelse(valid, APOS_doy - PPOS_doy + 1, NA),
    Senescence_Duration_days = ifelse(valid, EOS_doy - APOS_doy + 1, NA)
  )
}

run_hplm_folder <- function(input_dir, output_dir, years = YEARS_HIST) {
  dat <- load_daily_stack(input_dir)
  daily_dates <- parse_gosif_dates(dat$files)
  ord <- order(daily_dates)
  daily_dates <- daily_dates[ord]
  daily_matrix <- dat$matrix[, ord, drop = FALSE]
  template <- rast(dat$files[1])
  n_pixels <- nrow(daily_matrix)

  init_mat <- function() matrix(NA_real_, nrow = n_pixels, ncol = length(years))
  layers <- list(
    SOS_DOY = init_mat(), PPOS_DOY = init_mat(), APOS_DOY = init_mat(), EOS_DOY = init_mat(),
    Greenup_Duration_days = init_mat(), Plateau_Duration_days = init_mat(),
    Senescence_Duration_days = init_mat()
  )

  for (i in seq_along(years)) {
    idx <- which(year(daily_dates) == years[i])
    if (!length(idx)) next
    yr_out <- extract_hplm_year(daily_matrix[, idx, drop = FALSE], idx, daily_dates)
    for (nm in names(layers)) layers[[nm]][, i] <- yr_out[[nm]]
  }
  save_phenology_layers(template, output_dir, years, layers)
}

run_hplm_batch <- function(parent_input_dir, parent_output_dir, years = YEARS_HIST) {
  subfolders <- list.dirs(parent_input_dir, recursive = FALSE, full.names = TRUE)
  for (folder in subfolders) {
    out <- file.path(parent_output_dir, basename(folder))
    message("HPLM: ", basename(folder))
    run_hplm_folder(folder, out, years)
  }
  invisible(parent_output_dir)
}

if (sys.nframe() == 0L) {
  root <- get_data_root()
  run_hplm_batch(
    parent_input_dir = file.path(root, "gosif_sg"),
    parent_output_dir = file.path(root, "results", "phenology_hplm")
  )
}
