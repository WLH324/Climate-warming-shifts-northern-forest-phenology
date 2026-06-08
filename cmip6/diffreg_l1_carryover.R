# Pixel-wise first-difference regression for green-up (L1) variability:
#   dL1(t) ~ dCLI_L1(t) + lag(dL3)(t) + lag(dCLI_L3)(t)

suppressPackageStartupMessages({
  library(terra)
  source("R/config.R")
  source("R/paths.R")
})

overwrite <- FALSE
min_complete_obs <- 8L

build_cases <- function(root) {
  list(
    Hist = list(
      years = YEARS_HIST,
      phenology_dir = file.path(root, "phenology"),
      stage_climate_dir = file.path(root, "stage_climate"),
      out_dir = file.path(root, "results", "diffreg_l1", "historical"),
      future = FALSE
    ),
    SSP126 = list(
      years = YEARS_FUTURE,
      phenology_dir = file.path(root, "cmip6", "SSP1_2.6", "future_phenology"),
      stage_climate_dir = file.path(root, "cmip6", "SSP1_2.6", "stage_climate"),
      out_dir = file.path(root, "results", "diffreg_l1", "SSP1_2.6"),
      future = TRUE
    ),
    SSP245 = list(
      years = YEARS_FUTURE,
      phenology_dir = file.path(root, "cmip6", "SSP2_4.5", "future_phenology"),
      stage_climate_dir = file.path(root, "cmip6", "SSP2_4.5", "stage_climate"),
      out_dir = file.path(root, "results", "diffreg_l1", "SSP2_4.5"),
      future = TRUE
    ),
    SSP585 = list(
      years = YEARS_FUTURE,
      phenology_dir = file.path(root, "cmip6", "SSP5_8.5", "future_phenology"),
      stage_climate_dir = file.path(root, "cmip6", "SSP5_8.5", "stage_climate"),
      out_dir = file.path(root, "results", "diffreg_l1", "SSP5_8.5"),
      future = TRUE
    )
  )
}

r2_lm <- function(y, x) {
  ok <- is.finite(y) & apply(x, 1, function(z) all(is.finite(z)))
  if (sum(ok) < min_complete_obs) return(NA_real_)
  yy <- y[ok]; xx <- x[ok, , drop = FALSE]
  if (sd(yy) == 0) return(NA_real_)
  fit <- try(lm.fit(cbind(1, xx), yy), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  tss <- sum((yy - mean(yy))^2)
  if (!is.finite(tss) || tss <= 0) return(NA_real_)
  max(0, 1 - sum(fit$residuals^2) / tss)
}

coef_lm <- function(y, x) {
  ok <- is.finite(y) & apply(x, 1, function(z) all(is.finite(z)))
  p <- ncol(x)
  if (sum(ok) < max(min_complete_obs, p + 2L)) return(rep(NA_real_, p))
  yy <- y[ok]; xx <- x[ok, , drop = FALSE]
  if (sd(yy) == 0) return(rep(NA_real_, p))
  fit <- try(lm.fit(cbind(1, xx), yy), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  as.numeric(fit$coefficients[-1])
}

calc_pixel_l1 <- function(v, n_years, clim_names) {
  idx <- 0L
  pull <- function(n) { idx <<- idx + n; v[(idx - n + 1):idx] }

  L1 <- pull(n_years)
  L3 <- pull(n_years)
  cli_l1 <- lapply(seq_along(clim_names), function(i) pull(n_years))
  cli_l3 <- lapply(seq_along(clim_names), function(i) pull(n_years))

  dL1 <- diff(L1)
  dL3 <- diff(L3)
  d_cli_l1 <- lapply(cli_l1, diff)
  d_cli_l3 <- lapply(cli_l3, diff)

  y <- dL1[-1]
  x_cli_l1 <- do.call(cbind, lapply(d_cli_l1, function(z) z[-1]))
  x_carry <- matrix(dL3[-length(dL3)], ncol = 1)
  x_cli_l3prev <- do.call(cbind, lapply(d_cli_l3, function(z) z[-length(z)]))
  x_full <- cbind(x_cli_l1, x_carry, x_cli_l3prev)

  n_ok <- sum(is.finite(y) & apply(x_full, 1, function(z) all(is.finite(z))))
  r2_full <- r2_lm(y, x_full)
  if (!is.finite(r2_full)) {
    return(c(
      R2_L1 = NA, DeltaR2_CLI_L1 = NA, DeltaR2_carry_L1 = NA, DeltaR2_CLI_L3prev = NA,
      Share_CLI_L1 = NA, Share_carry_L1 = NA, Share_CLI_L3prev = NA,
      setNames(rep(NA, length(clim_names) * 2 + 1),
               c(paste0("Coef_", clim_names, "_L1"), "Coef_L3prev_to_L1",
                 paste0("Coef_", clim_names, "_L3prev"))),
      N_years_L1 = n_ok
    ))
  }

  r2_no_l1 <- r2_lm(y, cbind(x_carry, x_cli_l3prev))
  r2_no_carry <- r2_lm(y, cbind(x_cli_l1, x_cli_l3prev))
  r2_no_l3prev <- r2_lm(y, cbind(x_cli_l1, x_carry))
  d <- c(
    l1 = max(0, r2_full - r2_no_l1),
    carry = max(0, r2_full - r2_no_carry),
    l3prev = max(0, r2_full - r2_no_l3prev)
  )
  share <- function(x) if (r2_full > 0) x / r2_full else NA_real_
  coefs <- coef_lm(y, x_full)
  names(coefs) <- c(
    paste0("Coef_", clim_names, "_L1"), "Coef_L3prev_to_L1",
    paste0("Coef_", clim_names, "_L3prev")
  )
  c(
    R2_L1 = r2_full,
    DeltaR2_CLI_L1 = d["l1"], DeltaR2_carry_L1 = d["carry"], DeltaR2_CLI_L3prev = d["l3prev"],
    Share_CLI_L1 = share(d["l1"]), Share_carry_L1 = share(d["carry"]), Share_CLI_L3prev = share(d["l3prev"]),
    coefs, N_years_L1 = n_ok
  )
}

run_case <- function(cfg, forest, clim_names = names(CLIM_VARS)) {
  out_dir <- file.path(cfg$out_dir, forest)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(file.path(out_dir, "R2_L1.tif")) && !overwrite) return(invisible(NULL))

  years <- cfg$years
  n <- length(years)
  L1 <- read_stack(vapply(years, function(y) phenology_file(cfg$phenology_dir, "L1", y, forest, cfg$future), ""))
  L3 <- read_stack(vapply(years, function(y) phenology_file(cfg$phenology_dir, "L3", y, forest, cfg$future), ""))
  cli_l1 <- lapply(CLIM_VARS, function(vn) {
    read_stack(vapply(years, function(y) stage_climate_file(cfg$stage_climate_dir, vn, "L1", y, forest), ""))
  })
  cli_l3 <- lapply(CLIM_VARS, function(vn) {
    read_stack(vapply(years, function(y) stage_climate_file(cfg$stage_climate_dir, vn, "L3", y, forest), ""))
  })

  s <- c(L1, L3, rast(cli_l1), rast(cli_l3))
  result <- app(s, function(v) calc_pixel_l1(v, n, clim_names))
  for (nm in names(result)) {
    writeRaster(result[[nm]], file.path(out_dir, paste0(nm, ".tif")), overwrite = overwrite)
  }
  invisible(out_dir)
}

if (sys.nframe() == 0L) {
  root <- get_data_root()
  for (cfg in build_cases(root)) {
    for (forest in FOREST_TYPES) run_case(cfg, forest)
  }
}
