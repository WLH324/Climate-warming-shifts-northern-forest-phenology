# Pixel-wise first-difference regression for growing-season length:
#   dL_total(t) ~ dCLI_L1(t) + dCLI_L2(t) + dCLI_L3(t)

suppressPackageStartupMessages({
  library(terra)
  source("R/config.R")
  source("R/paths.R")
})

overwrite <- FALSE
min_complete_obs <- 8L
stages <- PHENO_STAGES

build_cases <- function(root) {
  list(
    Hist = list(
      years = YEARS_HIST,
      phenology_dir = file.path(root, "phenology"),
      stage_climate_dir = file.path(root, "stage_climate"),
      out_dir = file.path(root, "results", "diffreg_ltotal", "historical"),
      future = FALSE
    ),
    SSP126 = list(
      years = YEARS_FUTURE,
      phenology_dir = file.path(root, "cmip6", "SSP1_2.6", "future_phenology"),
      stage_climate_dir = file.path(root, "cmip6", "SSP1_2.6", "stage_climate"),
      out_dir = file.path(root, "results", "diffreg_ltotal", "SSP1_2.6"),
      future = TRUE
    ),
    SSP245 = list(
      years = YEARS_FUTURE,
      phenology_dir = file.path(root, "cmip6", "SSP2_4.5", "future_phenology"),
      stage_climate_dir = file.path(root, "cmip6", "SSP2_4.5", "stage_climate"),
      out_dir = file.path(root, "results", "diffreg_ltotal", "SSP2_4.5"),
      future = TRUE
    ),
    SSP585 = list(
      years = YEARS_FUTURE,
      phenology_dir = file.path(root, "cmip6", "SSP5_8.5", "future_phenology"),
      stage_climate_dir = file.path(root, "cmip6", "SSP5_8.5", "stage_climate"),
      out_dir = file.path(root, "results", "diffreg_ltotal", "SSP5_8.5"),
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

calc_pixel_ltotal <- function(v, n_years, clim_names, stage_names) {
  idx <- 0L
  pull <- function(n) { idx <<- idx + n; v[(idx - n + 1):idx] }

  ltotal <- pull(n_years)
  cli <- lapply(stage_names, function(st) lapply(clim_names, function(i) pull(n_years)))

  y <- diff(ltotal)
  x_blocks <- lapply(cli, function(st_layers) do.call(cbind, lapply(st_layers, diff)))
  x_full <- do.call(cbind, x_blocks)
  n_ok <- sum(is.finite(y) & apply(x_full, 1, function(z) all(is.finite(z))))
  r2_full <- r2_lm(y, x_full)

  coef_names <- unlist(lapply(stage_names, function(st) paste0("Coef_", clim_names, "_CLI_", st)))
  if (!is.finite(r2_full)) {
    return(c(
      R2_Ltotal = NA,
      DeltaR2_CLI_L1 = NA, DeltaR2_CLI_L2 = NA, DeltaR2_CLI_L3 = NA,
      Share_CLI_L1_on_Ltotal = NA, Share_CLI_L2_on_Ltotal = NA, Share_CLI_L3_on_Ltotal = NA,
      setNames(rep(NA, length(coef_names)), coef_names),
      Dominant_CLI_Ltotal = NA, N_years_Ltotal = n_ok
    ))
  }

  r2_without <- lapply(stage_names, function(st) {
    keep <- setdiff(stage_names, st)
    r2_lm(y, do.call(cbind, x_blocks[keep]))
  })
  delta <- sapply(stage_names, function(st) max(0, r2_full - r2_without[[st]]))
  names(delta) <- stage_names
  share <- function(d) if (r2_full > 0) d / r2_full else NA_real_
  coefs <- coef_lm(y, x_full)
  names(coefs) <- coef_names
  dominant <- if (max(delta, na.rm = TRUE) > 0) which.max(delta) else NA_real_

  c(
    R2_Ltotal = r2_full,
    DeltaR2_CLI_L1 = delta["L1"], DeltaR2_CLI_L2 = delta["L2"], DeltaR2_CLI_L3 = delta["L3"],
    Share_CLI_L1_on_Ltotal = share(delta["L1"]),
    Share_CLI_L2_on_Ltotal = share(delta["L2"]),
    Share_CLI_L3_on_Ltotal = share(delta["L3"]),
    coefs, Dominant_CLI_Ltotal = dominant, N_years_Ltotal = n_ok
  )
}

run_case <- function(cfg, forest, clim_names = names(CLIM_VARS)) {
  out_dir <- file.path(cfg$out_dir, forest)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(file.path(out_dir, "R2_Ltotal.tif")) && !overwrite) return(invisible(NULL))

  years <- cfg$years
  n <- length(years)
  ltotal <- if (cfg$future) {
    L1 <- read_stack(vapply(years, function(y) phenology_file(cfg$phenology_dir, "L1", y, forest, TRUE), ""))
    L2 <- read_stack(vapply(years, function(y) phenology_file(cfg$phenology_dir, "L2", y, forest, TRUE), ""))
    L3 <- read_stack(vapply(years, function(y) phenology_file(cfg$phenology_dir, "L3", y, forest, TRUE), ""))
    L1 + L2 + L3
  } else {
    read_stack(vapply(years, function(y) phenology_file(cfg$phenology_dir, "Ltotal", y, forest, FALSE), ""))
  }

  cli_stacks <- list()
  for (st in stages) {
    for (vn in CLIM_VARS) {
      cli_stacks[[paste(st, vn, sep = "_")]] <- read_stack(
        vapply(years, function(y) stage_climate_file(cfg$stage_climate_dir, vn, st, y, forest), "")
      )
    }
  }

  s <- c(ltotal, rast(cli_stacks))
  result <- app(s, function(v) calc_pixel_ltotal(v, n, clim_names, stages))
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
