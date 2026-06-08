# Random forest: yearly climate -> LGS; AI/PR stage decomposition; CMIP6 future projection.

suppressPackageStartupMessages({
  library(terra)
  library(ranger)
  library(dplyr)
  source("R/config.R")
})

CLIM_CFG <- list(
  Tmax  = list(folder = "Air_T_Max", prefix = "Tmax"),
  Tmin  = list(folder = "Air_T_Min", prefix = "Tmin"),
  Pre   = list(folder = "Pre", prefix = "Pre"),
  Srad  = list(folder = "Srad", prefix = "Srad"),
  SoilM = list(folder = "Soil_M", prefix = "Soil_M")
)

forest_clim_subdir <- function(forest) {
  switch(forest, NaturalForest = "yearly_natural_forest", PlantedForest = "yearly_planted_forest",
         stop("Unknown forest: ", forest))
}

resolve_clim_file <- function(dir, prefix, year) {
  for (p in c(paste0(prefix, "_", year, ".tif"), paste0(prefix, year, ".tif"))) {
    f <- file.path(dir, p)
    if (file.exists(f)) return(f)
  }
  hits <- list.files(dir, pattern = paste0(year, ".*\\.tif$"), full.names = TRUE)
  if (length(hits)) hits[1] else NA_character_
}

read_yearly_clim <- function(clim_root, var_key, years, forest) {
  cfg <- CLIM_CFG[[var_key]]
  dir <- file.path(clim_root, cfg$folder, forest_clim_subdir(forest))
  files <- vapply(years, function(y) resolve_clim_file(dir, cfg$prefix, y), "")
  if (any(is.na(files))) stop("Missing climate files in ", dir)
  rast(files)
}

build_pixel_year_df <- function(pheno_root, clim_root, forest, years) {
  lgs <- rast(file.path(pheno_root, "LGS", paste0("LGS_", years, "_", forest, ".tif")))
  mats <- setNames(lapply(names(CLIM_CFG), function(k) read_yearly_clim(clim_root, k, years, forest)), names(CLIM_CFG))
  bind_rows(lapply(seq_along(years), function(i) {
    data.frame(
      Year = years[i], L_total = values(lgs[[i]]),
      Tmax = values(mats$Tmax[[i]]), Tmin = values(mats$Tmin[[i]]),
      Pre = values(mats$Pre[[i]]), Srad = values(mats$Srad[[i]]),
      SoilM = values(mats$SoilM[[i]])
    )
  })) %>% filter(complete.cases(.))
}

calc_metrics <- function(obs, pred) {
  ok <- complete.cases(obs, pred)
  obs <- obs[ok]; pred <- pred[ok]
  if (length(obs) < 2) return(c(R2 = NA, RMSE = NA, MAE = NA, bias = NA, n = length(obs)))
  ss_res <- sum((obs - pred)^2); ss_tot <- sum((obs - mean(obs))^2)
  c(R2 = if (ss_tot > 0) 1 - ss_res / ss_tot else NA,
    RMSE = sqrt(mean((obs - pred)^2)), MAE = mean(abs(obs - pred)),
    bias = mean(pred - obs), n = length(obs))
}

save_importance_plot <- function(imp, out_png, main_title) {
  imp <- sort(imp, decreasing = TRUE)
  png(out_png, width = 2400, height = 1600, res = 300)
  par(mar = c(5, 9, 4, 2))
  barplot(imp, horiz = TRUE, las = 1, col = "#3182bd", border = NA,
          main = main_title, xlab = "Permutation importance")
  dev.off()
}

save_pred_obs_plot <- function(obs, pred, out_png, main_title) {
  ok <- complete.cases(obs, pred)
  obs <- obs[ok]; pred <- pred[ok]
  m <- calc_metrics(obs, pred)
  lim <- range(c(obs, pred), na.rm = TRUE)
  png(out_png, width = 2000, height = 2000, res = 300)
  par(mar = c(5, 5, 4, 2))
  plot(obs, pred, pch = 16, cex = 0.15, col = rgb(0, 0, 0, 0.08),
       xlab = "Observed LGS (days)", ylab = "Predicted LGS (days)",
       main = main_title, xlim = lim, ylim = lim)
  abline(0, 1, col = "red", lwd = 2)
  legend("topleft", sprintf("R2=%.3f RMSE=%.2f n=%s", m["R2"], m["RMSE"], format(m["n"], big.mark = ",")),
         bty = "n")
  dev.off()
}

train_rf_holdout <- function(df, train_years, valid_years) {
  df_tr <- filter(df, Year %in% train_years)
  df_va <- filter(df, Year %in% valid_years)
  rf <- ranger(L_total ~ Tmax + Tmin + Pre + Srad + SoilM,
               data = df_tr, num.trees = 500, mtry = 3, importance = "permutation", seed = 123)
  list(
    model = rf,
    metrics_train = calc_metrics(df_tr$L_total, predict(rf, df_tr)$predictions),
    metrics_valid = calc_metrics(df_va$L_total, predict(rf, df_va)$predictions),
    valid_by_year = df_va %>%
      mutate(pred = predict(rf, df_va)$predictions) %>%
      group_by(Year) %>%
      summarise(
        R2 = calc_metrics(L_total, pred)["R2"],
        RMSE = calc_metrics(L_total, pred)["RMSE"],
        MAE = calc_metrics(L_total, pred)["MAE"],
        bias = calc_metrics(L_total, pred)["bias"],
        n = calc_metrics(L_total, pred)["n"],
        .groups = "drop"
      )
  )
}

decompose_stages <- function(L_total_r, A_r, B_r) {
  A <- resample(A_r, L_total_r, "bilinear")
  B <- resample(B_r, L_total_r, "bilinear")
  L3 <- L_total_r / (1 + A + A * B)
  list(L1 = A * L3, L2 = B * (A * L3 + L3), L3 = L3, L_total = L_total_r)
}

read_future_climate_year <- function(future_clim_dir, year) {
  list(
    Tmax = rast(file.path(future_clim_dir, "tasmax", paste0("tasmax_bias_corrected_", year, ".tif"))),
    Tmin = rast(file.path(future_clim_dir, "tasmin", paste0("tasmin_bias_corrected_", year, ".tif"))),
    Pre  = rast(file.path(future_clim_dir, "pr", paste0("Pre_", year, ".tif"))),
    Srad = rast(file.path(future_clim_dir, "rsds", paste0("rsds_bias_corrected_", year, ".tif"))),
    SoilM = rast(file.path(future_clim_dir, "SoilM", paste0("SoilM_bias_corrected_", year, ".tif")))
  )
}

predict_future_year <- function(rf, future_clim_dir, year, A_r, B_r, skip_existing = TRUE, out_dir = NULL) {
  outs <- c(L_total = file.path(out_dir, paste0("L_total_", year, ".tif")),
            L1 = file.path(out_dir, paste0("L1_", year, ".tif")),
            L2 = file.path(out_dir, paste0("L2_", year, ".tif")),
            L3 = file.path(out_dir, paste0("L3_", year, ".tif")))
  if (skip_existing && all(file.exists(unlist(outs)))) return(invisible(outs))

  clim <- read_future_climate_year(future_clim_dir, year)
  newdata <- as.data.frame(lapply(clim, values))
  pred <- predict(rf, data = newdata)$predictions
  L_total_r <- clim$Tmax
  values(L_total_r) <- pred
  L_total_r <- mask(L_total_r, clim$Tmax)
  stages <- decompose_stages(L_total_r, A_r, B_r)

  if (!is.null(out_dir)) {
    writeRaster(stages$L_total, outs["L_total"], overwrite = TRUE)
    writeRaster(stages$L1, outs["L1"], overwrite = TRUE)
    writeRaster(stages$L2, outs["L2"], overwrite = TRUE)
    writeRaster(stages$L3, outs["L3"], overwrite = TRUE)
  }
  invisible(stages)
}

run_rf_training <- function(root, train_years = 2000:2019, valid_years = 2020:2024) {
  pheno_root <- file.path(root, "phenology")
  clim_root <- file.path(root, "climate")
  ai_pr_root <- file.path(root, "analysis", "ai_pr")
  model_dir <- file.path(root, "results", "rf_models")
  fig_dir <- file.path(model_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  metrics_all <- list()
  for (forest in FOREST_TYPES) {
    A_r <- rast(file.path(ai_pr_root, forest, "AI_25yr_Mean.tif"))
    B_r <- rast(file.path(ai_pr_root, forest, "PR_25yr_Mean.tif"))
    df <- build_pixel_year_df(pheno_root, clim_root, forest, YEARS_HIST)
    ho <- train_rf_holdout(df, train_years, valid_years)

    saveRDS(ho$model, file.path(model_dir, paste0("RF_yearly_", forest, ".rds")))
    saveRDS(A_r, file.path(model_dir, paste0("A_raster_", forest, ".rds")))
    saveRDS(B_r, file.path(model_dir, paste0("B_raster_", forest, ".rds")))
    write.csv(ho$valid_by_year, file.path(model_dir, paste0("holdout_by_year_", forest, ".csv")), row.names = FALSE)
    write.csv(data.frame(feature = names(ho$model$variable.importance),
                         importance = as.numeric(ho$model$variable.importance)),
              file.path(model_dir, paste0("importance_", forest, ".csv")), row.names = FALSE)

    fdir <- file.path(fig_dir, forest)
    dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
    save_importance_plot(ho$model$variable.importance, file.path(fdir, "importance.png"),
                         paste(forest, "RF importance"))
    df_va <- filter(df, Year %in% valid_years)
    save_pred_obs_plot(df_va$L_total, predict(ho$model, df_va)$predictions,
                       file.path(fdir, "pred_vs_obs_holdout.png"), paste(forest, "hold-out validation"))

    metrics_all[[forest]] <- data.frame(
      forest = forest, R2_train = ho$metrics_train["R2"], R2_valid = ho$metrics_valid["R2"],
      RMSE_valid = ho$metrics_valid["RMSE"], MAE_valid = ho$metrics_valid["MAE"]
    )
    message(forest, ": valid R2=", round(ho$metrics_valid["R2"], 3))
  }
  write.csv(bind_rows(metrics_all), file.path(model_dir, "training_metrics_summary.csv"), row.names = FALSE)
  invisible(model_dir)
}

run_rf_future_prediction <- function(root, ssp = "SSP2_4.5", years = YEARS_FUTURE, skip_existing = TRUE) {
  model_dir <- file.path(root, "results", "rf_models")
  output_base <- file.path(root, "cmip6", ssp, "future_phenology")
  future_base <- file.path(root, "cmip6", ssp)

  for (forest in FOREST_TYPES) {
    rf <- readRDS(file.path(model_dir, paste0("RF_yearly_", forest, ".rds")))
    A_r <- readRDS(file.path(model_dir, paste0("A_raster_", forest, ".rds")))
    B_r <- readRDS(file.path(model_dir, paste0("B_raster_", forest, ".rds")))
    forest_out <- file.path(output_base, forest)
    dir.create(forest_out, recursive = TRUE, showWarnings = FALSE)
    future_clim_dir <- file.path(future_base, paste0("yearly_climate_", tolower(forest)))

    for (yr in years) {
      message(forest, " ", ssp, " ", yr)
      predict_future_year(rf, future_clim_dir, yr, A_r, B_r, skip_existing, forest_out)
    }
  }
  invisible(output_base)
}

if (sys.nframe() == 0L) {
  root <- get_data_root()
  mode <- Sys.getenv("RF_MODE", unset = "train")
  if (mode == "train") {
    run_rf_training(root)
  } else if (mode == "predict") {
    run_rf_future_prediction(root, ssp = Sys.getenv("RF_SSP", unset = "SSP2_4.5"))
  } else {
    stop("Set RF_MODE to 'train' or 'predict'")
  }
}
