# Dominant PLS-PM path factor per stage (phenology carry-over vs climate pathways).

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  source("R/config.R")
})

stage_path_files <- function(plspm_dir, stage) {
  switch(stage,
    L1 = list(
      pheno_carry = "Path_L3prev_to_L1.tif",
      cli_direct = "Path_CLI_L1_to_L1.tif",
      cli_carry = "Path_CLI_L3prev_to_L1.tif"
    ),
    L2 = list(
      pheno_carry = "Path_L1_to_L2.tif",
      cli_direct = "Path_CLI_L2_to_L2.tif",
      cli_carry = "Path_CLI_L1_to_L2.tif"
    ),
    L3 = list(
      pheno_carry = "Path_L2_to_L3.tif",
      cli_direct = "Path_CLI_L3_to_L3.tif",
      cli_carry = "Path_CLI_L2_to_L3.tif"
    ),
    stop("Unknown stage: ", stage)
  )
}

stronger_env <- function(w) {
  if (is.na(w[1]) || is.na(w[2])) return(NA_real_)
  if (abs(w[1]) >= abs(w[2])) w[1] else w[2]
}

dominant_factor <- function(w) {
  if (is.na(w[1]) || is.na(w[2])) return(NA_real_)
  ap <- abs(w[1]); ae <- abs(w[2])
  if (ap > ae) return(ifelse(w[1] >= 0, 1, -1))
  if (ae > ap) return(ifelse(w[2] >= 0, 2, -2))
  0
}

detailed_dominant <- function(w) {
  if (any(is.na(w))) return(NA_real_)
  av <- abs(w)
  winners <- which(av == max(av))
  if (length(winners) != 1L) return(0)
  val <- w[winners]
  switch(winners,
    ifelse(val >= 0, 1, -1),
    ifelse(val >= 0, 2, -2),
    ifelse(val >= 0, 3, -3)
  )
}

run_stage <- function(plspm_dir, output_dir, stage) {
  files <- stage_path_files(plspm_dir, stage)
  paths <- file.path(plspm_dir, unlist(files))
  if (!all(file.exists(paths))) stop("Missing PLS-PM path rasters for ", stage)

  r_pheno <- rast(paths[1])
  r_direct <- rast(paths[2])
  r_carry <- rast(paths[3])

  env <- app(c(r_direct, r_carry), stronger_env)
  dominant <- app(c(r_pheno, env), dominant_factor)
  detailed <- app(c(r_pheno, r_direct, r_carry), detailed_dominant)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  writeRaster(dominant, file.path(output_dir, paste0("Dominant_Factor_", stage, ".tif")), overwrite = TRUE)
  writeRaster(detailed, file.path(output_dir, paste0("Detailed_Dominant_Factor_", stage, ".tif")), overwrite = TRUE)
  invisible(list(dominant = dominant, detailed = detailed))
}

if (sys.nframe() == 0L) {
  root <- get_data_root()
  plspm_dir <- file.path(root, "results", "plspm", "PlantedForest")
  output_dir <- file.path(plspm_dir, "dominant_path")
  for (st in PHENO_STAGES) run_stage(plspm_dir, output_dir, st)
}
