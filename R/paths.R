source("R/config.R")

phenology_file <- function(pheno_root, stage, year, forest, future = FALSE) {
  if (future) {
    return(file.path(pheno_root, forest, sprintf("%s_%d.tif", stage, year)))
  }
  if (stage == "Ltotal") {
    return(file.path(pheno_root, "LGS", sprintf("LGS_%d_%s.tif", year, forest)))
  }
  subdir <- switch(stage,
    L1 = "Greenup_Duration",
    L2 = "Plateau_Duration",
    L3 = "Senescence_Duration",
    stop("Unknown stage: ", stage)
  )
  prefix <- switch(stage,
    L1 = "Greenup_Duration_days",
    L2 = "Plateau_Duration_days",
    L3 = "Senescence_Duration_days"
  )
  file.path(pheno_root, subdir, sprintf("%s_%d_%s.tif", prefix, year, forest))
}

stage_climate_file <- function(clim_root, var, stage, year, forest) {
  file.path(clim_root, sprintf("%s_%s_%d_%s.tif", var, stage, year, forest))
}

read_stack <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing files:\n", paste(missing, collapse = "\n"))
  rast(paths)
}
