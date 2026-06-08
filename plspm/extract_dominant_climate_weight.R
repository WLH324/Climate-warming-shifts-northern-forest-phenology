# Dominant climate variable from PLS-PM outer weights (signed codes +/-1..5).

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  source("R/config.R")
})

WEIGHT_PREFIXES <- c("Weight_Tmax", "Weight_Tmin", "Weight_Pre", "Weight_SoilM", "Weight_Srad")
FACTOR_ID <- setNames(1:5, WEIGHT_PREFIXES)
FACTOR_LABELS <- c("Tmax", "Tmin", "Pre", "SoilM", "Srad")

stronger_weight <- function(w) {
  if (is.na(w[1]) && is.na(w[2])) return(NA_real_)
  if (is.na(w[1])) return(w[2])
  if (is.na(w[2])) return(w[1])
  if (abs(w[1]) >= abs(w[2])) w[1] else w[2]
}

dominant_climate <- function(w) {
  if (all(is.na(w))) return(NA_real_)
  aw <- abs(w)
  if (all(aw == 0 | is.na(aw))) return(0)
  idx <- which.max(aw)
  sign(w[idx]) * FACTOR_ID[idx]
}

run_stage <- function(plspm_dir, output_dir, stage) {
  if (stage == "L1") {
    curr <- file.path(plspm_dir, paste0(WEIGHT_PREFIXES, "_L1.tif"))
    prev <- file.path(plspm_dir, paste0(WEIGHT_PREFIXES, "_L3prev.tif"))
    if (!all(file.exists(c(curr, prev)))) stop("Missing L1/L3prev weight rasters")
    curr_stack <- rast(curr)
    prev_stack <- rast(prev)
    layers <- lapply(seq_along(WEIGHT_PREFIXES), function(i) {
      app(c(curr_stack[[i]], prev_stack[[i]]), stronger_weight)
    })
    weight_stack <- rast(layers)
  } else {
    files <- file.path(plspm_dir, paste0(WEIGHT_PREFIXES, "_", stage, ".tif"))
    if (!all(file.exists(files))) stop("Missing weight rasters for ", stage)
    weight_stack <- rast(files)
  }

  dominant <- app(weight_stack, dominant_climate)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out_tif <- file.path(output_dir, paste0("Dominant_Climate_", stage, ".tif"))
  writeRaster(dominant, out_tif, overwrite = TRUE)

  ft <- freq(dominant)
  total <- sum(ft$count[!is.na(ft$value)], na.rm = TRUE)
  stats <- data.frame(Factor = -5:5, Count = 0L, Percentage = 0, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(ft))) {
    val <- ft$value[i]
    if (!is.na(val)) {
      idx <- which(stats$Factor == val)
      if (length(idx)) stats$Count[idx] <- ft$count[i]
    }
  }
  if (total > 0) stats$Percentage <- round(stats$Count / total * 100, 2)
  stats$Label <- ifelse(stats$Factor == 0, "None", FACTOR_LABELS[abs(stats$Factor)])
  stats$Direction <- ifelse(stats$Factor == 0, "None",
                            ifelse(stats$Factor > 0, "Positive", "Negative"))
  write.csv(stats, file.path(output_dir, paste0("Dominant_Climate_Stats_", stage, ".csv")), row.names = FALSE)
  invisible(dominant)
}

if (sys.nframe() == 0L) {
  root <- get_data_root()
  plspm_dir <- file.path(root, "results", "plspm", "NaturalForest")
  output_dir <- file.path(plspm_dir, "dominant_climate")
  for (st in PHENO_STAGES) run_stage(plspm_dir, output_dir, st)
}
