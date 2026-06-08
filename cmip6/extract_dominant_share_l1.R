# Dominant driver from L1 difference-regression variance shares.
# Codes: 1 = CLI_L1, 2 = carry-over, 3 = CLI_L3prev, 0 = tie, NA = missing

suppressPackageStartupMessages({
  library(terra)
  source("R/config.R")
})

SHARE_FILES <- c(
  cli_l1 = "Share_CLI_L1.tif",
  carry = "Share_carry_L1.tif",
  cli_l3prev = "Share_CLI_L3prev.tif"
)

dominant_from_shares <- function(v) {
  if (length(v) != 3L || any(is.na(v))) return(NA_real_)
  imax <- which(v == max(v))
  if (length(imax) != 1L) return(0)
  as.numeric(imax)
}

summarize_dominant <- function(r, labels) {
  ft <- freq(r)
  total <- sum(ft$count[!is.na(ft$value)], na.rm = TRUE)
  stats <- data.frame(Code = 0:3, Count = 0L, Percentage = NA_real_, Label = labels, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(ft))) {
    val <- ft$value[i]
    if (!is.na(val)) {
      idx <- which(stats$Code == val)
      if (length(idx)) stats$Count[idx] <- ft$count[i]
    }
  }
  if (total > 0) stats$Percentage <- round(stats$Count / total * 100, 2)
  stats
}

run_extract_dominant_share_l1 <- function(input_dir, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(input_dir, SHARE_FILES)
  if (!all(file.exists(paths))) stop("Missing share rasters in: ", input_dir)

  shares <- rast(paths)
  dominant <- app(shares, dominant_from_shares)
  writeRaster(dominant, file.path(output_dir, "Dominant_Share_L1.tif"), overwrite = TRUE)

  labels <- c(
    "Tie",
    "Current-year L1 climate (Share_CLI_L1)",
    "Prior-year phenology (Share_carry_L1)",
    "Prior-year L3 climate (Share_CLI_L3prev)"
  )
  stats <- summarize_dominant(dominant, labels)
  write.csv(stats, file.path(output_dir, "Dominant_Share_L1_Stats.csv"), row.names = FALSE)
  invisible(list(raster = dominant, stats = stats))
}

if (sys.nframe() == 0L) {
  root <- get_data_root()
  for (forest in FOREST_TYPES) {
    input_dir <- file.path(root, "results", "diffreg_l1", "SSP1_2.6", forest)
    output_dir <- file.path(input_dir, "Dominant_Share_L1")
    run_extract_dominant_share_l1(input_dir, output_dir)
  }
}
