source("R/config.R")

load_daily_stack <- function(input_dir, pattern = "gosif_daily_.*\\.tif$") {
  files <- sort(list.files(input_dir, pattern = pattern, full.names = TRUE))
  if (!length(files)) stop("No daily GOSIF files in ", input_dir)
  list(files = files, stack = rast(files), matrix = as.matrix(values(rast(files))))
}

parse_gosif_dates <- function(files) {
  parts <- sub("\\.tif$", "", sub("^gosif_daily_", "", basename(files)))
  as.Date(paste(sub("_.*$", "", parts), sub("^.*_", "", parts)), format = "%Y %j")
}

save_phenology_layers <- function(template, output_dir, years, layers) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(layers)) {
    mat <- layers[[nm]]
    for (i in seq_along(years)) {
      r <- template
      values(r) <- mat[, i]
      writeRaster(r, file.path(output_dir, sprintf("%s_%d.tif", nm, years[i])),
                  overwrite = TRUE, gdal = c("COMPRESS=LZW"))
    }
    mean_vec <- rowMeans(mat, na.rm = TRUE)
    mean_vec[is.nan(mean_vec)] <- NA
    r_mean <- template
    values(r_mean) <- mean_vec
    writeRaster(r_mean, file.path(output_dir, sprintf("%s_%d_%d_mean.tif", nm, min(years), max(years))),
                overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  }
  invisible(output_dir)
}
