# Classify temporal stability of stage time allocation (constant / linear / nonlinear).

suppressPackageStartupMessages({
  library(terra)
  library(segmented)
  library(lmtest)
  source("R/config.R")
})

classify_pixel <- function(ts, years) {
  valid <- !is.na(ts)
  if (sum(valid) < 15) return(c(class = NA_integer_, pvalue = NA_real_))
  y <- ts[valid]; x <- years[valid]
  lm_mod <- lm(y ~ x)
  p_lin <- pf(summary(lm_mod)$fstatistic[1], summary(lm_mod)$fstatistic[2],
              summary(lm_mod)$fstatistic[3], lower.tail = FALSE)
  if (p_lin > 0.05) return(c(class = 1L, pvalue = p_lin))
  seg_mod <- try(segmented(lm_mod, seg.Z = ~x), silent = TRUE)
  if (inherits(seg_mod, "try-error")) return(c(class = 2L, pvalue = p_lin))
  p_seg <- davies.test(lm_mod, seg.Z = ~x)$p.value
  c(class = ifelse(p_seg <= 0.05, 3L, 2L), pvalue = p_lin)
}

run_classify_allocation <- function(greenup_dir, senescence_dir, plateau_dir, forest, years = YEARS_HIST) {
  vgu <- rast(file.path(greenup_dir, paste0("Greenup_Duration_days_", years, "_", forest, ".tif")))
  vss <- rast(file.path(senescence_dir, paste0("Senescence_Duration_days_", years, "_", forest, ".tif")))
  plt <- rast(file.path(plateau_dir, paste0("Plateau_Duration_days_", years, "_", forest, ".tif")))
  rp <- plt / (vgu + vss)
  rp[vgu + vss == 0] <- NA

  out <- app(rp, function(v) classify_pixel(v, years))
  names(out) <- c("Allocation_class", "Linear_pvalue")
  out
}

if (sys.nframe() == 0L) {
  root <- get_data_root()
  pheno <- file.path(root, "phenology")
  out <- run_classify_allocation(
    file.path(pheno, "Greenup_Duration"),
    file.path(pheno, "Senescence_Duration"),
    file.path(pheno, "Plateau_Duration"),
    "NaturalForest"
  )
  writeRaster(out, file.path(root, "results", "allocation_class.tif"), overwrite = TRUE)
}
