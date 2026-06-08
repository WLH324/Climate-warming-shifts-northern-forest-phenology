#' Project configuration. Set GOSIF_DATA_ROOT before running scripts.
#' @examples
#'   Sys.setenv(GOSIF_DATA_ROOT = "F:/GO_SIF")
get_data_root <- function() {
  root <- Sys.getenv("GOSIF_DATA_ROOT", unset = "data")
  normalizePath(root, winslash = "/", mustWork = FALSE)
}

FOREST_TYPES <- c("NaturalForest", "PlantedForest")
CLIM_VARS <- c(
  Tmax  = "Air_T_Max",
  Tmin  = "Air_T_Min",
  Pre   = "Pre",
  SoilM = "Soil_M",
  Srad  = "Srad"
)
PHENO_STAGES <- c("L1", "L2", "L3")
YEARS_HIST <- 2000:2024
YEARS_FUTURE <- 2025:2100
