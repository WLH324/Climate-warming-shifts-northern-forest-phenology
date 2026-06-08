# Pixel-wise three-stage PLS-PM with full carry-over paths and bootstrap significance.

suppressPackageStartupMessages({
  library(terra)
  library(plspm)
  source("R/config.R")
})

extract_cli_weights <- function(pls, block_name, var_order) {
  w <- setNames(rep(NA_real_, length(var_order)), var_order)
  src <- if (!is.null(pls$weights)) pls$weights else pls$outer_model
  if (is.null(src) || nrow(src) == 0) return(w)
  colnames(src) <- tolower(colnames(src))
  rows <- src[tolower(src$block) == tolower(block_name), , drop = FALSE]
  wcol <- if ("weight" %in% colnames(rows)) "weight" else if ("loading" %in% colnames(rows)) "loading" else NULL
  ncol_nm <- if ("name" %in% colnames(rows)) "name" else if ("manifest" %in% colnames(rows)) "manifest" else NULL
  if (is.null(wcol) || is.null(ncol_nm)) return(w)
  for (i in seq_along(var_order)) {
    hit <- rows[tolower(rows[[ncol_nm]]) == tolower(var_order[i]), , drop = FALSE]
    if (nrow(hit)) w[i] <- hit[[wcol]][1]
  }
  w
}

boot_path_sig <- function(pls, pred, resp) {
  if (is.null(pls$boot) || is.null(pls$boot$paths)) return(NA_integer_)
  bp <- as.data.frame(pls$boot$paths)
  rn <- rownames(bp)
  targets <- c(paste0(pred, " -> ", resp), paste0(pred, "->", resp))
  hit <- which(rn %in% targets)
  if (!length(hit)) {
    pat <- paste0("^", pred, "\\s*(->|~)\\s*", resp, "$")
    hit <- grep(pat, rn)
  }
  if (length(hit) != 1L) return(NA_integer_)
  lo_col <- grep("0\\.025|\\.025", colnames(bp), value = TRUE)[1]
  hi_col <- grep("0\\.975|\\.975", colnames(bp), value = TRUE)[1]
  if (is.na(lo_col) || is.na(hi_col)) return(NA_integer_)
  lo <- suppressWarnings(as.numeric(bp[hit, lo_col]))
  hi <- suppressWarnings(as.numeric(bp[hit, hi_col]))
  as.integer(is.finite(lo) && is.finite(hi) && lo * hi > 0)
}

clean_path <- function(x) if (is.na(x) || abs(x) > 1) NA_real_ else x

read_plspm_inputs <- function(clim_root, pheno_root, forest, years = YEARS_HIST) {
  read_stage <- function(var, stage) {
    rast(file.path(clim_root, paste0(var, "_", stage, "_", years, "_", forest, ".tif")))
  }
  L1_dur <- rast(file.path(pheno_root, "Greenup_Duration",
                           paste0("Greenup_Duration_days_", years, "_", forest, ".tif")))
  L2_dur <- rast(file.path(pheno_root, "Plateau_Duration",
                           paste0("Plateau_Duration_days_", years, "_", forest, ".tif")))
  L3_dur <- rast(file.path(pheno_root, "Senescence_Duration",
                           paste0("Senescence_Duration_days_", years, "_", forest, ".tif")))
  lag_stack <- function(r) c(r[[1]] * NA, r[[1:(nlyr(r) - 1)]])

  c(
    read_stage("Air_T_Max", "L1"), read_stage("Air_T_Min", "L1"), read_stage("Pre", "L1"),
    read_stage("Soil_M", "L1"), read_stage("Srad", "L1"), L1_dur,
    read_stage("Air_T_Max", "L2"), read_stage("Air_T_Min", "L2"), read_stage("Pre", "L2"),
    read_stage("Soil_M", "L2"), read_stage("Srad", "L2"), L2_dur,
    read_stage("Air_T_Max", "L3"), read_stage("Air_T_Min", "L3"), read_stage("Pre", "L3"),
    read_stage("Soil_M", "L3"), read_stage("Srad", "L3"), L3_dur,
    lag_stack(read_stage("Air_T_Max", "L3")), lag_stack(read_stage("Air_T_Min", "L3")),
    lag_stack(read_stage("Pre", "L3")), lag_stack(read_stage("Soil_M", "L3")),
    lag_stack(read_stage("Srad", "L3")), lag_stack(L3_dur)
  )
}

build_plspm_model <- function() {
  latent <- c("CLI_L1", "CLI_L2", "CLI_L3", "CLI_L3_prev", "L3_prev", "L1", "L2", "L3")
  inner <- matrix(0, 8, 8, dimnames = list(latent, latent))
  inner["L1", c("CLI_L1", "L3_prev", "CLI_L3_prev")] <- 1
  inner["L2", c("L1", "CLI_L1", "CLI_L2")] <- 1
  inner["L3", c("L2", "CLI_L2", "CLI_L3")] <- 1
  outer <- list(
    CLI_L1 = c("L1_Air_T_Max", "L1_Air_T_Min", "L1_Pre", "L1_Soil_M", "L1_Srad"),
    CLI_L2 = c("L2_Air_T_Max", "L2_Air_T_Min", "L2_Pre", "L2_Soil_M", "L2_Srad"),
    CLI_L3 = c("L3_Air_T_Max", "L3_Air_T_Min", "L3_Pre", "L3_Soil_M", "L3_Srad"),
    CLI_L3_prev = c("L3_prev_Air_T_Max", "L3_prev_Air_T_Min", "L3_prev_Pre",
                    "L3_prev_Soil_M", "L3_prev_Srad"),
    L3_prev = "L3_prev", L1 = "L1", L2 = "L2", L3 = "L3"
  )
  list(inner = inner, outer = outer, modes = c("B", "B", "B", "B", "A", "A", "A", "A"))
}

make_pixel_plspm_fun <- function(model, boot = TRUE, br = 200L) {
  function(v) {
    n_var <- 24L; n_out <- 44L
    result <- rep(NA_real_, n_out)
    if (sum(!is.na(v)) < 480L || length(v) %% n_var != 0L) return(result)
    n_y <- as.integer(length(v) %/% n_var)
    if (n_y < 20L) return(result)

    mat <- matrix(v, nrow = n_y, ncol = n_var)
    df <- as.data.frame(mat)
    colnames(df) <- c(
      "L1_Air_T_Max", "L1_Air_T_Min", "L1_Pre", "L1_Soil_M", "L1_Srad", "L1",
      "L2_Air_T_Max", "L2_Air_T_Min", "L2_Pre", "L2_Soil_M", "L2_Srad", "L2",
      "L3_Air_T_Max", "L3_Air_T_Min", "L3_Pre", "L3_Soil_M", "L3_Srad", "L3",
      "L3_prev_Air_T_Max", "L3_prev_Air_T_Min", "L3_prev_Pre",
      "L3_prev_Soil_M", "L3_prev_Srad", "L3_prev"
    )
    df <- df[complete.cases(df), , drop = FALSE]
    if (nrow(df) < 20L || any(apply(df, 2L, sd, na.rm = TRUE) == 0, na.rm = TRUE)) return(result)

    pls <- tryCatch(
      plspm(df, model$inner, model$outer, model$modes,
            scaled = TRUE, boot.val = boot, br = br, maxiter = 1000, tol = 1e-6),
      error = function(e) NULL
    )
    if (is.null(pls)) return(result)

    pc <- pls$path_coefs
    inner_sum <- as.data.frame(pls$inner_summary)
    result[1:9] <- c(
      clean_path(pc["L1", "CLI_L1"]), clean_path(pc["L1", "L3_prev"]), clean_path(pc["L1", "CLI_L3_prev"]),
      clean_path(pc["L2", "L1"]), clean_path(pc["L2", "CLI_L1"]), clean_path(pc["L2", "CLI_L2"]),
      clean_path(pc["L3", "L2"]), clean_path(pc["L3", "CLI_L2"]), clean_path(pc["L3", "CLI_L3"])
    )
    r2 <- function(nm) if (nm %in% rownames(inner_sum)) inner_sum[nm, "R2"] else NA_real_
    w1 <- extract_cli_weights(pls, "CLI_L1", model$outer$CLI_L1)
    w2 <- extract_cli_weights(pls, "CLI_L2", model$outer$CLI_L2)
    w3 <- extract_cli_weights(pls, "CLI_L3", model$outer$CLI_L3)
    w3p <- extract_cli_weights(pls, "CLI_L3_prev", model$outer$CLI_L3_prev)
    result[10:35] <- c(pls$gof, r2("L1"), r2("L2"), r2("L3"), w1, w2, w3, w3p,
                       mean(pls$outer_model$communality, na.rm = TRUE),
                       mean(pls$outer_model$redundancy, na.rm = TRUE))
    result[36:44] <- c(
      boot_path_sig(pls, "CLI_L1", "L1"), boot_path_sig(pls, "L3_prev", "L1"), boot_path_sig(pls, "CLI_L3_prev", "L1"),
      boot_path_sig(pls, "L1", "L2"), boot_path_sig(pls, "CLI_L1", "L2"), boot_path_sig(pls, "CLI_L2", "L2"),
      boot_path_sig(pls, "L2", "L3"), boot_path_sig(pls, "CLI_L2", "L3"), boot_path_sig(pls, "CLI_L3", "L3")
    )
    result
  }
}

OUTPUT_NAMES <- c(
  "Path_CLI_L1_to_L1", "Path_L3prev_to_L1", "Path_CLI_L3prev_to_L1",
  "Path_L1_to_L2", "Path_CLI_L1_to_L2", "Path_CLI_L2_to_L2",
  "Path_L2_to_L3", "Path_CLI_L2_to_L3", "Path_CLI_L3_to_L3",
  "GOF", "R2_L1", "R2_L2", "R2_L3",
  paste0(rep(c("Weight_Tmax", "Weight_Tmin", "Weight_Pre", "Weight_SoilM", "Weight_Srad"), 4),
         rep(c("_L1", "_L2", "_L3", "_L3prev"), each = 5)),
  "Avg_Communality", "Avg_Redundancy",
  paste0("Sig_", c("CLI_L1_to_L1", "L3prev_to_L1", "CLI_L3prev_to_L1",
                   "L1_to_L2", "CLI_L1_to_L2", "CLI_L2_to_L2",
                   "L2_to_L3", "CLI_L2_to_L3", "CLI_L3_to_L3"))
)

run_plspm_pixelwise <- function(clim_root, pheno_root, output_dir, forest,
                                extent = NULL, boot = TRUE, br = 200L, cores = 1L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  model <- build_plspm_model()
  stack <- read_plspm_inputs(clim_root, pheno_root, forest)
  if (!is.null(extent)) stack <- crop(stack, extent)

  pls_out <- app(stack, make_pixel_plspm_fun(model, boot, br), cores = cores)
  names(pls_out) <- OUTPUT_NAMES
  for (i in seq_len(nlyr(pls_out))) {
    writeRaster(pls_out[[i]], file.path(output_dir, paste0(names(pls_out)[i], ".tif")), overwrite = TRUE)
  }
  invisible(pls_out)
}

if (sys.nframe() == 0L) {
  root <- get_data_root()
  run_plspm_pixelwise(
    clim_root = file.path(root, "stage_climate"),
    pheno_root = file.path(root, "phenology"),
    output_dir = file.path(root, "results", "plspm", "NaturalForest"),
    forest = "NaturalForest",
    extent = NULL,
    boot = TRUE,
    br = 200L
  )
}
