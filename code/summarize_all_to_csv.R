#!/usr/bin/env Rscript
# summarize_all_to_csv.R
# Aggregate per-replicate RDS outputs into a flat CSV table per scenario.
#
# Usage:
#   Rscript summarize_all_to_csv.R [results_dir] [output_dir] [scenarios_csv]
#
# Defaults:
#   results_dir   = results
#   output_dir    = csv_output
#   scenarios_csv = SCENARIOS.csv  (in the working directory or one level up)
#
# Inputs:
#   <results_dir>/<scenario_id>/rep_*.rds   (one R list per replicate)
#   <scenarios_csv>                          (true values per scenario)
#
# Outputs (one row per row described):
#   sim_meta.csv      one row per scenario   (n, beta1, censoring, t0, true values)
#   sim_params.csv    one row per (scenario, parameter)  (bias, ESE, ASE, SER, ECP)
#   sim_xc.csv        one row per (scenario, method)     (Xc estimates, SE, ECP)
#   sim_cindex.csv    one row per (scenario, method)     (C-index train and test)
#   sim_survival.csv  one row per (scenario, method)     (S(t0) at reference covariate)
#   sim_cr.csv        one row per (scenario, method)     (critical region bounds)
#   sim_boot_diag.csv one row per scenario               (bootstrap diagnostics)
#
# The summarizer caps each scenario at the first 1000 successful replicates by default.
# Modify N_TARGET below for a different cap.

N_TARGET <- 1000

args <- commandArgs(trailingOnly = TRUE)
basedir       <- if (length(args) >= 1) args[1] else "results"
outdir        <- if (length(args) >= 2) args[2] else "csv_output"
scenarios_csv <- if (length(args) >= 3) args[3] else NULL

if (is.null(scenarios_csv)) {
  for (cand in c("SCENARIOS.csv", file.path("..", "SCENARIOS.csv"))) {
    if (file.exists(cand)) { scenarios_csv <- cand; break }
  }
}
if (is.null(scenarios_csv) || !file.exists(scenarios_csv)) {
  stop("SCENARIOS.csv not found. Pass its path as the third argument.")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Load scenarios and build a lookup for true values
scen_tab <- read.csv(scenarios_csv, stringsAsFactors = FALSE)
scenarios <- scen_tab$scenario_id

true_vals <- function(sid) {
  row <- scen_tab[scen_tab$scenario_id == sid, ]
  if (nrow(row) != 1) stop(sprintf("Scenario %s not found in %s", sid, scenarios_csv))
  list(s_t0 = row$s_t0_true, crlb = row$crlb_true, crub = row$crub_true)
}

pnames <- c("beta0", "log_b1", "alpha1", "alpha2")

meta_rows    <- list()
param_rows   <- list()
xc_rows      <- list()
cindex_rows  <- list()
surv_rows    <- list()
cr_rows      <- list()
boot_rows    <- list()

# Helper: safely extract a competitor field
sf <- function(res_list, method, field) {
  sapply(res_list, function(r) {
    v <- r$competitors[[method]][[field]]
    if (is.null(v)) NA else v
  })
}

for (sid in scenarios) {
  files <- list.files(file.path(basedir, sid), pattern = "[.]rds$", full.names = TRUE)
  if (length(files) == 0) {
    cat(sprintf("%s: no files found, skipping\n", sid))
    next
  }
  files <- files[1:min(N_TARGET, length(files))]
  cat(sprintf("Reading %s: %d files ...\n", sid, length(files)))
  res  <- lapply(files, readRDS)
  meta <- res[[1]]
  nrep <- length(res)

  # Provenance checks
  prov_ok <- TRUE
  for (pf in c("ei", "G_type", "t0_S", "t0_CR")) {
    vals <- sapply(res, function(r) as.character(r[[pf]]))
    if (length(unique(vals)) > 1) {
      cat(sprintf("  *** %s PROVENANCE MISMATCH on %s: %s ***\n",
                  sid, pf, paste(unique(vals), collapse = ", ")))
      prov_ok <- FALSE
    }
  }
  b1_vals <- sapply(res, function(r) if (!is.null(r$beta1)) r$beta1 else 1)
  if (length(unique(b1_vals)) > 1) {
    cat(sprintf("  *** %s PROVENANCE MISMATCH on beta1 ***\n", sid))
    prov_ok <- FALSE
  }
  if (!prov_ok) {
    cat(sprintf("  *** SKIPPING %s ***\n", sid))
    next
  }

  b1 <- unique(b1_vals)
  tv <- true_vals(sid)

  # ---- 1. Metadata ----
  meta_rows[[sid]] <- data.frame(
    scenario   = sid,
    n          = meta$n,
    G_type     = meta$G_type,
    ei         = meta$ei,
    beta1      = b1,
    cens_rate  = round(mean(sapply(res, function(r) r$censor_rate)), 3),
    cp_true    = meta$cp_true,
    s_t0_true  = tv$s_t0,
    crlb_true  = tv$crlb,
    crub_true  = tv$crub,
    t0_S       = meta$t0_S,
    t0_CR      = meta$t0_CR,
    nrep       = nrep,
    stringsAsFactors = FALSE
  )

  # ---- 2. MCE Parameter Estimation ----
  par_mat  <- t(sapply(res, function(r) r$est_par))
  bias_mat <- t(sapply(res, function(r) r$bias))
  bse_mat  <- t(sapply(res, function(r) r$boot_se))
  cov_mat  <- t(sapply(res, function(r) r$cov_par))

  for (j in 1:4) {
    ese <- sd(par_mat[, j], na.rm = TRUE)
    ase <- mean(bse_mat[, j], na.rm = TRUE)
    param_rows[[paste0(sid, "_", pnames[j])]] <- data.frame(
      scenario = sid,
      param    = pnames[j],
      bias     = round(mean(bias_mat[, j], na.rm = TRUE), 6),
      ESE      = round(ese, 6),
      ASE      = round(ase, 6),
      SE_ratio = round(ase / ese, 4),
      ECP      = round(mean(cov_mat[, j], na.rm = TRUE), 4),
      stringsAsFactors = FALSE
    )
  }

  # ---- 3. Critical Point Xc ----
  cp_vec    <- sapply(res, function(r) r$est_cpoint)
  cp_bias   <- sapply(res, function(r) r$bias_cp)
  cp_se_d   <- sapply(res, function(r) r$boot_cp_se)
  cp_se_dir <- sapply(res, function(r) r$boot_cp_se_direct)
  cp_cov_d  <- sapply(res, function(r) r$cov_cp)
  ese_cp    <- sd(cp_vec, na.rm = TRUE)

  ci_dir  <- cbind(cp_vec - qnorm(0.975) * cp_se_dir,
                   cp_vec + qnorm(0.975) * cp_se_dir)
  cov_dir <- meta$cp_true > ci_dir[, 1] & meta$cp_true < ci_dir[, 2]

  cp_cov_pct <- sapply(res, function(r) {
    v <- r$cov_cp_pct
    if (is.null(v)) NA else v
  })

  xc_rows[[paste0(sid, "_mce")]] <- data.frame(
    scenario = sid, method = "MCE",
    bias = round(mean(cp_bias, na.rm = TRUE), 6),
    ESE  = round(ese_cp, 6),
    ASE_delta  = round(mean(cp_se_d, na.rm = TRUE), 6),
    ASE_direct = round(mean(cp_se_dir, na.rm = TRUE), 6),
    SE_ratio_delta  = round(mean(cp_se_d, na.rm = TRUE) / ese_cp, 4),
    SE_ratio_direct = round(mean(cp_se_dir, na.rm = TRUE) / ese_cp, 4),
    ECP_delta  = round(mean(cp_cov_d, na.rm = TRUE), 4),
    ECP_direct = round(mean(cov_dir, na.rm = TRUE), 4),
    ECP_percentile = round(mean(cp_cov_pct, na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )

  for (method in c("spline", "quadratic", "quad_interact", "spline_interact")) {
    xc_vec  <- sf(res, method, "xc")
    bia_vec <- sf(res, method, "xc_bias")
    cov_vec <- sf(res, method, "xc_cov")
    if (all(is.na(xc_vec))) next
    xc_rows[[paste0(sid, "_", method)]] <- data.frame(
      scenario = sid, method = method,
      bias = round(mean(bia_vec, na.rm = TRUE), 6),
      ESE  = round(sd(xc_vec, na.rm = TRUE), 6),
      ASE_delta  = NA, ASE_direct = NA,
      SE_ratio_delta = NA, SE_ratio_direct = NA,
      ECP_delta = NA, ECP_direct = NA,
      ECP_percentile = round(mean(cov_vec, na.rm = TRUE), 4),
      stringsAsFactors = FALSE
    )
  }

  # ---- 4. C-index ----
  c_train_mce <- sapply(res, function(r) r$c_train)
  c_test_mce  <- sapply(res, function(r) r$c_test)
  cindex_rows[[paste0(sid, "_mce")]] <- data.frame(
    scenario = sid, method = "MCE",
    c_train = round(mean(c_train_mce, na.rm = TRUE), 4),
    c_train_sd = round(sd(c_train_mce, na.rm = TRUE), 4),
    c_test  = round(mean(c_test_mce,  na.rm = TRUE), 4),
    c_test_sd = round(sd(c_test_mce, na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )
  for (method in c("spline", "quadratic", "quad_interact", "spline_interact", "oracle")) {
    ct <- sf(res, method, "cindex_train")
    cv <- sf(res, method, "cindex")
    if (all(is.na(ct))) next
    cindex_rows[[paste0(sid, "_", method)]] <- data.frame(
      scenario = sid, method = method,
      c_train = round(mean(ct, na.rm = TRUE), 4),
      c_train_sd = round(sd(ct, na.rm = TRUE), 4),
      c_test  = round(mean(cv, na.rm = TRUE), 4),
      c_test_sd = round(sd(cv, na.rm = TRUE), 4),
      stringsAsFactors = FALSE
    )
  }

  # ---- 5. Survival S(t0) ----
  s_train_mce <- sapply(res, function(r) r$s_train)
  s_test_mce  <- sapply(res, function(r) r$s_test)
  surv_rows[[paste0(sid, "_mce")]] <- data.frame(
    scenario = sid, method = "MCE",
    s_t0_train = round(mean(s_train_mce, na.rm = TRUE), 4),
    s_t0_test  = round(mean(s_test_mce,  na.rm = TRUE), 4),
    s_t0_true  = tv$s_t0,
    bias_train = round(mean(s_train_mce, na.rm = TRUE) - tv$s_t0, 4),
    bias_test  = round(mean(s_test_mce,  na.rm = TRUE) - tv$s_t0, 4),
    stringsAsFactors = FALSE
  )
  for (method in c("spline", "quadratic", "quad_interact", "spline_interact", "oracle")) {
    s_survfit <- sf(res, method, "s_t0")
    s_km      <- sf(res, method, "s_t0_km")
    if (all(is.na(s_survfit)) && all(is.na(s_km))) next
    surv_rows[[paste0(sid, "_", method)]] <- data.frame(
      scenario = sid, method = method,
      s_t0_train = round(mean(s_km, na.rm = TRUE), 4),
      s_t0_test  = NA,
      s_t0_true  = tv$s_t0,
      bias_train = round(mean(s_km, na.rm = TRUE) - tv$s_t0, 4),
      bias_test  = NA,
      stringsAsFactors = FALSE
    )
  }

  # ---- 6. Critical Region ----
  crlb_mce <- sapply(res, function(r) r$est_crlb)
  crub_mce <- sapply(res, function(r) r$est_crub)
  cr_rows[[paste0(sid, "_mce")]] <- data.frame(
    scenario = sid, method = "MCE",
    crlb_mean = round(mean(crlb_mce, na.rm = TRUE), 4),
    crub_mean = round(mean(crub_mce, na.rm = TRUE), 4),
    crlb_bias = round(mean(crlb_mce, na.rm = TRUE) - tv$crlb, 4),
    crub_bias = round(mean(crub_mce, na.rm = TRUE) - tv$crub, 4),
    crlb_true = tv$crlb,
    crub_true = tv$crub,
    stringsAsFactors = FALSE
  )
  for (method in c("spline", "quadratic", "quad_interact", "spline_interact", "oracle")) {
    crlb_v <- sf(res, method, "crlb_km")
    crub_v <- sf(res, method, "crub_km")
    if (all(is.na(crlb_v)) && all(is.na(crub_v))) next
    cr_rows[[paste0(sid, "_", method)]] <- data.frame(
      scenario = sid, method = method,
      crlb_mean = round(mean(crlb_v, na.rm = TRUE), 4),
      crub_mean = round(mean(crub_v, na.rm = TRUE), 4),
      crlb_bias = round(mean(crlb_v, na.rm = TRUE) - tv$crlb, 4),
      crub_bias = round(mean(crub_v, na.rm = TRUE) - tv$crub, 4),
      crlb_true = tv$crlb,
      crub_true = tv$crub,
      stringsAsFactors = FALSE
    )
  }

  # ---- 7. Bootstrap diagnostics ----
  nboot_used_vec <- sapply(res, function(r) if (!is.null(r$nboot_used)) r$nboot_used else NA)
  nboot_run_vec  <- sapply(res, function(r) if (!is.null(r$nboot_run))  r$nboot_run  else NA)
  boot_method    <- res[[1]]$boot_method
  boot_rows[[sid]] <- data.frame(
    scenario = sid,
    boot_method      = if (is.null(boot_method)) NA else boot_method,
    nboot_used_mean  = round(mean(nboot_used_vec, na.rm = TRUE), 1),
    nboot_used_min   = min(nboot_used_vec, na.rm = TRUE),
    nboot_run        = if (length(unique(nboot_run_vec)) == 1) nboot_run_vec[1] else NA,
    stringsAsFactors = FALSE
  )

  cat(sprintf("  %s done (%d reps)\n", sid, nrep))
}

# Write all CSVs
write_table <- function(rows, fname) {
  if (length(rows) == 0) return(invisible(NULL))
  df <- do.call(rbind, rows)
  path <- file.path(outdir, fname)
  write.csv(df, path, row.names = FALSE)
  cat(sprintf("Wrote %s (%d rows)\n", path, nrow(df)))
}

write_table(meta_rows,    "sim_meta.csv")
write_table(param_rows,   "sim_params.csv")
write_table(xc_rows,      "sim_xc.csv")
write_table(cindex_rows,  "sim_cindex.csv")
write_table(surv_rows,    "sim_survival.csv")
write_table(cr_rows,      "sim_cr.csv")
write_table(boot_rows,    "sim_boot_diag.csv")

cat(sprintf("\nDone. All CSVs written to %s\n", outdir))
