#!/usr/bin/env Rscript
# summarize_all_to_csv.R — Summarize all 23 simulation scenarios to CSV
# Usage: Rscript summarize_all_to_csv.R [results_dir] [output_dir]
# Output: <output_dir>/*.csv  (default: csv_output)

args <- commandArgs(trailingOnly = TRUE)
basedir <- if (length(args) >= 1) args[1] else "results"
outdir  <- if (length(args) >= 2) args[2] else "csv_output"
dir.create(outdir, showWarnings = FALSE)

# ---- Scenario definitions ----
scenarios <- paste0("S", sprintf("%02d", 1:23))

# True values by (G, epsilon, beta1) — same as summarize_results.R
true_vals <- list(
  "logistic_norm_2" = list(s_t0 = 0.4341, crlb = -2.90, crub = 3.70),
  "logistic_ev_2"   = list(s_t0 = 0.4988, crlb = -3.27, crub = 3.88),
  "exp_norm_2"      = list(s_t0 = 0.5670, crlb = -3.00, crub = 3.75),
  "exp_ev_2"        = list(s_t0 = 0.6357, crlb = -3.35, crub = 3.92),
  "logistic_ev_1"   = list(s_t0 = 0.4990, crlb = -3.27, crub = 7.77)
)

pnames <- c("beta0", "log_b1", "alpha1", "alpha2")

# ---- Accumulators ----
meta_rows    <- list()
param_rows   <- list()
xc_rows      <- list()
cindex_rows  <- list()
surv_rows    <- list()
cr_rows      <- list()
boot_rows    <- list()

# ---- Helper: safely extract competitor field ----
sf <- function(res_list, method, field) {
  sapply(res_list, function(r) {
    v <- r$competitors[[method]][[field]]
    if (is.null(v)) NA else v
  })
}

# ---- Process each scenario ----
for (sid in scenarios) {
  files <- list.files(file.path(basedir, sid), pattern = "[.]rds$", full.names = TRUE)
  if (length(files) == 0) {
    cat(sprintf("%s: no files found, skipping\n", sid))
    next
  }
  files <- files[1:min(1000, length(files))]
  cat(sprintf("Reading %s: %d files ...\n", sid, length(files)))
  res  <- lapply(files, readRDS)
  meta <- res[[1]]
  nrep <- length(res)

  # ---- Provenance checks (same as summarize_results.R) ----
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
  required_fields <- c("boot_cp_ci_pct", "cov_cp_pct", "beta1")
  required_comp   <- c("spline_interact")
  n_old <- sum(sapply(res, function(r) {
    any(!required_fields %in% names(r)) ||
      any(!required_comp %in% names(r$competitors))
  }))
  if (n_old > 0) {
    cat(sprintf("  *** %s: %d/%d files missing expected fields ***\n", sid, n_old, nrep))
    prov_ok <- FALSE
  }
  if (!prov_ok) {
    cat(sprintf("  *** SKIPPING %s ***\n", sid))
    next
  }

  # ---- Lookup true values ----
  b1  <- unique(b1_vals)
  key <- paste0(meta$G_type, "_", meta$ei, "_", b1)
  tv  <- true_vals[[key]]
  if (is.null(tv)) {
    cat(sprintf("%s: no true values for key '%s', skipping\n", sid, key))
    next
  }

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

  # Direct bootstrap CI coverage
  ci_dir  <- cbind(cp_vec - qnorm(0.975) * cp_se_dir,
                   cp_vec + qnorm(0.975) * cp_se_dir)
  cov_dir <- meta$cp_true > ci_dir[, 1] & meta$cp_true < ci_dir[, 2]

  # Percentile CI coverage
  cp_cov_pct <- sapply(res, function(r) {
    v <- r$cov_cp_pct
    if (is.null(v)) NA else v
  })

  # MCE row
  xc_rows[[paste0(sid, "_MCE")]] <- data.frame(
    scenario     = sid,
    method       = "MCE",
    xc_mean      = round(mean(cp_vec, na.rm = TRUE), 6),
    bias         = round(mean(cp_bias, na.rm = TRUE), 6),
    RMSE         = round(sqrt(mean(cp_bias^2, na.rm = TRUE)), 6),
    ESE          = round(ese_cp, 6),
    CI_cov       = round(mean(cp_cov_d, na.rm = TRUE), 4),
    CI_cov_delta   = round(mean(cp_cov_d, na.rm = TRUE), 4),
    CI_cov_direct  = round(mean(cov_dir, na.rm = TRUE), 4),
    CI_cov_pct     = round(mean(cp_cov_pct, na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )

  # Competitor Xc rows
  comp_names <- names(res[[1]]$competitors)
  for (m in comp_names) {
    xb    <- sf(res, m, "xc_bias")
    xv    <- sapply(res, function(r) r$competitors[[m]]$xc)
    xc    <- sf(res, m, "xc_cov")
    ese_m <- sd(xv, na.rm = TRUE)
    xc_rows[[paste0(sid, "_", m)]] <- data.frame(
      scenario     = sid,
      method       = m,
      xc_mean      = round(mean(xv, na.rm = TRUE), 6),
      bias         = round(mean(xb, na.rm = TRUE), 6),
      RMSE         = round(sqrt(mean(xb^2, na.rm = TRUE)), 6),
      ESE          = round(ese_m, 6),
      CI_cov       = round(mean(as.numeric(xc), na.rm = TRUE), 4),
      CI_cov_delta   = NA_real_,
      CI_cov_direct  = NA_real_,
      CI_cov_pct     = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  # ---- 4. C-index ----
  c_tr <- sapply(res, function(r) r$c_train)
  c_te <- sapply(res, function(r) r$c_test)
  cindex_rows[[paste0(sid, "_MCE")]] <- data.frame(
    scenario     = sid,
    method       = "MCE",
    c_train_mean = round(mean(c_tr, na.rm = TRUE), 6),
    c_train_sd   = round(sd(c_tr, na.rm = TRUE), 6),
    c_test_mean  = round(mean(c_te, na.rm = TRUE), 6),
    c_test_sd    = round(sd(c_te, na.rm = TRUE), 6),
    stringsAsFactors = FALSE
  )
  for (m in comp_names) {
    ct <- sf(res, m, "cindex_train")
    ce <- sf(res, m, "cindex")
    cindex_rows[[paste0(sid, "_", m)]] <- data.frame(
      scenario     = sid,
      method       = m,
      c_train_mean = round(mean(ct, na.rm = TRUE), 6),
      c_train_sd   = round(sd(ct, na.rm = TRUE), 6),
      c_test_mean  = round(mean(ce, na.rm = TRUE), 6),
      c_test_sd    = round(sd(ce, na.rm = TRUE), 6),
      stringsAsFactors = FALSE
    )
  }

  # ---- 5. Survival S(t0|X=1) ----
  s_mce <- sapply(res, function(r) r$s_test)
  surv_rows[[paste0(sid, "_MCE")]] <- data.frame(
    scenario = sid,
    method   = "MCE",
    s_t0_km  = round(mean(s_mce, na.rm = TRUE), 6),
    bias_km  = round(mean(s_mce, na.rm = TRUE) - tv$s_t0, 6),
    s_t0_sf  = NA_real_,
    bias_sf  = NA_real_,
    stringsAsFactors = FALSE
  )
  for (m in comp_names) {
    st_sf <- sf(res, m, "s_t0")
    st_km <- sf(res, m, "s_t0_km")
    surv_rows[[paste0(sid, "_", m)]] <- data.frame(
      scenario = sid,
      method   = m,
      s_t0_km  = if (all(is.na(st_km))) NA_real_ else round(mean(st_km, na.rm = TRUE), 6),
      bias_km  = if (all(is.na(st_km))) NA_real_ else round(mean(st_km, na.rm = TRUE) - tv$s_t0, 6),
      s_t0_sf  = if (all(is.na(st_sf))) NA_real_ else round(mean(st_sf, na.rm = TRUE), 6),
      bias_sf  = if (all(is.na(st_sf))) NA_real_ else round(mean(st_sf, na.rm = TRUE) - tv$s_t0, 6),
      stringsAsFactors = FALSE
    )
  }

  # ---- 6. CR Bounds (KM-based) ----
  crlb_mce <- sapply(res, function(r) r$est_crlb)
  crub_mce <- sapply(res, function(r) r$est_crub)
  cr_rows[[paste0(sid, "_MCE")]] <- data.frame(
    scenario  = sid,
    method    = "MCE",
    crlb_mean = round(mean(crlb_mce, na.rm = TRUE), 4),
    crlb_bias = round(mean(crlb_mce, na.rm = TRUE) - tv$crlb, 4),
    crlb_rmse = round(sqrt(mean((crlb_mce - tv$crlb)^2, na.rm = TRUE)), 4),
    crub_mean = round(mean(crub_mce, na.rm = TRUE), 4),
    crub_bias = round(mean(crub_mce, na.rm = TRUE) - tv$crub, 4),
    crub_rmse = round(sqrt(mean((crub_mce - tv$crub)^2, na.rm = TRUE)), 4),
    n_nonNA   = sum(!is.na(crlb_mce)),
    stringsAsFactors = FALSE
  )
  for (m in comp_names) {
    cl_km <- sf(res, m, "crlb_km")
    cu_km <- sf(res, m, "crub_km")
    cr_rows[[paste0(sid, "_", m)]] <- data.frame(
      scenario  = sid,
      method    = m,
      crlb_mean = if (all(is.na(cl_km))) NA_real_ else round(mean(cl_km, na.rm = TRUE), 4),
      crlb_bias = if (all(is.na(cl_km))) NA_real_ else round(mean(cl_km, na.rm = TRUE) - tv$crlb, 4),
      crlb_rmse = if (all(is.na(cl_km))) NA_real_ else round(sqrt(mean((cl_km - tv$crlb)^2, na.rm = TRUE)), 4),
      crub_mean = if (all(is.na(cu_km))) NA_real_ else round(mean(cu_km, na.rm = TRUE), 4),
      crub_bias = if (all(is.na(cu_km))) NA_real_ else round(mean(cu_km, na.rm = TRUE) - tv$crub, 4),
      crub_rmse = if (all(is.na(cu_km))) NA_real_ else round(sqrt(mean((cu_km - tv$crub)^2, na.rm = TRUE)), 4),
      n_nonNA   = sum(!is.na(cl_km)),
      stringsAsFactors = FALSE
    )
  }

  # ---- 7. Bootstrap diagnostics ----
  nboot_used <- sapply(res, function(r) r$nboot_used)
  nboot_run  <- res[[1]]$nboot_run
  boot_rows[[sid]] <- data.frame(
    scenario        = sid,
    nboot_used_mean = round(mean(nboot_used, na.rm = TRUE), 1),
    nboot_run       = nboot_run,
    convergence_rate = round(mean(nboot_used, na.rm = TRUE) / nboot_run, 4),
    stringsAsFactors = FALSE
  )

  cat(sprintf("  %s done (%d reps)\n", sid, nrep))
}

# ---- Write CSVs ----
write_csv <- function(rows, fname) {
  df <- do.call(rbind, rows)
  f  <- file.path(outdir, fname)
  write.csv(df, f, row.names = FALSE)
  cat(sprintf("Wrote %s (%d rows)\n", f, nrow(df)))
}

write_csv(meta_rows,   "sim_meta.csv")
write_csv(param_rows,  "sim_params.csv")
write_csv(xc_rows,     "sim_xc.csv")
write_csv(cindex_rows, "sim_cindex.csv")
write_csv(surv_rows,   "sim_survival.csv")
write_csv(cr_rows,     "sim_cr.csv")
write_csv(boot_rows,   "sim_boot_diag.csv")

cat("\nDone. All CSVs written to ", outdir, "\n")
