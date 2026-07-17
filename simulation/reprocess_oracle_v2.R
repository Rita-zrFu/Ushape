#!/usr/bin/env Rscript
# reprocess_oracle_v2.R
#
# Recompute oracle Cox with the corrected Z-specific 7-parameter specification:
#   - Per-subject true critical point X_c(Z_i) = (alpha1*Z1_i - alpha2*Z2_i - beta0) / (1+beta1)
#   - Branch indicator B_R added to the Cox formula (7 parameters total)
#   - Under min-EV error and beta0=0, this parameterization spans the true Cox PH
#     linear predictor exactly; under general beta0 it remains correctly specified
#     once a B_R coefficient is allowed.
#
# Loads existing .rds files, regenerates training/test data from the saved seed
# (using identical logic to simulation_one_rep.R), refits ONLY the oracle with
# the new spec, and stores the result as `r$competitors$oracle_v2` while
# preserving the original `r$competitors$oracle` field for audit.
#
# Does NOT touch MCE, bootstrap, or any other competitor.
#
# Usage:
#   Rscript reprocess_oracle_v2.R <SCENARIO_ID> [RESULTS_DIR]
# Example:
#   Rscript reprocess_oracle_v2.R S01 results
#   Rscript reprocess_oracle_v2.R S08 /path/to/results

suppressPackageStartupMessages({
  library(survival)
  library(SurvMetrics)  # loaded for backward-compat only; not used for Cindex
})

# Fast C-index via survival::concordance — identical to SurvMetrics::Cindex but
# ~3000x faster for large n_test. Convention: higher `marker` = longer survival
# (concordant). Equivalent to SurvMetrics::Cindex(surv_obj, marker).
fast_cindex <- function(surv_obj, marker) {
  survival::concordance(surv_obj ~ marker)$concordance
}

args <- commandArgs(trailingOnly = TRUE)
sid     <- if (length(args) >= 1) args[1] else "S01"
basedir <- if (length(args) >= 2) args[2] else "results_local"

# === Fixed DGM parameters (must match simulation_one_rep.R) ===
beta0  <- 0
alpha1 <- 3
alpha2 <- -3
x_lo   <- -5
x_hi   <- 10
z1_ref <- 0.5
z2_ref <- 1
n_test <- 5000
euler  <- 0.5772156649
gumbel_sd <- pi / sqrt(6)

# === Helpers (match simulation_one_rep.R) ===
H_fun <- function(x_vec, z1_vec, z2_vec, theta) {
  pmax(-x_vec + theta[3] * z1_vec,
       theta[1] + exp(theta[2]) * x_vec + theta[4] * z2_vec)
}

# Kernel-weighted KM on linear predictor (matches simulation_one_rep.R:287)
km_S_on_lp <- function(lp_all, y_all, d_all, lp_target, t0_val, bw) {
  w <- exp(-0.5 * ((lp_all - lp_target) / bw)^2)
  event_times <- sort(unique(y_all[d_all == 1]))
  event_times <- event_times[event_times <= t0_val]
  if (length(event_times) == 0) return(1)
  S <- 1
  for (yk in event_times) {
    at_risk <- sum(w * (y_all >= yk))
    events  <- sum(w * (y_all == yk) * d_all)
    if (at_risk > 0) S <- S * (1 - events / at_risk)
  }
  return(max(S, 0))
}

# === Locate .rds files ===
sdir  <- file.path(basedir, sid)
files <- list.files(sdir, pattern = "[.]rds$", full.names = TRUE)
if (length(files) == 0) {
  stop(sprintf("No .rds files found in %s", sdir))
}
cat(sprintf("Reprocessing oracle for %s: %d files in %s\n", sid, length(files), sdir))

n_done <- 0; n_fail <- 0; n_skip <- 0

for (i in seq_along(files)) {
  r <- tryCatch(readRDS(files[i]), error = function(e) NULL)
  if (is.null(r)) { n_skip <- n_skip + 1; next }

  # Scenario-level parameters from the saved rec
  n        <- r$n
  ei       <- r$ei
  G_type   <- r$G_type
  seed     <- r$seed
  cl       <- r$censor_lo
  cu       <- r$censor_hi
  t0_S     <- r$t0_S
  t0_CR    <- r$t0_CR
  # `beta1` is not stored as its own field; recover from theta_true = (beta0, log(beta1), alpha1, alpha2)
  theta_true <- r$theta_true
  beta1      <- exp(theta_true[2])
  cp_true    <- (alpha1 * z1_ref - alpha2 * z2_ref - beta0) / (1 + beta1)
  G_fun      <- if (G_type == "logistic") function(s) 10 * plogis(s, scale = 5) else function(s) exp(s / 5)

  # === Regenerate training data (matches simulation_one_rep.R:99-122) ===
  set.seed(seed)
  z1 <- rnorm(n, 0, 1)
  z2 <- rbinom(n, 1, 0.5)
  x  <- runif(n, x_lo, x_hi)
  H  <- H_fun(x, z1, z2, theta_true)
  epsilon <- if (ei == "norm") {
    rnorm(n, 0, 3)
  } else {
    -((-log(-log(runif(n))) - euler) * (3 / gumbel_sd))
  }
  T_surv <- G_fun(-H + epsilon)
  censor <- runif(n, cl, cu)
  y     <- pmin(T_surv, censor)
  delta <- as.numeric(T_surv < censor)

  # === Regenerate test data (matches simulation_one_rep.R:322-338) ===
  set.seed(seed + 100000)
  z1_test <- rnorm(n_test); z2_test <- rbinom(n_test, 1, 0.5)
  x_test  <- runif(n_test, x_lo, x_hi)
  H_test  <- H_fun(x_test, z1_test, z2_test, theta_true)
  eps_test <- if (ei == "norm") {
    rnorm(n_test, 0, 3)
  } else {
    -((-log(-log(runif(n_test))) - euler) * (3 / gumbel_sd))
  }
  T_test      <- G_fun(-H_test + eps_test)
  censor_test <- runif(n_test, cl, cu)
  y_test      <- pmin(T_test, censor_test)
  delta_test  <- as.numeric(T_test < censor_test)

  df_train <- data.frame(X = x,      Z1 = z1,      Z2 = z2)
  df_test  <- data.frame(X = x_test, Z1 = z1_test, Z2 = z2_test)
  x_grid   <- seq(x_lo, x_hi, by = 0.25)

  # === Fit Z-specific 7-parameter oracle Cox ===
  ora_v2 <- tryCatch({
    # Per-subject critical points (train)
    cp_true_i_train <- (alpha1 * df_train$Z1 - alpha2 * df_train$Z2 - beta0) / (1 + beta1)
    B_L <- as.numeric(df_train$X < cp_true_i_train)
    B_R <- 1 - B_L
    X_L <- (df_train$X - cp_true_i_train) * B_L
    X_R <- (df_train$X - cp_true_i_train) * B_R
    df_ora <- data.frame(y = y, delta = delta,
      B_R = B_R,
      X_L = X_L, X_R = X_R,
      Z1_L = df_train$Z1 * B_L, Z1_R = df_train$Z1 * B_R,
      Z2_L = df_train$Z2 * B_L, Z2_R = df_train$Z2 * B_R)

    fit_ora <- coxph(Surv(y, delta) ~ B_R + X_L + X_R + Z1_L + Z1_R + Z2_L + Z2_R, data = df_ora)

    # Train C-index
    lp_ora_train <- predict(fit_ora, type = "lp")
    c_ora_train  <- fast_cindex(Surv(y, delta), -lp_ora_train)

    # Test C-index (per-subject X_c on the test set)
    cp_true_i_test <- (alpha1 * df_test$Z1 - alpha2 * df_test$Z2 - beta0) / (1 + beta1)
    B_L_te <- as.numeric(df_test$X < cp_true_i_test)
    B_R_te <- 1 - B_L_te
    df_ora_te <- data.frame(
      B_R = B_R_te,
      X_L = (df_test$X - cp_true_i_test) * B_L_te,
      X_R = (df_test$X - cp_true_i_test) * B_R_te,
      Z1_L = df_test$Z1 * B_L_te, Z1_R = df_test$Z1 * B_R_te,
      Z2_L = df_test$Z2 * B_L_te, Z2_R = df_test$Z2 * B_R_te)
    lp_ora_test <- predict(fit_ora, newdata = df_ora_te, type = "lp")
    c_ora_test  <- fast_cindex(Surv(y_test, delta_test), -lp_ora_test)

    # Survfit-based S(t0) at reference subject (Z=z_ref, X=1; left branch since 1<cp_true=1.5 or 2.25)
    ref_ora <- data.frame(
      B_R = 0,
      X_L = (1 - cp_true), X_R = 0,
      Z1_L = z1_ref, Z1_R = 0,
      Z2_L = z2_ref, Z2_R = 0)
    s_ora <- tryCatch({
      sf <- survfit(fit_ora, newdata = ref_ora)
      idx <- max(which(sf$time <= t0_S)); sf$surv[idx]
    }, error = function(e) NA)

    # Survfit-based CR (over X grid, reference Z fixed)
    cr_ora <- tryCatch({
      s_g <- sapply(x_grid, function(xv) {
        bl <- as.numeric(xv < cp_true); br <- 1 - bl
        nd <- data.frame(
          B_R = br,
          X_L = (xv - cp_true) * bl, X_R = (xv - cp_true) * br,
          Z1_L = z1_ref * bl, Z1_R = z1_ref * br,
          Z2_L = z2_ref * bl, Z2_R = z2_ref * br)
        sf <- survfit(fit_ora, newdata = nd)
        idx <- max(which(sf$time <= t0_CR), 0)
        if (idx == 0) 1 else sf$surv[idx]
      })
      ic <- !is.na(s_g) & s_g >= 0.80
      if (any(ic)) list(lo = min(x_grid[ic]), hi = max(x_grid[ic])) else list(lo = NA, hi = NA)
    }, error = function(e) list(lo = NA, hi = NA))

    # KM-based S(t0) and CR
    bw_ora <- 1.06 * sd(lp_ora_train) * n^(-1/5)
    lp_ref_ora <- predict(fit_ora, newdata = ref_ora, type = "lp")
    s_ora_km <- km_S_on_lp(lp_ora_train, y, delta, lp_ref_ora, t0_S, bw_ora)
    cr_ora_km <- tryCatch({
      s_g <- sapply(x_grid, function(xv) {
        bl <- as.numeric(xv < cp_true); br <- 1 - bl
        nd <- data.frame(
          B_R = br,
          X_L = (xv - cp_true) * bl, X_R = (xv - cp_true) * br,
          Z1_L = z1_ref * bl, Z1_R = z1_ref * br,
          Z2_L = z2_ref * bl, Z2_R = z2_ref * br)
        lp_x <- predict(fit_ora, newdata = nd, type = "lp")
        km_S_on_lp(lp_ora_train, y, delta, lp_x, t0_CR, bw_ora)
      })
      ic <- !is.na(s_g) & s_g >= 0.80
      if (any(ic)) list(lo = min(x_grid[ic]), hi = max(x_grid[ic])) else list(lo = NA, hi = NA)
    }, error = function(e) list(lo = NA, hi = NA))

    list(
      xc = cp_true, xc_bias = 0, xc_ci = c(NA, NA), xc_cov = NA,
      cindex_train = c_ora_train, cindex = c_ora_test,
      s_t0 = s_ora, crlb = cr_ora$lo, crub = cr_ora$hi,
      s_t0_km = s_ora_km, crlb_km = cr_ora_km$lo, crub_km = cr_ora_km$hi,
      spec = "z_specific_7param_with_BR")
  }, error = function(e) {
    cat(sprintf("  rep %d (seed %d) oracle_v2 FAILED: %s\n", i, seed, e$message))
    list(xc = NA, xc_bias = NA, xc_ci = c(NA, NA), xc_cov = NA,
         cindex_train = NA, cindex = NA, s_t0 = NA, crlb = NA, crub = NA,
         s_t0_km = NA, crlb_km = NA, crub_km = NA,
         spec = "z_specific_7param_with_BR")
  })

  # Attach as a new field; preserve original `oracle`
  if (is.null(r$competitors)) r$competitors <- list()
  r$competitors$oracle_v2 <- ora_v2

  saveRDS(r, files[i])
  n_done <- n_done + 1
  if (i %% 100 == 0) cat(sprintf("  %d/%d done (failures so far: %d)\n", i, length(files), n_fail))
}

cat(sprintf("Done. Reprocessed %d files (skipped %d, failed %d) for scenario %s.\n",
            n_done, n_skip, n_fail, sid))
