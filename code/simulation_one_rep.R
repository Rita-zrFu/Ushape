#!/usr/bin/env Rscript
# simulation_one_rep.R
# One replicate of the USR simulation study (companion to the manuscript).
#
# Data-generating mechanism (see Sections 2 and 6 of the manuscript):
#   T = G(-H + epsilon), H(X, Z) = max(-X + alpha1*Z1, beta0 + beta1*X + alpha2*Z2)
#   X ~ U(-5, 10), Z1 ~ N(0, 1), Z2 ~ Bernoulli(0.5)
#   alpha1 = 3, alpha2 = -3, beta0 = 0, beta1 in {1, 2}
#   epsilon: N(0, 9) ('norm') or centered min-extreme-value with sd=3 ('ev')
#   G(s): logistic (10 * plogis(s, scale = 5)) or exponential (exp((s+1)/5))
#   Censoring: independent U(censor_lo, censor_hi); bounds set per scenario
#
# Estimation:
#   Point: DEoptim on the hard C-index objective (compiled via Rcpp)
#   Bootstrap: L-BFGS-B on the sigmoid-smoothed C-index
#     sigma_h = round(n^(1/4))  ->  {4, 5, 6} for n in {200, 500, 1000}
#     525 bootstrap runs, first 500 converged replicates used
#   Xc (critical point) standard error: delta method on bootstrap covariance
#   Xc percentile CI: 2.5%/97.5% bootstrap quantiles (saved alongside Wald)
#
# Competitors (when run_competitors = TRUE):
#   Cox + natural spline; Cox + quadratic; Cox + quadratic with X*Z interactions;
#   Cox + spline with interactions; oracle piecewise Cox at the true Xc.
#
# Usage:
#   Rscript simulation_one_rep.R <scenario_id> <seed> <n> <ei> <G_type>
#                                 <censor_lo> <censor_hi> <run_competitors> [beta1]
#
# Example (S07, seed 1, n=1000, normal errors, logistic G, beta1=2):
#   Rscript simulation_one_rep.R S07 1 1000 norm logistic 0.3000 8.3000 TRUE 2
#
# Output:
#   ${RESULTS_DIR:-results}/<scenario_id>/rep_<seed>_j<batch_tag>.rds
#   (One R list per replicate; see SAVE RESULTS section at the bottom for fields.)
#
# Environment:
#   MYCINDEX_CPP_PATH : absolute path to myCindex.cpp; defaults to ../code/myCindex.cpp
#   RESULTS_DIR       : base directory for output (default: results)
#   SLURM_ARRAY_JOB_ID: optional; used as batch_tag for provenance

# === PARSE ARGUMENTS ===
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop("Usage: Rscript simulation_one_rep.R scenario_id seed n ei G_type censor_lo censor_hi run_competitors [beta1]")
}
scenario_id     <- args[1]
seed            <- as.integer(args[2])
n               <- as.integer(args[3])
ei              <- args[4]       # "norm" or "ev"
G_type          <- args[5]       # "logistic" or "exp"
censor_lo       <- as.numeric(args[6])
censor_hi       <- as.numeric(args[7])
run_competitors <- as.logical(args[8])
beta1_arg       <- if (length(args) >= 9) as.numeric(args[9]) else 1

cat(sprintf("=== %s | seed=%d | n=%d | ei=%s | G=%s | cens=U(%.4f,%.4f) | comp=%s | b1=%s ===\n",
            scenario_id, seed, n, ei, G_type, censor_lo, censor_hi, run_competitors, beta1_arg))
cat(sprintf("Start: %s\n", Sys.time()))

# === PACKAGES ===
suppressPackageStartupMessages({
  library(DEoptim)
  library(survival)
  library(SurvMetrics)
  library(Rcpp)
  library(splines)
})

# Compile C++ kernel (myCindex + myCindex_smooth).
# Resolution order: (1) MYCINDEX_CPP_PATH env var, (2) ./myCindex.cpp (same dir as this script),
# (3) ./code/myCindex.cpp (relative to working directory).
cpp_candidates <- c(
  Sys.getenv("MYCINDEX_CPP_PATH", unset = ""),
  "myCindex.cpp",
  file.path("code", "myCindex.cpp")
)
cpp_path <- cpp_candidates[file.exists(cpp_candidates)][1]
if (is.na(cpp_path)) {
  stop("myCindex.cpp not found. Set MYCINDEX_CPP_PATH or place file in code/.")
}
sourceCpp(cpp_path)

# === FIXED PARAMETERS ===
beta0  <- 0
beta1  <- beta1_arg
alpha1 <- 3
alpha2 <- -3
x_lo   <- -5
x_hi   <- 10

theta_true <- c(beta0, log(beta1), alpha1, alpha2)

# Reference covariates for Xc, S(t0), and CR evaluation
z1_ref <- 0.5
z2_ref <- 1
cp_true <- (alpha1 * z1_ref - alpha2 * z2_ref - beta0) / (1 + beta1)

# Survival evaluation time t0_S (S(t0_S|Xc) ~ 0.5);
# critical-region evaluation time t0_CR (S(t0_CR|Xc) ~ 0.9, threshold 0.80)
t0_S  <- if (G_type == "logistic") 5.0 else 1.0
t0_CR <- if (G_type == "logistic") 2.0 else 0.3

# G functions
G_fun <- if (G_type == "logistic") {
  function(s) 10 * plogis(s, scale = 5)
} else {
  function(s) exp((s + 1) / 5)
}

# Extreme value helpers (used for min-EV epsilon generation)
euler <- 0.5772156649
gumbel_sd <- pi / sqrt(6)

# H function
H_fun <- function(x_vec, z1_vec, z2_vec, theta) {
  pmax(-x_vec + theta[3] * z1_vec,
       theta[1] + exp(theta[2]) * x_vec + theta[4] * z2_vec)
}

# === SECTION 1: DATA GENERATION ===
set.seed(seed)

z1 <- rnorm(n, 0, 1)
z2 <- rbinom(n, 1, 0.5)
x  <- runif(n, x_lo, x_hi)

H <- H_fun(x, z1, z2, theta_true)

epsilon <- if (ei == "norm") {
  rnorm(n, 0, 3)
} else if (ei == "ev") {
  # Minimum extreme value (centered, sd=3). PH holds for any G.
  -((-log(-log(runif(n))) - euler) * (3 / gumbel_sd))
} else {
  stop("Unknown ei: ", ei)
}

T_surv <- G_fun(-H + epsilon)

censor <- runif(n, censor_lo, censor_hi)
y     <- pmin(T_surv, censor)
delta <- as.numeric(T_surv < censor)
censor_rate <- 1 - mean(delta)

mymat <- cbind(x, z1, z2, y, delta)

# === SECTION 2: MCE POINT ESTIMATION (DEoptim on hard C-index) ===
# myCindex is used ONLY for optimization. For reporting C-index, use SurvMetrics::Cindex.

cindex_DE <- function(theta) {
  htemp <- H_fun(mymat[,1], mymat[,2], mymat[,3], theta)
  return(-myCindex(htemp, mymat[,4], mymat[,5]))
}

controlDE <- list(reltol = 1e-6, steptol = 100, itermax = 1000,
                  trace = FALSE, parallelType = 0)

fit_DE <- DEoptim(fn = cindex_DE,
                  lower = theta_true - 5, upper = theta_true + 5,
                  control = controlDE)

est_par    <- as.numeric(fit_DE$optim$bestmem)
est_cindex <- -fit_DE$optim$bestval

est_cpoint <- (-est_par[1] + est_par[3] * z1_ref - est_par[4] * z2_ref) /
              (1 + exp(est_par[2]))

bias    <- est_par - theta_true
bias_cp <- est_cpoint - cp_true

cat(sprintf("MCE complete. C=%.4f, Xc_hat=%.3f (true=%.3f)\n",
            est_cindex, est_cpoint, cp_true))

# === SECTION 3: MCE BOOTSTRAP — Smooth C-index + L-BFGS-B ===
# Sigmoid-smoothed myCindex with sigma_h = round(n^(1/4)) -> {4, 5, 6} for n in {200, 500, 1000}.
# L-BFGS-B exploits the smooth gradient landscape; warm-started at the DEoptim estimate.
# Satisfies Cattaneo, Jansson & Nagasawa (2020) consistency condition for cube-root smoothed
# bootstrap.

nboot_target <- 500
nboot_run    <- 525   # extra runs guard against non-convergence
bw      <- 3.0
sigma_h <- round(n^(1/4))
lower_b <- est_par - bw
upper_b <- est_par + bw

boot_pars_all <- matrix(NA, nboot_run, length(est_par))

for (b in 1:nboot_run) {
  set.seed(b)
  idx <- sample(n, replace = TRUE)

  obj_smooth <- function(theta) {
    htemp <- H_fun(mymat[idx, 1], mymat[idx, 2], mymat[idx, 3], theta)
    -myCindex_smooth(htemp, mymat[idx, 4], mymat[idx, 5], sigma_h = sigma_h)
  }

  fit_b <- tryCatch(
    optim(par = est_par, fn = obj_smooth, method = "L-BFGS-B",
          lower = lower_b, upper = upper_b,
          control = list(maxit = 500)),
    error = function(e) list(par = rep(NA, length(est_par)))
  )
  boot_pars_all[b, ] <- fit_b$par
}

boot_valid_idx <- which(complete.cases(boot_pars_all))
if (length(boot_valid_idx) >= nboot_target) {
  boot_pars <- boot_pars_all[boot_valid_idx[1:nboot_target], ]
} else {
  boot_pars <- boot_pars_all[boot_valid_idx, ]
  cat(sprintf("WARNING: only %d/%d bootstrap replicates converged\n",
              length(boot_valid_idx), nboot_target))
}
nboot_used <- nrow(boot_pars)

boot_se  <- apply(boot_pars, 2, sd)
boot_cov <- cov(boot_pars)

# Xc SE via delta method.
# Xc = (-theta1 + theta3 * z1_ref - theta4 * z2_ref) / (1 + exp(theta2))
denom <- 1 + exp(est_par[2])
numer <- -est_par[1] + est_par[3] * z1_ref - est_par[4] * z2_ref
grad_xc <- c(
  -1 / denom,                              # d/d(beta0)
  -numer * exp(est_par[2]) / denom^2,      # d/d(log beta1)
  z1_ref / denom,                          # d/d(alpha1)
  -z2_ref / denom                          # d/d(alpha2)
)
boot_cp_se_delta <- as.numeric(sqrt(t(grad_xc) %*% boot_cov %*% grad_xc))

# Direct bootstrap SE for Xc (saved for comparison; CI uses the delta-method SE)
boot_cp <- apply(boot_pars, 1, function(p) {
  (-p[1] + p[3] * z1_ref - p[4] * z2_ref) / (1 + exp(p[2]))
})
boot_cp_se_direct <- sd(boot_cp, na.rm = TRUE)
boot_cp_se <- boot_cp_se_delta

cat(sprintf("Xc SE: delta=%.4f, direct=%.4f\n", boot_cp_se_delta, boot_cp_se_direct))

# Wald CIs (primary)
ci_par <- cbind(est_par - qnorm(0.975) * boot_se,
                est_par + qnorm(0.975) * boot_se)
ci_cp  <- c(est_cpoint - qnorm(0.975) * boot_cp_se,
             est_cpoint + qnorm(0.975) * boot_cp_se)

# Bootstrap percentile CI for Xc (saved for sensitivity, not reported as primary)
boot_cp_ci_pct <- quantile(boot_cp, c(0.025, 0.975), na.rm = TRUE)
cov_cp_pct <- cp_true > boot_cp_ci_pct[1] & cp_true < boot_cp_ci_pct[2]

cov_par <- theta_true > ci_par[, 1] & theta_true < ci_par[, 2]
cov_cp  <- cp_true > ci_cp[1] & cp_true < ci_cp[2]

cat(sprintf("Bootstrap complete. boot_SE(Xc)=%.4f\n", boot_cp_se))

# === SECTION 4: S(t0) AND CR ESTIMATION (training set) ===
# Kernel-weighted Kaplan-Meier (see manuscript Section 4).
# S_h(t | X, Z) = prod_{y_k <= t} { 1 - weighted_events_k / weighted_at_risk_k }
# Weights: Gaussian kernel on H values; bandwidth via Silverman's rule.

est_H_train <- H_fun(x, z1, z2, est_par)
bw_h <- 1.06 * sd(est_H_train) * length(est_H_train)^(-1/5)

estimate_S_at_X <- function(x_val, z1_val, z2_val, H_all, y_all, d_all, t0_val, bw) {
  H_target <- H_fun(x_val, z1_val, z2_val, est_par)
  w <- exp(-0.5 * ((H_all - H_target) / bw)^2)
  event_times <- sort(unique(y_all[d_all == 1]))
  event_times <- event_times[event_times <= t0_val]
  if (length(event_times) == 0) return(1)
  S <- 1
  for (yk in event_times) {
    at_risk <- w * (y_all >= yk)
    events  <- w * (y_all == yk) * d_all
    if (sum(at_risk) > 0) S <- S * (1 - sum(events) / sum(at_risk))
  }
  max(S, 0)
}

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
  max(S, 0)
}

s_train <- estimate_S_at_X(1, z1_ref, z2_ref, est_H_train, y, delta, t0_S, bw_h)

x_cr_grid <- seq(x_lo, x_hi, by = 0.25)
s_cr_grid <- sapply(x_cr_grid, function(xv) {
  estimate_S_at_X(xv, z1_ref, z2_ref, est_H_train, y, delta, t0_CR, bw_h)
})
in_cr <- s_cr_grid >= 0.80
if (any(in_cr)) {
  est_crlb <- min(x_cr_grid[in_cr])
  est_crub <- max(x_cr_grid[in_cr])
} else {
  est_crlb <- NA
  est_crub <- NA
}

# === SECTION 5: HOLDOUT TEST SET ===
n_test <- 5000
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
censor_test <- runif(n_test, censor_lo, censor_hi)
y_test      <- pmin(T_test, censor_test)
delta_test  <- as.numeric(T_test < censor_test)

h_est_test <- H_fun(x_test, z1_test, z2_test, est_par)
c_test     <- Cindex(Surv(y_test, delta_test), -h_est_test)

est_H_eval <- H_fun(x, z1, z2, est_par)
c_train    <- Cindex(Surv(y, delta), -est_H_eval)

bw_h_test <- 1.06 * sd(h_est_test) * length(h_est_test)^(-1/5)
s_test <- estimate_S_at_X(1, z1_ref, z2_ref, h_est_test, y_test, delta_test, t0_S, bw_h_test)

# === SECTION 6: COMPETITOR METHODS ===
comp_results <- list()

if (run_competitors) {
  df_train <- data.frame(y = y, delta = delta, X = x, Z1 = z1, Z2 = z2)
  df_test  <- data.frame(y = y_test, delta = delta_test,
                         X = x_test, Z1 = z1_test, Z2 = z2_test)
  x_grid <- seq(x_lo, x_hi, by = 0.25)

  # ---- 6a. Cox + natural spline ----
  tryCatch({
    fit_sp <- coxph(Surv(y, delta) ~ ns(X, df = 4) + Z1 + Z2, data = df_train)
    pred_df <- data.frame(X = x_grid, Z1 = z1_ref, Z2 = z2_ref)
    lp_grid <- predict(fit_sp, newdata = pred_df, type = "lp")
    xc_sp <- x_grid[which.min(lp_grid)]

    lp_train_sp <- predict(fit_sp, type = "lp")
    c_sp_train <- Cindex(Surv(y, delta), -lp_train_sp)
    lp_test_sp <- predict(fit_sp, newdata = df_test, type = "lp")
    c_sp_test <- Cindex(Surv(y_test, delta_test), -lp_test_sp)

    s_sp <- tryCatch({
      sf <- survfit(fit_sp, newdata = data.frame(X = 1, Z1 = z1_ref, Z2 = z2_ref))
      idx_t0 <- max(which(sf$time <= t0_S)); sf$surv[idx_t0]
    }, error = function(e) NA)

    cr_sp <- tryCatch({
      s_grid <- sapply(x_grid, function(xv) {
        sf <- survfit(fit_sp, newdata = data.frame(X = xv, Z1 = z1_ref, Z2 = z2_ref))
        idx_t0 <- max(which(sf$time <= t0_CR), 0)
        if (idx_t0 == 0) 1 else sf$surv[idx_t0]
      })
      in_cr <- s_grid >= 0.80
      if (any(in_cr)) list(lo = min(x_grid[in_cr]), hi = max(x_grid[in_cr]))
      else list(lo = NA, hi = NA)
    }, error = function(e) list(lo = NA, hi = NA))

    B_comp <- 500
    boot_xc_sp <- numeric(B_comp)
    for (bb in 1:B_comp) {
      set.seed(bb + 200000)
      df_bb <- df_train[sample(n, replace = TRUE), ]
      fit_bb <- tryCatch(
        coxph(Surv(y, delta) ~ ns(X, df = 4) + Z1 + Z2, data = df_bb),
        error = function(e) NULL)
      if (!is.null(fit_bb)) {
        lp_bb <- predict(fit_bb, newdata = pred_df, type = "lp")
        boot_xc_sp[bb] <- x_grid[which.min(lp_bb)]
      } else {
        boot_xc_sp[bb] <- NA
      }
    }
    ci_xc_sp <- quantile(boot_xc_sp, c(0.025, 0.975), na.rm = TRUE)

    lp_tr_sp <- predict(fit_sp, type = "lp")
    bw_sp <- 1.06 * sd(lp_tr_sp) * length(lp_tr_sp)^(-1/5)
    lp_ref_sp <- predict(fit_sp, newdata = data.frame(X = 1, Z1 = z1_ref, Z2 = z2_ref), type = "lp")
    s_sp_km <- km_S_on_lp(lp_tr_sp, y, delta, lp_ref_sp, t0_S, bw_sp)
    cr_sp_km <- tryCatch({
      s_g <- sapply(x_grid, function(xv) {
        lp_x <- predict(fit_sp, newdata = data.frame(X = xv, Z1 = z1_ref, Z2 = z2_ref), type = "lp")
        km_S_on_lp(lp_tr_sp, y, delta, lp_x, t0_CR, bw_sp)
      })
      ic <- s_g >= 0.80
      if (any(ic)) list(lo = min(x_grid[ic]), hi = max(x_grid[ic])) else list(lo = NA, hi = NA)
    }, error = function(e) list(lo = NA, hi = NA))

    comp_results$spline <- list(
      xc = xc_sp, xc_bias = xc_sp - cp_true,
      xc_ci = as.numeric(ci_xc_sp),
      xc_cov = cp_true > ci_xc_sp[1] & cp_true < ci_xc_sp[2],
      cindex_train = c_sp_train, cindex = c_sp_test,
      s_t0 = s_sp, crlb = cr_sp$lo, crub = cr_sp$hi,
      s_t0_km = s_sp_km, crlb_km = cr_sp_km$lo, crub_km = cr_sp_km$hi)
    cat(sprintf("Spline: Xc=%.2f, C_train=%.4f, C_test=%.4f\n",
                xc_sp, c_sp_train, c_sp_test))
  }, error = function(e) {
    cat(sprintf("Spline FAILED: %s\n", e$message))
    comp_results$spline <<- list(xc=NA, xc_bias=NA, xc_ci=c(NA,NA), xc_cov=NA,
      cindex_train=NA, cindex=NA, s_t0=NA, crlb=NA, crub=NA,
      s_t0_km=NA, crlb_km=NA, crub_km=NA)
  })

  # ---- 6b. Cox + quadratic ----
  tryCatch({
    fit_q <- coxph(Surv(y, delta) ~ X + I(X^2) + Z1 + Z2, data = df_train)
    coefs <- coef(fit_q)
    xc_q <- -coefs["X"] / (2 * coefs["I(X^2)"])

    V <- vcov(fit_q)
    grad <- c(-1/(2*coefs["I(X^2)"]), coefs["X"]/(2*coefs["I(X^2)"]^2), 0, 0)
    xc_q_se <- as.numeric(sqrt(t(grad) %*% V %*% grad))
    ci_xc_q <- c(xc_q - qnorm(0.975) * xc_q_se, xc_q + qnorm(0.975) * xc_q_se)

    lp_train_q <- predict(fit_q, type = "lp")
    c_q_train <- Cindex(Surv(y, delta), -lp_train_q)
    lp_test_q <- predict(fit_q, newdata = df_test, type = "lp")
    c_q_test <- Cindex(Surv(y_test, delta_test), -lp_test_q)

    s_q <- tryCatch({
      sf <- survfit(fit_q, newdata = data.frame(X=1, Z1=z1_ref, Z2=z2_ref))
      idx_t0 <- max(which(sf$time <= t0_S)); sf$surv[idx_t0]
    }, error = function(e) NA)

    cr_q <- tryCatch({
      s_g <- sapply(x_grid, function(xv) {
        sf <- survfit(fit_q, newdata=data.frame(X=xv, Z1=z1_ref, Z2=z2_ref))
        idx <- max(which(sf$time <= t0_CR), 0)
        if (idx == 0) 1 else sf$surv[idx]
      })
      ic <- s_g >= 0.80
      if (any(ic)) list(lo=min(x_grid[ic]), hi=max(x_grid[ic])) else list(lo=NA, hi=NA)
    }, error=function(e) list(lo=NA, hi=NA))

    lp_tr_q <- predict(fit_q, type="lp")
    bw_q <- 1.06 * sd(lp_tr_q) * length(lp_tr_q)^(-1/5)
    lp_ref_q <- predict(fit_q, newdata=data.frame(X=1, Z1=z1_ref, Z2=z2_ref), type="lp")
    s_q_km <- km_S_on_lp(lp_tr_q, y, delta, lp_ref_q, t0_S, bw_q)
    cr_q_km <- tryCatch({
      s_g <- sapply(x_grid, function(xv) {
        lp_x <- predict(fit_q, newdata=data.frame(X=xv, Z1=z1_ref, Z2=z2_ref), type="lp")
        km_S_on_lp(lp_tr_q, y, delta, lp_x, t0_CR, bw_q)
      })
      ic <- s_g >= 0.80
      if (any(ic)) list(lo=min(x_grid[ic]), hi=max(x_grid[ic])) else list(lo=NA, hi=NA)
    }, error=function(e) list(lo=NA, hi=NA))

    comp_results$quadratic <- list(
      xc = as.numeric(xc_q), xc_bias = as.numeric(xc_q - cp_true),
      xc_ci = as.numeric(ci_xc_q),
      xc_cov = cp_true > ci_xc_q[1] & cp_true < ci_xc_q[2],
      cindex_train = c_q_train, cindex = c_q_test,
      s_t0 = s_q, crlb = cr_q$lo, crub = cr_q$hi,
      s_t0_km = s_q_km, crlb_km = cr_q_km$lo, crub_km = cr_q_km$hi)
    cat(sprintf("Quad: Xc=%.2f, C_train=%.4f, C_test=%.4f\n",
                xc_q, c_q_train, c_q_test))
  }, error = function(e) {
    cat(sprintf("Quad FAILED: %s\n", e$message))
    comp_results$quadratic <<- list(xc=NA, xc_bias=NA, xc_ci=c(NA,NA), xc_cov=NA,
      cindex_train=NA, cindex=NA, s_t0=NA, crlb=NA, crub=NA,
      s_t0_km=NA, crlb_km=NA, crub_km=NA)
  })

  # ---- 6c. Oracle piecewise Cox (known change point) ----
  tryCatch({
    B_L <- as.numeric(df_train$X < cp_true)
    B_R <- as.numeric(df_train$X >= cp_true)
    df_ora <- data.frame(y = y, delta = delta,
      X_L = (df_train$X - cp_true) * B_L,
      X_R = (df_train$X - cp_true) * B_R,
      Z1_L = df_train$Z1 * B_L, Z1_R = df_train$Z1 * B_R,
      Z2_L = df_train$Z2 * B_L, Z2_R = df_train$Z2 * B_R)

    fit_ora <- coxph(Surv(y, delta) ~ X_L + X_R + Z1_L + Z1_R + Z2_L + Z2_R, data = df_ora)

    lp_ora_train <- predict(fit_ora, type = "lp")
    c_ora_train <- Cindex(Surv(y, delta), -lp_ora_train)

    B_L_te <- as.numeric(df_test$X < cp_true)
    B_R_te <- as.numeric(df_test$X >= cp_true)
    df_ora_te <- data.frame(
      X_L = (df_test$X - cp_true) * B_L_te, X_R = (df_test$X - cp_true) * B_R_te,
      Z1_L = df_test$Z1 * B_L_te, Z1_R = df_test$Z1 * B_R_te,
      Z2_L = df_test$Z2 * B_L_te, Z2_R = df_test$Z2 * B_R_te)
    lp_ora_test <- predict(fit_ora, newdata = df_ora_te, type = "lp")
    c_ora_test <- Cindex(Surv(y_test, delta_test), -lp_ora_test)

    ref_ora <- data.frame(
      X_L = (1 - cp_true), X_R = 0,
      Z1_L = z1_ref, Z1_R = 0,
      Z2_L = z2_ref, Z2_R = 0)
    s_ora <- tryCatch({
      sf <- survfit(fit_ora, newdata = ref_ora)
      idx <- max(which(sf$time <= t0_S)); sf$surv[idx]
    }, error = function(e) NA)

    cr_ora <- tryCatch({
      s_g <- sapply(x_grid, function(xv) {
        bl <- as.numeric(xv < cp_true); br <- 1 - bl
        nd <- data.frame(
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

    bw_ora <- 1.06 * sd(lp_ora_train) * n^(-1/5)
    lp_ref_ora <- predict(fit_ora, newdata = ref_ora, type = "lp")
    s_ora_km <- km_S_on_lp(lp_ora_train, y, delta, lp_ref_ora, t0_S, bw_ora)
    cr_ora_km <- tryCatch({
      s_g <- sapply(x_grid, function(xv) {
        bl <- as.numeric(xv < cp_true); br <- 1 - bl
        nd <- data.frame(
          X_L = (xv - cp_true) * bl, X_R = (xv - cp_true) * br,
          Z1_L = z1_ref * bl, Z1_R = z1_ref * br,
          Z2_L = z2_ref * bl, Z2_R = z2_ref * br)
        lp_x <- predict(fit_ora, newdata = nd, type = "lp")
        km_S_on_lp(lp_ora_train, y, delta, lp_x, t0_CR, bw_ora)
      })
      ic <- !is.na(s_g) & s_g >= 0.80
      if (any(ic)) list(lo = min(x_grid[ic]), hi = max(x_grid[ic])) else list(lo = NA, hi = NA)
    }, error = function(e) list(lo = NA, hi = NA))

    comp_results$oracle <- list(
      xc = cp_true, xc_bias = 0, xc_ci = c(NA, NA), xc_cov = NA,
      cindex_train = c_ora_train, cindex = c_ora_test,
      s_t0 = s_ora, crlb = cr_ora$lo, crub = cr_ora$hi,
      s_t0_km = s_ora_km, crlb_km = cr_ora_km$lo, crub_km = cr_ora_km$hi)
    cat(sprintf("Oracle: C_train=%.4f, C_test=%.4f\n", c_ora_train, c_ora_test))
  }, error = function(e) {
    cat(sprintf("Oracle FAILED: %s\n", e$message))
    comp_results$oracle <<- list(xc=NA, xc_bias=NA, xc_ci=c(NA,NA), xc_cov=NA,
      cindex_train=NA, cindex=NA, s_t0=NA, crlb=NA, crub=NA,
      s_t0_km=NA, crlb_km=NA, crub_km=NA)
  })

  # ---- 6d. Cox + quadratic with X*Z interactions ----
  tryCatch({
    fit_qi <- coxph(Surv(y, delta) ~ X + I(X^2) + X:Z1 + X:Z2 + I(X^2):Z1 + I(X^2):Z2 + Z1 + Z2,
                    data = df_train)
    co <- coef(fit_qi)
    b_X <- co["X"]; b_X2 <- co["I(X^2)"]
    b_XZ1 <- co["X:Z1"]; b_XZ2 <- co["X:Z2"]
    b_X2Z1 <- co["I(X^2):Z1"]; b_X2Z2 <- co["I(X^2):Z2"]
    N_val <- as.numeric(b_X + b_XZ1 * z1_ref + b_XZ2 * z2_ref)
    D_val <- as.numeric(2 * (b_X2 + b_X2Z1 * z1_ref + b_X2Z2 * z2_ref))
    xc_qi <- -N_val / D_val

    grad_qi <- rep(0, length(co)); names(grad_qi) <- names(co)
    grad_qi["X"] <- -1 / D_val
    grad_qi["I(X^2)"] <- 2 * N_val / D_val^2
    grad_qi["X:Z1"] <- -z1_ref / D_val
    grad_qi["X:Z2"] <- -z2_ref / D_val
    grad_qi["I(X^2):Z1"] <- 2 * N_val * z1_ref / D_val^2
    grad_qi["I(X^2):Z2"] <- 2 * N_val * z2_ref / D_val^2
    V_qi <- vcov(fit_qi)
    xc_qi_se <- as.numeric(sqrt(t(grad_qi) %*% V_qi %*% grad_qi))
    ci_xc_qi <- c(xc_qi - qnorm(0.975) * xc_qi_se, xc_qi + qnorm(0.975) * xc_qi_se)

    lp_tr_qi <- predict(fit_qi, type = "lp")
    c_qi_train <- Cindex(Surv(y, delta), -lp_tr_qi)
    lp_te_qi <- predict(fit_qi, newdata = df_test, type = "lp")
    c_qi_test <- Cindex(Surv(y_test, delta_test), -lp_te_qi)

    s_qi <- tryCatch({
      sf <- survfit(fit_qi, newdata = data.frame(X = 1, X2 = 1, Z1 = z1_ref, Z2 = z2_ref))
      idx <- max(which(sf$time <= t0_S)); sf$surv[idx]
    }, error = function(e) NA)

    cr_qi <- tryCatch({
      s_g <- sapply(x_grid, function(xv) {
        sf <- survfit(fit_qi, newdata = data.frame(X = xv, X2 = xv^2, Z1 = z1_ref, Z2 = z2_ref))
        idx <- max(which(sf$time <= t0_CR), 0)
        if (idx == 0) 1 else sf$surv[idx]
      })
      ic <- s_g >= 0.80
      if (any(ic)) list(lo = min(x_grid[ic]), hi = max(x_grid[ic])) else list(lo = NA, hi = NA)
    }, error = function(e) list(lo = NA, hi = NA))

    bw_qi <- 1.06 * sd(lp_tr_qi) * length(lp_tr_qi)^(-1/5)
    lp_ref_qi <- predict(fit_qi, newdata = data.frame(X = 1, X2 = 1, Z1 = z1_ref, Z2 = z2_ref), type = "lp")
    s_qi_km <- km_S_on_lp(lp_tr_qi, y, delta, lp_ref_qi, t0_S, bw_qi)
    cr_qi_km <- tryCatch({
      s_g <- sapply(x_grid, function(xv) {
        lp_x <- predict(fit_qi, newdata = data.frame(X = xv, X2 = xv^2, Z1 = z1_ref, Z2 = z2_ref), type = "lp")
        km_S_on_lp(lp_tr_qi, y, delta, lp_x, t0_CR, bw_qi)
      })
      ic <- s_g >= 0.80
      if (any(ic)) list(lo = min(x_grid[ic]), hi = max(x_grid[ic])) else list(lo = NA, hi = NA)
    }, error = function(e) list(lo = NA, hi = NA))

    comp_results$quad_interact <- list(
      xc = xc_qi, xc_bias = xc_qi - cp_true,
      xc_ci = as.numeric(ci_xc_qi),
      xc_cov = cp_true > ci_xc_qi[1] & cp_true < ci_xc_qi[2],
      cindex_train = c_qi_train, cindex = c_qi_test,
      s_t0 = s_qi, crlb = cr_qi$lo, crub = cr_qi$hi,
      s_t0_km = s_qi_km, crlb_km = cr_qi_km$lo, crub_km = cr_qi_km$hi)
    cat(sprintf("Quad+int: Xc=%.2f, C_train=%.4f, C_test=%.4f\n",
                xc_qi, c_qi_train, c_qi_test))
  }, error = function(e) {
    cat(sprintf("Quad+int FAILED: %s\n", e$message))
    comp_results$quad_interact <<- list(xc=NA, xc_bias=NA, xc_ci=c(NA,NA), xc_cov=NA,
      cindex_train=NA, cindex=NA, s_t0=NA, crlb=NA, crub=NA,
      s_t0_km=NA, crlb_km=NA, crub_km=NA)
  })

  # ---- 6e. Cox + spline with interactions ----
  tryCatch({
    fit_si <- coxph(Surv(y, delta) ~ ns(X, df = 4) * Z1 + ns(X, df = 4) * Z2,
                    data = df_train)
    lp_at_zref_si <- function(xv) {
      predict(fit_si, newdata = data.frame(X = xv, Z1 = z1_ref, Z2 = z2_ref), type = "lp")
    }
    opt_si <- optimize(lp_at_zref_si, interval = c(x_lo, x_hi))
    xc_si <- opt_si$minimum

    co_si <- coef(fit_si)
    V_si <- vcov(fit_si)
    h_eps <- 1e-5
    grad_si <- numeric(length(co_si))
    for (jj in seq_along(co_si)) {
      co_plus <- co_minus <- co_si
      co_plus[jj] <- co_plus[jj] + h_eps
      co_minus[jj] <- co_minus[jj] - h_eps
      mm_fn <- function(xv, cc) {
        nd <- data.frame(X = xv, Z1 = z1_ref, Z2 = z2_ref)
        mm <- model.matrix(delete.response(terms(fit_si)), data = nd)
        if ("(Intercept)" %in% colnames(mm)) mm <- mm[, -1, drop = FALSE]
        as.numeric(mm %*% cc)
      }
      xc_plus <- optimize(function(xv) mm_fn(xv, co_plus), interval = c(x_lo, x_hi))$minimum
      xc_minus <- optimize(function(xv) mm_fn(xv, co_minus), interval = c(x_lo, x_hi))$minimum
      grad_si[jj] <- (xc_plus - xc_minus) / (2 * h_eps)
    }
    xc_si_se <- as.numeric(sqrt(t(grad_si) %*% V_si %*% grad_si))
    ci_xc_si <- c(xc_si - qnorm(0.975) * xc_si_se, xc_si + qnorm(0.975) * xc_si_se)

    lp_train_si <- predict(fit_si, type = "lp")
    c_si_train <- Cindex(Surv(y, delta), -lp_train_si)
    lp_test_si <- predict(fit_si, newdata = df_test, type = "lp")
    c_si_test <- Cindex(Surv(y_test, delta_test), -lp_test_si)

    s_si <- tryCatch({
      sf <- survfit(fit_si, newdata = data.frame(X = 1, Z1 = z1_ref, Z2 = z2_ref))
      idx_t0 <- max(which(sf$time <= t0_S)); sf$surv[idx_t0]
    }, error = function(e) NA)

    cr_si <- tryCatch({
      s_grid <- sapply(x_grid, function(xv) {
        sf <- survfit(fit_si, newdata = data.frame(X = xv, Z1 = z1_ref, Z2 = z2_ref))
        idx_t0 <- max(which(sf$time <= t0_CR), 0)
        if (idx_t0 == 0) 1 else sf$surv[idx_t0]
      })
      in_cr <- s_grid >= 0.80
      if (any(in_cr)) list(lo = min(x_grid[in_cr]), hi = max(x_grid[in_cr]))
      else list(lo = NA, hi = NA)
    }, error = function(e) list(lo = NA, hi = NA))

    lp_tr_si <- predict(fit_si, type = "lp")
    bw_si <- 1.06 * sd(lp_tr_si) * length(lp_tr_si)^(-1/5)
    lp_ref_si <- predict(fit_si, newdata = data.frame(X = 1, Z1 = z1_ref, Z2 = z2_ref), type = "lp")
    s_si_km <- km_S_on_lp(lp_tr_si, y, delta, lp_ref_si, t0_S, bw_si)
    cr_si_km <- tryCatch({
      s_g <- sapply(x_grid, function(xv) {
        lp_x <- predict(fit_si, newdata = data.frame(X = xv, Z1 = z1_ref, Z2 = z2_ref), type = "lp")
        km_S_on_lp(lp_tr_si, y, delta, lp_x, t0_CR, bw_si)
      })
      ic <- s_g >= 0.80
      if (any(ic)) list(lo = min(x_grid[ic]), hi = max(x_grid[ic])) else list(lo = NA, hi = NA)
    }, error = function(e) list(lo = NA, hi = NA))

    comp_results$spline_interact <- list(
      xc = xc_si, xc_bias = xc_si - cp_true,
      xc_ci = as.numeric(ci_xc_si),
      xc_cov = cp_true > ci_xc_si[1] & cp_true < ci_xc_si[2],
      cindex_train = c_si_train, cindex = c_si_test,
      s_t0 = s_si, crlb = cr_si$lo, crub = cr_si$hi,
      s_t0_km = s_si_km, crlb_km = cr_si_km$lo, crub_km = cr_si_km$hi)
    cat(sprintf("Spline+int: Xc=%.2f, C_train=%.4f, C_test=%.4f\n",
                xc_si, c_si_train, c_si_test))
  }, error = function(e) {
    cat(sprintf("Spline+int FAILED: %s\n", e$message))
    comp_results$spline_interact <<- list(xc=NA, xc_bias=NA, xc_ci=c(NA,NA), xc_cov=NA,
      cindex_train=NA, cindex=NA, s_t0=NA, crlb=NA, crub=NA,
      s_t0_km=NA, crlb_km=NA, crub_km=NA)
  })
}

# === SECTION 7: SAVE RESULTS ===
batch_tag <- Sys.getenv("SLURM_ARRAY_JOB_ID", unset = "")
if (batch_tag == "") batch_tag <- format(Sys.time(), "%Y%m%d%H%M%S")

result <- list(
  # Provenance
  scenario_id = scenario_id, seed = seed, batch_tag = batch_tag,
  timestamp = as.character(Sys.time()),
  n = n, ei = ei, G_type = G_type, beta1 = beta1,
  censor_lo = censor_lo, censor_hi = censor_hi,
  t0_S = t0_S, t0_CR = t0_CR,

  # True values
  theta_true = theta_true, cp_true = cp_true,

  # Observed
  censor_rate = censor_rate,

  # MCE estimates
  est_par = est_par,
  est_cindex = est_cindex,
  c_train = c_train,
  est_cpoint = est_cpoint,
  bias = bias, bias_cp = bias_cp,

  # Bootstrap
  boot_se = boot_se,
  boot_cp_se = boot_cp_se,
  boot_cp_se_direct = boot_cp_se_direct,
  boot_cp_ci_pct = as.numeric(boot_cp_ci_pct),
  boot_cov = boot_cov,
  boot_method = sprintf("smooth_lbfgsb_sh%d", sigma_h),
  nboot_used = nboot_used,
  nboot_run = nboot_run,
  ci_par = ci_par, ci_cp = ci_cp,
  cov_par = cov_par, cov_cp = cov_cp,
  cov_cp_pct = cov_cp_pct,

  # CR estimates
  est_crlb = est_crlb, est_crub = est_crub,

  # S(t0) and C-index (train + test)
  s_train = s_train, s_test = s_test,
  c_test = c_test,

  # Competitors
  competitors = comp_results
)

results_base <- Sys.getenv("RESULTS_DIR", "results")
out_dir <- file.path(results_base, scenario_id)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_file <- file.path(out_dir, sprintf("rep_%04d_j%s.rds", seed, batch_tag))
saveRDS(result, out_file)

cat(sprintf("Saved: %s\nEnd: %s\n", out_file, Sys.time()))
