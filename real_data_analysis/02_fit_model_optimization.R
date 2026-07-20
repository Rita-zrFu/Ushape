## =========================
## 02_fit_model_optimization.R
## JITTER (COBYLA) + OPTIONAL UNBOUNDED REFINE (accept only if C-index improves)
##
## Workflow:
## 1) Read derived analysis dataset
## 2) Build loess-pair-based initial values
## 3) Run jittered bounded COBYLA optimization and keep the best solution
## 4) Refine the best point estimate using unbounded Nelder-Mead,
##    but accept only if C-index improves
## 5) Run bootstrap with the same strategy
## 6) Save fitted object and summary tables
## =========================

if (!requireNamespace("mycpp", quietly = TRUE)) {
  stop("Package 'mycpp' is required. Please install it first, e.g. remotes::install_local('mycpp').")
}

library(mycpp)
library(nloptr)
library(future)
library(future.apply)

## =========================
## Paths and settings
## =========================
derived_data_dir <- "data/derived"
results_dir <- "results"

input_file <- file.path(derived_data_dir, "df_for_analysis.csv")
fit_file   <- file.path(results_dir, "realdata_fit.rds")

if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

n_workers <- max(1, parallel::detectCores() - 1)
SCALE <- 1

## manuscript-level settings
J_point  <- 50
J_boot   <- 10
B        <- 500
sd_bmi   <- 0.20
frac_band <- 0.25
cp_lo <- 18
cp_hi <- 35

## =========================
## Read data
## =========================
df <- read.csv(input_file, header = TRUE)

required_cols <- c("bmi", "BMIc", "gender", "old", "time", "d")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in input data: ", paste(missing_cols, collapse = ", "))
}

## =========================
## 0) utilities
## =========================
cpoint_BMI <- function(a, b0, b1, SCALE = 1) {
  xstar_c <- (a - b0) / (1 + b1)
  25 + SCALE * xstar_c
}

get_par_g <- function(theta, g) {
  if (g == 1) return(list(a = 0,        b0 = theta[4],  b1 = theta[8]))
  if (g == 2) return(list(a = theta[1], b0 = theta[5],  b1 = theta[9]))
  if (g == 3) return(list(a = theta[2], b0 = theta[6],  b1 = theta[10]))
  if (g == 4) return(list(a = theta[3], b0 = theta[7],  b1 = theta[11]))
  stop("g must be 1, 2, 3, or 4")
}

H_from_theta <- function(theta, x, z1, z2) {
  g <- 1L + z1 * 2L + z2
  H <- numeric(length(x))
  for (gi in 1:4) {
    idx <- (g == gi)
    par <- get_par_g(theta, gi)
    H[idx] <- pmax(-x[idx] + par$a, par$b0 + par$b1 * x[idx])
  }
  H
}

## objective: maximize C-index => minimize negative C-index
obj_fun <- function(theta, x, z1, z2, y, d) {
  H <- H_from_theta(theta, x, z1, z2)
  -myCindex(H, y, d)
}

## bounds
lb <- c(rep(-20, 7), rep(0, 4))
ub <- c(rep(40, 7),  rep(4, 4))

proj_to_bounds <- function(x, lb, ub, eps = 1e-8) {
  pmin(pmax(x, lb + eps), ub - eps)
}

fit_once_with_init <- function(dat, theta_init) {
  y  <- dat$time
  d  <- dat$d
  x  <- dat$BMIc
  z1 <- dat$gender
  z2 <- dat$old
  
  theta0 <- proj_to_bounds(as.numeric(theta_init), lb, ub)
  
  fit <- tryCatch(
    nloptr(
      x0     = theta0,
      eval_f = function(th) obj_fun(th, x, z1, z2, y, d),
      lb     = lb,
      ub     = ub,
      opts   = list(
        algorithm   = "NLOPT_LN_COBYLA",
        maxeval     = 2000,
        xtol_rel    = 1e-6,
        print_level = 0
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit) || is.null(fit$solution)) {
    return(list(
      theta = rep(NA_real_, length(theta0)),
      cindex = NA_real_,
      status = NA_integer_
    ))
  }
  
  th  <- fit$solution
  H   <- H_from_theta(th, x, z1, z2)
  cdx <- myCindex(H, y, d)
  
  list(theta = th, cindex = cdx, status = fit$status)
}

.approx_y2x <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 2L) {
    mx <- if (length(x) == 0L) NA_real_ else mean(x)
    return(function(h) rep(mx, length(h)))
  }
  
  ord <- order(y, x)
  x <- x[ord]
  y <- y[ord]
  
  keep <- !duplicated(y)
  x <- x[keep]
  y <- y[keep]
  
  if (length(x) < 2L) {
    mx <- mean(x)
    return(function(h) rep(mx, length(h)))
  }
  
  approxfun(y, x, rule = 2)
}

cps_from_theta <- function(theta, SCALE = 1) {
  sapply(1:4, function(g) {
    par <- get_par_g(theta, g)
    25 + (par$a - par$b0) / (1 + par$b1) * SCALE
  })
}

## =========================
## 1) loess pairs initial
## =========================
.build_binned <- function(df, bin_breaks = seq(15, 40, by = 2.5)) {
  stopifnot(all(c("bmi", "time", "d", "old", "gender") %in% names(df)))
  
  mids <- head(bin_breaks, -1) + diff(bin_breaks) / 2
  
  df$old_f    <- ifelse(df$old == 1L, "≥65", "<65")
  df$gender_f <- ifelse(df$gender == 1L, "Male", "Female")
  df$bmi_bin  <- cut(df$bmi, breaks = bin_breaks, include.lowest = TRUE, right = FALSE)
  df$bmi_mid  <- mids[as.integer(df$bmi_bin)]
  
  key <- interaction(df$old_f, df$gender_f, df$bmi_bin, drop = TRUE)
  N   <- as.vector(tapply(df$d,    key, length))
  Ev  <- as.vector(tapply(df$d,    key, sum, na.rm = TRUE))
  PY  <- as.vector(tapply(df$time, key, sum, na.rm = TRUE))
  
  idx <- match(levels(key), key)
  kk  <- df[idx, c("old_f", "gender_f", "bmi_bin", "bmi_mid")]
  
  out <- data.frame(kk, N = N, Events = Ev, PY = PY, IR_1000PY = Ev / PY * 1000)
  out[is.finite(out$bmi_mid) & out$PY > 0, ]
}

.group_pairs_and_inits <- function(df_group,
                                   SCALE = 1,
                                   span_loess = 0.6,
                                   n_levels   = 12,
                                   frac_band  = 0.40) {
  xBMI <- df_group$bmi_mid
  y    <- df_group$IR_1000PY
  
  lo  <- loess(y ~ xBMI, span = span_loess, degree = 1)
  y_s <- predict(lo, xBMI)
  i0  <- which.min(y_s)
  
  xL <- xBMI[seq_len(i0)]
  yL <- y_s[seq_len(i0)]
  xR <- xBMI[i0:length(xBMI)]
  yR <- y_s[i0:length(y_s)]
  
  y0   <- y_s[i0]
  y_hi <- y0 + frac_band * (max(y_s, na.rm = TRUE) - y0)
  hs   <- seq(y0, y_hi, length.out = n_levels)
  
  xl_BMI <- .approx_y2x(xL, yL)(hs)
  xr_BMI <- .approx_y2x(xR, yR)(hs)
  
  xl_c <- (xl_BMI - 25) / SCALE
  xr_c <- (xr_BMI - 25) / SCALE
  fit  <- lm(xl_c ~ xr_c)
  
  beta1    <- -unname(coef(fit)[2])
  x_star_c <-  unname(coef(fit)[1]) / (1 + beta1)
  alpha    <- x_star_c
  beta0    <- -beta1 * x_star_c
  
  cp <- cpoint_BMI(alpha, beta0, beta1, SCALE)
  
  list(
    x_pairs = data.frame(xl = xl_BMI, xr = xr_BMI, y = hs),
    alpha = alpha,
    beta0 = beta0,
    beta1 = beta1,
    cp = cp
  )
}

build_base_pairs <- function(df,
                             SCALE = 1,
                             span_loess = 0.6,
                             n_levels = 12,
                             frac_band = 0.40,
                             bin_breaks = seq(15, 40, by = 2.5)) {
  db <- .build_binned(df, bin_breaks)
  
  pick <- function(of, gf) {
    gdf <- db[db$old_f == of & db$gender_f == gf, ]
    .group_pairs_and_inits(gdf, SCALE, span_loess, n_levels, frac_band)
  }
  
  g1 <- pick("<65", "Female")
  g2 <- pick("≥65", "Female")
  g3 <- pick("<65", "Male")
  g4 <- pick("≥65", "Male")
  
  list(g1 = g1, g2 = g2, g3 = g3, g4 = g4)
}

inits_from_groups <- function(G) {
  theta <- c(
    a01   = G$g2$alpha,
    a10   = G$g3$alpha,
    a11   = G$g4$alpha,
    b0_00 = G$g1$beta0,
    b0_01 = G$g2$beta0,
    b0_10 = G$g3$beta0,
    b0_11 = G$g4$beta0,
    b1_00 = G$g1$beta1,
    b1_01 = G$g2$beta1,
    b1_10 = G$g3$beta1,
    b1_11 = G$g4$beta1
  )
  
  cps <- c(g1 = G$g1$cp, g2 = G$g2$cp, g3 = G$g3$cp, g4 = G$g4$cp)
  
  list(theta_init = as.numeric(theta), cp = cps)
}

## =========================
## 2) jitter pairs
## =========================
jitter_pairs_once <- function(x_pairs, sd_bmi = 0.15, minBMI = 15, maxBMI = 40) {
  n <- nrow(x_pairs)
  xp <- x_pairs
  xp$xl <- pmin(pmax(xp$xl + rnorm(n, 0, sd_bmi), minBMI), maxBMI)
  xp$xr <- pmin(pmax(xp$xr + rnorm(n, 0, sd_bmi), minBMI), maxBMI)
  xp
}

group_init_from_pairs <- function(x_pairs, SCALE = 1) {
  xl_c <- (x_pairs$xl - 25) / SCALE
  xr_c <- (x_pairs$xr - 25) / SCALE
  fit  <- lm(xl_c ~ xr_c)
  
  beta1    <- -unname(coef(fit)[2])
  x_star_c <-  unname(coef(fit)[1]) / (1 + beta1)
  alpha    <- x_star_c
  beta0    <- -beta1 * x_star_c
  cp       <- cpoint_BMI(alpha, beta0, beta1, SCALE)
  
  list(alpha = alpha, beta0 = beta0, beta1 = beta1, cp = cp)
}

cp_in_range <- function(cps, lo, hi) {
  all(is.finite(cps)) && all(cps >= lo & cps <= hi)
}

one_jitter_try <- function(BASE, df, SCALE, sd_bmi, cp_lo, cp_hi,
                           span_loess, n_levels, frac_band) {
  JP <- list(
    g1 = jitter_pairs_once(BASE$g1$x_pairs, sd_bmi),
    g2 = jitter_pairs_once(BASE$g2$x_pairs, sd_bmi),
    g3 = jitter_pairs_once(BASE$g3$x_pairs, sd_bmi),
    g4 = jitter_pairs_once(BASE$g4$x_pairs, sd_bmi)
  )
  
  G <- list(
    g1 = group_init_from_pairs(JP$g1, SCALE),
    g2 = group_init_from_pairs(JP$g2, SCALE),
    g3 = group_init_from_pairs(JP$g3, SCALE),
    g4 = group_init_from_pairs(JP$g4, SCALE)
  )
  
  init_cp <- c(G$g1$cp, G$g2$cp, G$g3$cp, G$g4$cp)
  if (!cp_in_range(init_cp, cp_lo, cp_hi)) return(NULL)
  
  init <- inits_from_groups(G)$theta_init
  fit  <- fit_once_with_init(df, init)
  if (!is.finite(fit$cindex)) return(NULL)
  
  final_cp <- cps_from_theta(fit$theta, SCALE)
  
  list(
    theta = fit$theta,
    cindex = fit$cindex,
    init_cp = init_cp,
    final_cp = final_cp,
    init = init
  )
}

run_jitter_cobyla_parallel <- function(df, J = 50,
                                       sd_bmi = 0.15,
                                       cp_lo = 18, cp_hi = 35,
                                       SCALE = 1,
                                       span_loess = 0.6,
                                       n_levels = 12,
                                       frac_band = 0.40,
                                       bin_breaks = seq(15, 40, by = 2.5),
                                       workers = 4) {
  if (!"BMIc" %in% names(df)) {
    df$BMIc <- (df$bmi - 25) / SCALE
  }
  
  BASE <- build_base_pairs(df, SCALE, span_loess, n_levels, frac_band, bin_breaks)
  
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(multisession, workers = workers)
  
  tries <- future.apply::future_lapply(
    seq_len(J),
    function(j) {
      library(nloptr)
      library(mycpp)
      one_jitter_try(BASE, df, SCALE, sd_bmi, cp_lo, cp_hi,
                     span_loess, n_levels, frac_band)
    },
    future.seed = TRUE
  )
  
  kept <- Filter(Negate(is.null), tries)
  
  if (length(kept) == 0L) {
    return(list(
      best = NULL,
      kept = list(),
      base_pairs = BASE,
      settings = list(
        J = J,
        sd_bmi = sd_bmi,
        cp_range = c(cp_lo, cp_hi),
        span_loess = span_loess,
        n_levels = n_levels,
        frac_band = frac_band
      )
    ))
  }
  
  idx_best <- which.max(sapply(kept, `[[`, "cindex"))
  best <- kept[[idx_best]]
  
  list(
    best = best,
    kept = kept,
    base_pairs = BASE,
    settings = list(
      J = J,
      sd_bmi = sd_bmi,
      cp_range = c(cp_lo, cp_hi),
      span_loess = span_loess,
      n_levels = n_levels,
      frac_band = frac_band
    )
  )
}

## =========================
## 3) unbounded refine
## =========================
refine_unbounded_accept_if_improve <- function(df, theta_start,
                                               SCALE = 1,
                                               maxit = 1500,
                                               reltol = 1e-8) {
  if (!"BMIc" %in% names(df)) {
    df$BMIc <- (df$bmi - 25) / SCALE
  }
  
  x  <- df$BMIc
  z1 <- df$gender
  z2 <- df$old
  y  <- df$time
  d  <- df$d
  
  H0 <- H_from_theta(theta_start, x, z1, z2)
  c0 <- myCindex(H0, y, d)
  
  fwrap <- function(th) obj_fun(th, x, z1, z2, y, d)
  
  opt <- optim(
    par     = as.numeric(theta_start),
    fn      = fwrap,
    method  = "Nelder-Mead",
    control = list(maxit = maxit, reltol = reltol)
  )
  
  th1 <- opt$par
  H1  <- H_from_theta(th1, x, z1, z2)
  c1  <- myCindex(H1, y, d)
  
  if (is.finite(c1) && is.finite(c0) && c1 > c0) {
    list(
      theta = th1,
      cindex = c1,
      improved = TRUE,
      conv = opt$convergence,
      feval = opt$counts["function"]
    )
  } else {
    list(
      theta = theta_start,
      cindex = c0,
      improved = FALSE,
      conv = opt$convergence,
      feval = opt$counts["function"]
    )
  }
}

## =========================
## 4) point estimate
## =========================
set.seed(456)

res_jit <- run_jitter_cobyla_parallel(
  df = df,
  J = J_point,
  sd_bmi = sd_bmi,
  frac_band = frac_band,
  cp_lo = cp_lo,
  cp_hi = cp_hi,
  SCALE = SCALE,
  workers = n_workers
)

if (is.null(res_jit$best)) {
  stop("Point-estimate optimization failed: no valid jittered solution was found.")
}

theta_hat0 <- res_jit$best$theta
cindex0    <- res_jit$best$cindex

ref_pt <- refine_unbounded_accept_if_improve(
  df = df,
  theta_start = theta_hat0,
  SCALE = SCALE,
  maxit = 1500
)

theta_hat  <- ref_pt$theta
cindex_hat <- ref_pt$cindex

cat(sprintf(
  "Point estimate refine improved? %s | C-index: %.6f -> %.6f\n",
  ref_pt$improved, cindex0, cindex_hat
))

## =========================
## 5) bootstrap
## =========================
set.seed(456)
ind_big <- replicate(B, sample.int(nrow(df), replace = TRUE), simplify = FALSE)

boot_estimate_one <- function(df_boot,
                              J_boot = 10,
                              sd_bmi = 0.20,
                              cp_lo = 18, cp_hi = 35,
                              SCALE = 1,
                              span_loess = 0.6,
                              n_levels = 12,
                              frac_band = 0.25,
                              bin_breaks = seq(15, 40, by = 2.5),
                              refine_maxit = 1500) {
  if (!"BMIc" %in% names(df_boot)) {
    df_boot$BMIc <- (df_boot$bmi - 25) / SCALE
  }
  
  BASE <- build_base_pairs(df_boot, SCALE, span_loess, n_levels, frac_band, bin_breaks)
  
  kept <- vector("list", J_boot)
  k <- 0L
  
  for (j in seq_len(J_boot)) {
    out <- one_jitter_try(BASE, df_boot, SCALE, sd_bmi, cp_lo, cp_hi,
                          span_loess, n_levels, frac_band)
    if (!is.null(out) && is.finite(out$cindex)) {
      k <- k + 1L
      kept[[k]] <- out
    }
  }
  
  kept <- kept[seq_len(k)]
  
  if (length(kept) == 0L) {
    return(list(
      theta = rep(NA_real_, 11),
      cindex = NA_real_,
      refined = FALSE,
      improved = FALSE
    ))
  }
  
  idx_best <- which.max(sapply(kept, `[[`, "cindex"))
  best0 <- kept[[idx_best]]
  
  th0 <- best0$theta
  
  ref <- refine_unbounded_accept_if_improve(
    df = df_boot,
    theta_start = th0,
    SCALE = SCALE,
    maxit = refine_maxit
  )
  
  list(
    theta = ref$theta,
    cindex = ref$cindex,
    refined = TRUE,
    improved = ref$improved
  )
}

old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)
future::plan(multisession, workers = n_workers)

t1 <- proc.time()

set.seed(23)
boot_list <- future_lapply(
  ind_big,
  function(ind) {
    library(mycpp)
    library(nloptr)
    
    dfb <- df[ind, , drop = FALSE]
    
    res <- boot_estimate_one(
      df_boot = dfb,
      J_boot = J_boot,
      sd_bmi = sd_bmi,
      frac_band = frac_band,
      cp_lo = cp_lo,
      cp_hi = cp_hi,
      SCALE = SCALE,
      refine_maxit = 1500
    )
    
    c(res$theta, res$cindex)
  },
  future.seed = TRUE
)

t2 <- proc.time()
print(t2 - t1)

boot_mat <- do.call(rbind, boot_list)
colnames(boot_mat) <- c(paste0("param_", 1:11), "cindex")

saveRDS(
  list(
    theta_hat = theta_hat,
    boot_mat = boot_mat
  ),
  fit_file
)

## =========================
## 6) summaries
## =========================
param_names <- c(
  "a01", "a10", "a11",
  "b0_00", "b0_01", "b0_10", "b0_11",
  "b1_00", "b1_01", "b1_10", "b1_11"
)

names(theta_hat) <- param_names

boot_par <- boot_mat[, 1:11, drop = FALSE]
colnames(boot_par) <- param_names

se_par <- apply(boot_par, 2, sd, na.rm = TRUE)
ci_perc <- t(apply(boot_par, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
colnames(ci_perc) <- c("perc_lwr", "perc_upr")

ci_norm <- cbind(
  norm_lwr = theta_hat - 1.96 * se_par,
  norm_upr = theta_hat + 1.96 * se_par
)

summary_par <- data.frame(
  param     = param_names,
  est       = as.numeric(theta_hat),
  se        = se_par,
  ci_perc,
  ci_norm,
  row.names = NULL
)

print(summary_par)

cp_hat <- cps_from_theta(theta_hat, SCALE = SCALE)
names(cp_hat) <- c("cp_g1_<65_F", "cp_g2_≥65_F", "cp_g3_<65_M", "cp_g4_≥65_M")
print(cp_hat)

cp_boot <- t(apply(boot_par, 1, cps_from_theta, SCALE = SCALE))
colnames(cp_boot) <- names(cp_hat)

cp_ci_perc <- t(apply(cp_boot, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
colnames(cp_ci_perc) <- c("perc_lwr", "perc_upr")

cp_se <- apply(cp_boot, 2, sd, na.rm = TRUE)
cp_ci_norm <- cbind(
  norm_lwr = cp_hat - 1.96 * cp_se,
  norm_upr = cp_hat + 1.96 * cp_se
)

summary_cp <- data.frame(
  group     = names(cp_hat),
  est       = as.numeric(cp_hat),
  se        = cp_se,
  cp_ci_perc,
  cp_ci_norm,
  row.names = NULL
)

print(summary_cp)

cindex_boot <- boot_mat[, "cindex"]
cindex_se   <- sd(cindex_boot, na.rm = TRUE)
cindex_ci   <- quantile(cindex_boot, probs = c(0.025, 0.975), na.rm = TRUE)

cindex_summary <- data.frame(
  est = cindex_hat,
  se = cindex_se,
  perc_lwr = as.numeric(cindex_ci[1]),
  perc_upr = as.numeric(cindex_ci[2]),
  norm_lwr = cindex_hat - 1.96 * cindex_se,
  norm_upr = cindex_hat + 1.96 * cindex_se
)

print(cindex_summary)

## optional: save summary tables
write.csv(summary_par, file.path(results_dir, "summary_parameters.csv"), row.names = FALSE)
write.csv(summary_cp, file.path(results_dir, "summary_change_points.csv"), row.names = FALSE)
write.csv(cindex_summary, file.path(results_dir, "summary_cindex.csv"), row.names = FALSE)

cat("Saved fitted object to: ", fit_file, "\n", sep = "")
cat("Saved summary tables to: ", results_dir, "\n", sep = "")