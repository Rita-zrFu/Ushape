## =========================================================
## 04_plot_fix_time.R
## Plot risk versus BMI at fixed times with bootstrap CI bands
## Standalone script:
## - reads derived analysis data
## - reads fitted model results from results/realdata_fit.rds
## - uses saved bootstrap estimates from 02
## - saves figure to figures/risk_vs_bmi_two_times.pdf
## =========================================================

library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(scales)
library(patchwork)
library(grid)

## -----------------------------
## Paths and settings
## -----------------------------
derived_data_path <- "data/derived/df_for_analysis.csv"
fit_result_path   <- "results/realdata_fit.rds"
fig_dir           <- "figures"

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

SCALE    <- 1
h_kernel <- 1.0

bmi_grid <- seq(18, 35, by = 0.2)
B_use    <- 200

## -----------------------------
## Read data and fit result
## -----------------------------
df <- read.csv(derived_data_path, header = TRUE)
fit_obj <- readRDS(fit_result_path)

theta_hat <- fit_obj$theta_hat
boot_mat  <- fit_obj$boot_mat

if (!"BMIc" %in% names(df)) {
  df$BMIc <- (df$bmi - 25) / SCALE
}

set.seed(1)
boot_idx <- sample.int(nrow(boot_mat), min(B_use, nrow(boot_mat)), replace = FALSE)
boot_par_use <- boot_mat[boot_idx, 1:11, drop = FALSE]

## -----------------------------
## Helper functions copied from model-fitting logic
## -----------------------------
get_par_g <- function(theta, g){
  if (g == 1) return(list(a = 0,        b0 = theta[4],  b1 = theta[8]))
  if (g == 2) return(list(a = theta[1], b0 = theta[5],  b1 = theta[9]))
  if (g == 3) return(list(a = theta[2], b0 = theta[6],  b1 = theta[10]))
  if (g == 4) return(list(a = theta[3], b0 = theta[7],  b1 = theta[11]))
  stop("g must be 1, 2, 3, or 4")
}

H_from_theta <- function(theta, x, z1, z2){
  g <- 1L + z1 * 2L + z2
  H <- numeric(length(x))
  for (gi in 1:4){
    idx <- (g == gi)
    par <- get_par_g(theta, gi)
    H[idx] <- pmax(-x[idx] + par$a, par$b0 + par$b1 * x[idx])
  }
  H
}

cps_from_theta <- function(theta, SCALE = 1){
  sapply(1:4, function(g){
    par <- get_par_g(theta, g)
    25 + (par$a - par$b0) / (1 + par$b1) * SCALE
  })
}

event_times <- function(time, status){
  sort(unique(time[status == 1]))
}

S_from_weights_fast <- function(time, status, w, t_grid){
  ord <- order(time)
  time   <- time[ord]
  status <- status[ord]
  w      <- w[ord]
  
  ut <- sort(unique(time))
  risk_w <- rev(cumsum(rev(w)))
  first_idx <- match(ut, time)
  r_ut <- risk_w[first_idx]
  
  d_map <- tapply(w[status == 1], time[status == 1], sum)
  d_ut <- as.numeric(d_map[match(ut, names(d_map))])
  d_ut[is.na(d_ut)] <- 0
  
  haz <- ifelse(r_ut > 0, d_ut / r_ut, 0)
  S_ut <- cumprod(1 - haz)
  
  idx <- match(t_grid, ut)
  S_grid <- S_ut[idx]
  
  list(time = t_grid, S = S_grid)
}

## -----------------------------
## CP helper
## -----------------------------
cp_df_from_theta <- function(theta_hat, SCALE = 1){
  cp_hat <- cps_from_theta(theta_hat, SCALE = SCALE)
  group_labels <- c("Female, <65", "Female, ≥65", "Male, <65", "Male, ≥65")
  data.frame(Group = group_labels, cp = as.numeric(cp_hat))
}

## -----------------------------
## H(BMI) curve data
## -----------------------------
H_df_from_theta <- function(theta_hat, SCALE = 1, bmi_grid = seq(18, 35, 0.2)){
  group_labels <- c("Female, <65", "Female, ≥65", "Male, <65", "Male, ≥65")
  Hb <- (bmi_grid - 25) / SCALE
  
  bind_rows(lapply(1:4, function(g){
    par <- get_par_g(theta_hat, g)
    Hx  <- pmax(-Hb + par$a, par$b0 + par$b1 * Hb)
    tibble(BMI = bmi_grid, H = Hx, Group = group_labels[g])
  }))
}

## -----------------------------
## Precompute per-group caches
## -----------------------------
make_group_caches <- function(df, theta, SCALE = 1, t_fixed = 15){
  if (!"BMIc" %in% names(df)) {
    df$BMIc <- (df$bmi - 25) / SCALE
  }
  
  G <- list(
    df[df$gender == 0 & df$old == 0, ],  # Female, <65
    df[df$gender == 0 & df$old == 1, ],  # Female, ≥65
    df[df$gender == 1 & df$old == 0, ],  # Male, <65
    df[df$gender == 1 & df$old == 1, ]   # Male, ≥65
  )
  
  group_labels <- c("Female, <65", "Female, ≥65", "Male, <65", "Male, ≥65")
  caches <- vector("list", 4)
  
  for (g in 1:4){
    gdat <- G[[g]]
    if (nrow(gdat) == 0) {
      stop(paste0("Empty group: ", group_labels[g]))
    }
    
    x <- gdat$BMIc
    H_i_raw <- H_from_theta(theta, x, gdat$gender, gdat$old)
    H_i <- -H_i_raw
    
    t_ev <- sort(unique(gdat$time[gdat$d == 1]))
    if (length(t_ev) == 0) {
      stop(paste0("No events in group: ", group_labels[g]))
    }
    
    t_star <- t_ev[which.min(abs(t_ev - t_fixed))]
    par <- get_par_g(theta, g)
    
    caches[[g]] <- list(
      label  = group_labels[g],
      dat    = gdat,
      H_i    = H_i,
      t_ev   = t_ev,
      t_star = t_star,
      par    = par
    )
  }
  
  caches
}

## -----------------------------
## Compute risk curve for one group
## -----------------------------
risk_curve_one_group_fast <- function(cache, theta, SCALE = 1, bmi_grid, h_kernel = 1.0){
  gdat   <- cache$dat
  H_i    <- cache$H_i
  t_ev   <- cache$t_ev
  t_star <- cache$t_star
  par    <- cache$par
  
  Hb <- (bmi_grid - 25) / SCALE
  Hq_raw <- pmax(-Hb + par$a, par$b0 + par$b1 * Hb)
  Hq_tilde <- -Hq_raw
  
  k_gauss <- function(u, h) {
    (1 / (sqrt(2 * pi) * h)) * exp(-(u * u) / (2 * h * h))
  }
  
  S_star <- numeric(length(bmi_grid))
  for (j in seq_along(bmi_grid)){
    w  <- k_gauss(Hq_tilde[j] - H_i, h_kernel)
    Sj <- S_from_weights_fast(gdat$time, gdat$d, w, t_ev)
    idx <- match(t_star, Sj$time)
    S_star[j] <- if (!is.na(idx)) Sj$S[idx] else tail(Sj$S, 1)
  }
  
  ok <- is.finite(Hq_tilde) & is.finite(S_star)
  if (sum(ok) < 2L) stop("Not enough finite points in risk curve.")
  
  H_ok_tilde <- Hq_tilde[ok]
  S_ok       <- S_star[ok]
  BMI_ok     <- bmi_grid[ok]
  
  ord <- order(H_ok_tilde)
  fit <- isoreg(H_ok_tilde[ord], S_ok[ord])
  H_mono <- H_ok_tilde[ord]
  S_mono <- fit$yf
  
  ordB <- order(BMI_ok)
  BMI_plot    <- BMI_ok[ordB]
  Htilde_on_B <- H_ok_tilde[ordB]
  S_on_B      <- splinefun(H_mono, S_mono, method = "hyman")(Htilde_on_B)
  
  tibble(
    BMI = BMI_plot,
    Risk = 1 - S_on_B,
    t_star = t_star
  )
}

## -----------------------------
## Point estimate curves
## -----------------------------
risk_point_df <- function(df, theta_hat, SCALE = 1, t_fixed = 15, bmi_grid, h_kernel = 1.0){
  caches_hat <- make_group_caches(df, theta_hat, SCALE, t_fixed)
  bind_rows(lapply(1:4, function(g){
    out <- risk_curve_one_group_fast(caches_hat[[g]], theta_hat, SCALE, bmi_grid, h_kernel)
    out$Group <- caches_hat[[g]]$label
    out
  }))
}

## -----------------------------
## Bootstrap CI band
## -----------------------------
risk_ci_df <- function(df, theta_hat, boot_par_use,
                       SCALE = 1, t_fixed = 15, bmi_grid, h_kernel = 1.0){
  
  base_caches <- make_group_caches(df, theta_hat, SCALE, t_fixed)
  
  Gdat   <- lapply(base_caches, `[[`, "dat")
  Tev    <- lapply(base_caches, `[[`, "t_ev")
  Tstar  <- sapply(base_caches, `[[`, "t_star")
  labels <- sapply(base_caches, `[[`, "label")
  
  make_cache_for_theta <- function(theta_b){
    caches_b <- vector("list", 4)
    for (g in 1:4){
      gdat <- Gdat[[g]]
      x <- if (!"BMIc" %in% names(gdat)) (gdat$bmi - 25) / SCALE else gdat$BMIc
      
      H_i <- -H_from_theta(theta_b, x, gdat$gender, gdat$old)
      
      caches_b[[g]] <- list(
        label  = labels[g],
        dat    = gdat,
        H_i    = H_i,
        t_ev   = Tev[[g]],
        t_star = Tstar[g],
        par    = get_par_g(theta_b, g)
      )
    }
    caches_b
  }
  
  boot_long <- bind_rows(lapply(seq_len(nrow(boot_par_use)), function(b){
    th <- as.numeric(boot_par_use[b, ])
    caches_b <- make_cache_for_theta(th)
    
    bind_rows(lapply(1:4, function(g){
      out <- risk_curve_one_group_fast(caches_b[[g]], th, SCALE, bmi_grid, h_kernel)
      out$Group <- caches_b[[g]]$label
      out$b <- b
      out
    }))
  }))
  
  boot_long %>%
    group_by(Group, BMI) %>%
    summarise(
      lwr = quantile(Risk, 0.025, na.rm = TRUE),
      upr = quantile(Risk, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

## -----------------------------
## Plot helpers
## -----------------------------
plot_risk_panel <- function(df, theta_hat, boot_par_use,
                            SCALE = 1, t_fixed = 15, bmi_grid, h_kernel = 1.0,
                            y_max = 0.45, show_legend = TRUE, show_y = TRUE,
                            col_title = NULL){
  
  col_map <- c(
    "Female, <65" = "#D55E00",
    "Female, ≥65" = "#0072B2",
    "Male, <65"   = "#009E73",
    "Male, ≥65"   = "#CC79A7"
  )
  
  df_point <- risk_point_df(df, theta_hat, SCALE, t_fixed, bmi_grid, h_kernel)
  df_ci    <- risk_ci_df(df, theta_hat, boot_par_use, SCALE, t_fixed, bmi_grid, h_kernel)
  
  ggplot() +
    geom_ribbon(
      data = df_ci,
      aes(x = BMI, ymin = lwr, ymax = upr, fill = Group),
      alpha = 0.12,
      color = NA,
      show.legend = FALSE
    ) +
    geom_line(
      data = df_point,
      aes(x = BMI, y = Risk, color = Group, linetype = Group),
      linewidth = 1.2,
      show.legend = TRUE
    ) +
    scale_color_manual(values = col_map, name = "Group") +
    scale_fill_manual(values = col_map, guide = "none") +
    scale_linetype_manual(
      values = c(
        "Female, <65" = "solid",
        "Female, ≥65" = "solid",
        "Male, <65"   = "dashed",
        "Male, ≥65"   = "dashed"
      ),
      name = "Group"
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      labels = percent_format(accuracy = 1)
    ) +
    labs(
      title = col_title,
      x = "BMI",
      y = "Estimated event risk"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = if (show_legend) "top" else "none",
      legend.justification = "left",
      legend.key.width = unit(1.6, "cm"),
      legend.key.height = unit(0.6, "cm"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    {
      if (!show_y) {
        theme(
          axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank()
        )
      } else {
        theme()
      }
    }
}

plot_H_panel <- function(theta_hat, SCALE = 1, bmi_grid, show_y = TRUE){
  col_map <- c(
    "Female, <65" = "#D55E00",
    "Female, ≥65" = "#0072B2",
    "Male, <65"   = "#009E73",
    "Male, ≥65"   = "#CC79A7"
  )
  
  df_H <- H_df_from_theta(theta_hat, SCALE, bmi_grid)
  group_levels <- c("Female, <65", "Female, ≥65", "Male, <65", "Male, ≥65")
  df_H$Group <- factor(df_H$Group, levels = group_levels)
  
  p <- ggplot(df_H, aes(x = BMI, y = H, color = Group, linetype = Group)) +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = col_map, breaks = group_levels, guide = "none") +
    scale_linetype_manual(
      values = c(
        "Female, <65" = "solid",
        "Female, ≥65" = "solid",
        "Male, <65"   = "dashed",
        "Male, ≥65"   = "dashed"
      ),
      breaks = group_levels,
      guide = "none"
    ) +
    labs(
      x = "BMI",
      y = if (show_y) "Estimated risk score" else NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
  
  if (!show_y){
    p <- p + theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  
  p
}

## -----------------------------
## Build final four-panel figure
## -----------------------------
p_risk15 <- plot_risk_panel(
  df = df,
  theta_hat = theta_hat,
  boot_par_use = boot_par_use,
  SCALE = SCALE,
  t_fixed = 15,
  bmi_grid = bmi_grid,
  h_kernel = h_kernel,
  y_max = 0.45,
  show_legend = TRUE,
  show_y = TRUE,
  col_title = "15-year event risk"
)

p_H15 <- plot_H_panel(
  theta_hat = theta_hat,
  SCALE = SCALE,
  bmi_grid = bmi_grid,
  show_y = TRUE
)

p_risk10 <- plot_risk_panel(
  df = df,
  theta_hat = theta_hat,
  boot_par_use = boot_par_use,
  SCALE = SCALE,
  t_fixed = 10,
  bmi_grid = bmi_grid,
  h_kernel = h_kernel,
  y_max = 0.45,
  show_legend = TRUE,
  show_y = FALSE,
  col_title = "10-year event risk"
)

p_H10 <- plot_H_panel(
  theta_hat = theta_hat,
  SCALE = SCALE,
  bmi_grid = bmi_grid,
  show_y = FALSE
)

p_all <- (p_risk15 | p_risk10) / (p_H15 | p_H10) +
  plot_layout(guides = "collect", heights = c(2.0, 1.3)) +
  plot_annotation() &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.justification = "center"
  )

print(p_all)

ggsave(
  filename = file.path(fig_dir, "risk_vs_bmi_two_times.pdf"),
  plot = p_all,
  width = 12,
  height = 7
)