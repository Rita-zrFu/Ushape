## =========================
## 03_plot_fix_bmi.R
## Plot risk versus time for fixed BMI values
## Standalone script:
## - reads derived analysis data
## - reads fitted model results from results/realdata_fit.rds
## - saves figure to figures/risk_vs_time_fixed_bmi.png
## =========================

library(dplyr)
library(tibble)
library(ggplot2)
library(scales)
library(grid)

## -----------------------------
## Paths and settings
## -----------------------------
derived_data_path <- "data/derived/df_for_analysis.csv"
fit_result_path   <- "results/realdata_fit.rds"
fig_dir           <- "figures"

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

SCALE    <- 1
h_kernel <- 0.7
bmi_vec  <- c(22, 27)

## -----------------------------
## Read data and fit result
## -----------------------------
df <- read.csv(derived_data_path, header = TRUE)
fit_obj <- readRDS(fit_result_path)

theta_hat <- fit_obj$theta_hat

if (!"BMIc" %in% names(df)) {
  df$BMIc <- (df$bmi - 25) / SCALE
}

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
## Risk curve over time for one subgroup at one fixed BMI
## -----------------------------
risk_curve_over_time_one_group <- function(gdat, theta, b0, SCALE = 1, h_kernel = 0.7){
  x <- if (!"BMIc" %in% names(gdat)) (gdat$bmi - 25) / SCALE else gdat$BMIc
  H_i <- H_from_theta(theta, x, gdat$gender, gdat$old)
  
  g_idx <- 1L + unique(gdat$gender) * 2L + unique(gdat$old)
  par   <- get_par_g(theta, g_idx)
  
  Hb0 <- (b0 - 25) / SCALE
  H_query <- pmax(-Hb0 + par$a, par$b0 + par$b1 * Hb0)
  
  k_gauss <- function(u, h) {
    (1 / (sqrt(2 * pi) * h)) * exp(-(u * u) / (2 * h * h))
  }
  
  w <- k_gauss(H_query - H_i, h_kernel)
  
  t_ev <- event_times(gdat$time, gdat$d)
  Sj <- S_from_weights_fast(gdat$time, gdat$d, w, t_ev)
  
  tibble(
    time = Sj$time,
    risk = 1 - Sj$S
  )
}

## -----------------------------
## Build multi-group / multi-BMI data
## -----------------------------
risk_vs_time_multi_bmi <- function(df, theta_hat,
                                   bmi_vec = c(22, 27, 32),
                                   SCALE = 1, h_kernel = 0.7){
  if (!"BMIc" %in% names(df)) {
    df$BMIc <- (df$bmi - 25) / SCALE
  }
  
  G <- list(
    df[df$gender == 0 & df$old == 0, ],  # Female, <65
    df[df$gender == 0 & df$old == 1, ],  # Female, ≥65
    df[df$gender == 1 & df$old == 0, ],  # Male, <65
    df[df$gender == 1 & df$old == 1, ]   # Male, ≥65
  )
  
  labs <- c("Female, <65", "Female, ≥65", "Male, <65", "Male, ≥65")
  
  bind_rows(lapply(bmi_vec, function(b0){
    bind_rows(lapply(1:4, function(g){
      risk_curve_over_time_one_group(
        gdat = G[[g]],
        theta = theta_hat,
        b0 = b0,
        SCALE = SCALE,
        h_kernel = h_kernel
      ) %>%
        mutate(
          BMI = b0,
          Group = labs[g]
        )
    }))
  }))
}

## -----------------------------
## Plotting
## -----------------------------
plot_risk_vs_time_facet_bmi <- function(df_rt){
  col_map <- c(
    "Female, <65" = "#D55E00",
    "Female, ≥65" = "#0072B2",
    "Male, <65"   = "#009E73",
    "Male, ≥65"   = "#CC79A7"
  )
  
  lt_map <- c(
    "Female, <65" = "solid",
    "Female, ≥65" = "solid",
    "Male, <65"   = "dashed",
    "Male, ≥65"   = "dashed"
  )
  
  ggplot(df_rt, aes(x = time, y = risk, color = Group, linetype = Group)) +
    geom_line(linewidth = 1.15) +
    scale_color_manual(values = col_map, name = "Group") +
    scale_linetype_manual(values = lt_map, guide = "none") +
    scale_y_continuous(
      limits = c(0, NA),
      labels = scales::percent_format(accuracy = 1)
    ) +
    facet_wrap(
      ~BMI,
      nrow = 1,
      labeller = labeller(BMI = function(x) paste0("BMI = ", x))
    ) +
    labs(
      x = "Time (years)",
      y = "Risk"
    ) +
    guides(
      color = guide_legend(
        nrow = 1,
        override.aes = list(
          linewidth = 1.6,
          linetype = unname(lt_map)
        )
      )
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.justification = "center",
      legend.key.width = unit(1.8, "cm"),
      legend.key.height = unit(0.7, "cm"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

## -----------------------------
## Run
## -----------------------------
df_rt <- risk_vs_time_multi_bmi(
  df = df,
  theta_hat = theta_hat,
  bmi_vec = bmi_vec,
  SCALE = SCALE,
  h_kernel = h_kernel
)

p2 <- plot_risk_vs_time_facet_bmi(df_rt)
print(p2)

ggsave(
  filename = file.path(fig_dir, "risk_vs_time_fixed_bmi.png"),
  plot = p2,
  width = 8.5,
  height = 3
)