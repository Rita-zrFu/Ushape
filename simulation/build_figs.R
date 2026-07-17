#!/usr/bin/env Rscript
# build_figs.R
# Reproduce the two USR simulation figures from the aggregated CSVs:
#   figures/fig_sim_xc.pdf     -- Xc bias, mean +/- 1.96 * ESE, as filled boxes
#   figures/fig_sim_cindex.pdf -- test-set Harrell C-index, mean +/- 1.96 * SD
#
# Panel layout (matches the manuscript captions):
#   (a) Main grid: beta1 = 2 and 30% censoring, faceted 2x2 by (G, epsilon),
#       x-axis = sample size n.
#   (b) Censoring sensitivity at epsilon = N(0,9) and beta1 = 2: the additional
#       15% and 50% censoring scenarios (target_cens_rate 0.15 and 0.40) for
#       each G. Faceted by G; x-axis is (n, censoring%) combinations.
#   (c) Milder U-shape at beta1 = 1 (logistic G, min-EV epsilon). Single panel,
#       x-axis = n.
#
# Method colors and shapes are shared across the two figures.
#
# Inputs (defaults):
#   csv_output/sim_xc.csv, csv_output/sim_cindex.csv
#
# The 23-scenario table is embedded below so that this script has no external
# CSV dependency other than the aggregated MCE/competitor CSVs.
#
# Usage:
#   Rscript build_figs.R [csv_dir] [out_dir]

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
csv_dir <- if (length(args) >= 1) args[1] else "csv_output"
out_dir <- if (length(args) >= 2) args[2] else "figures"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Scenario table (matches cluster/dispatch_full.sh) ---
scen <- read.table(text = "
scenario_id n    ei      G_type   censor_lo censor_hi beta1 target_cens_rate
S01         200  norm    logistic 1.9800    9.9800    2     0.15
S02         200  norm    logistic 0.3000    8.3000    2     0.30
S03         200  norm    logistic 0.1100    4.9100    2     0.40
S04         500  norm    logistic 1.9800    9.9800    2     0.15
S05         500  norm    logistic 0.3000    8.3000    2     0.30
S06         500  norm    logistic 0.1100    4.9100    2     0.40
S07         1000 norm    logistic 0.3000    8.3000    2     0.30
S08         200  ev      logistic 0.3300    8.3300    2     0.30
S09         500  ev      logistic 0.3300    8.3300    2     0.30
S10         1000 ev      logistic 0.3300    8.3300    2     0.30
S11         200  norm    exp      0.0400    4.0400    2     0.15
S12         200  norm    exp      0.0100    1.9100    2     0.30
S13         200  norm    exp      0.0300    0.8300    2     0.40
S14         500  norm    exp      0.0400    4.0400    2     0.15
S15         500  norm    exp      0.0100    1.9100    2     0.30
S16         500  norm    exp      0.0300    0.8300    2     0.40
S17         1000 norm    exp      0.0100    1.9100    2     0.30
S18         200  ev      exp      0.0400    1.8400    2     0.30
S19         500  ev      exp      0.0400    1.8400    2     0.30
S20         1000 ev      exp      0.0400    1.8400    2     0.30
S21         200  ev      logistic 1.9000    7.9000    1     0.30
S22         500  ev      logistic 1.9000    7.9000    1     0.30
S23         1000 ev      logistic 1.9000    7.9000    1     0.30
",
header = TRUE, stringsAsFactors = FALSE)

xc <- read.csv(file.path(csv_dir, "sim_xc.csv"),     stringsAsFactors = FALSE)
ci <- read.csv(file.path(csv_dir, "sim_cindex.csv"), stringsAsFactors = FALSE)

# --- Column-name compatibility ---
# Older csv_output/ files use xc_mean/CI_cov_*; current summarizer uses
# bias/ECP_*.  Only 'bias' and 'ESE' are needed here (both formats have them).
if (!all(c("bias", "ESE") %in% names(xc))) {
  stop("sim_xc.csv missing required columns: bias and ESE")
}
# C-index: current summarizer emits c_test/c_test_sd; older format used
# c_test_mean/c_test_sd. Normalize to c_test/c_test_sd.
if ("c_test_mean" %in% names(ci) && !"c_test" %in% names(ci)) {
  ci$c_test <- ci$c_test_mean
}
if (!all(c("c_test", "c_test_sd") %in% names(ci))) {
  stop("sim_cindex.csv missing required columns: c_test and c_test_sd")
}

# --- Method display (order, color, shape, label) ---
methods_xc <- c("MCE", "spline", "quadratic", "spline_interact", "quad_interact")
methods_ci <- c("MCE", "spline", "quadratic", "spline_interact", "quad_interact", "oracle")
labels_all <- c(MCE = "MCE",
                spline = "Spline",
                quadratic = "Quadratic",
                spline_interact = "Spline×Z",
                quad_interact = "Quad.×Z",
                oracle = "Oracle")
# Palette approximates the manuscript figure.
colors_all <- c(MCE             = "#2CA02C",
                spline          = "#E58139",
                quadratic       = "#9E9AC8",
                spline_interact = "#E377C2",
                quad_interact   = "#8C6D31",
                oracle          = "#7F7F7F")
shapes_all <- c(MCE             = 16,
                spline          = 17,
                quadratic       = 15,
                spline_interact = 18,
                quad_interact   = 8,
                oracle          = 4)

# --- Merge CSVs with scenario metadata ---
merge_scen <- function(df) {
  df <- merge(df, scen, by.x = "scenario", by.y = "scenario_id", all.x = TRUE)
  df$G_lab  <- ifelse(df$G_type == "logistic", "logistic", "exp")
  df$ei_lab <- ifelse(df$ei == "norm", "N(0,9)", "min-EV")
  df$GE_lab <- paste0(df$G_lab, ", ", df$ei_lab)
  df$G_facet <- paste0("G = ", df$G_lab)
  # Manuscript figure labels 0.40 as "50%"; match this here.
  df$cens_lab_pct <- ifelse(df$target_cens_rate == 0.40, "50%",
                            paste0(round(100 * df$target_cens_rate), "%"))
  df$nc_lab <- paste0("n=", df$n, ", ", df$cens_lab_pct)
  df
}
xc <- merge_scen(xc)
ci <- merge_scen(ci)

# --- Panel assignment ---
panel_a_df <- function(df) df[df$beta1 == 2 & df$target_cens_rate == 0.30, ]
panel_b_df <- function(df) df[df$beta1 == 2 & df$target_cens_rate != 0.30 & df$ei == "norm", ]
panel_c_df <- function(df) df[df$beta1 == 1, ]

GE_LEVELS <- c("logistic, N(0,9)", "logistic, min-EV", "exp, N(0,9)", "exp, min-EV")
G_LEVELS  <- c("G = logistic", "G = exp")

# --- Shared plotting primitives ---
theme_biom <- theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(size = 9),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.key.width = unit(1.2, "lines"),
        plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8))

# Dodge width for numeric x (used in panels (a) and (c)).
DODGE_N <- position_dodge(width = 55)
DODGE_C <- position_dodge(width = 0.6)  # categorical x (panel (b))

make_plot <- function(df, methods, y_var, y_sd, y_lab,
                       facet_type = c("2x2", "1x2", "none"),
                       x_var, dodge, include_hline = FALSE, title = NULL) {
  facet_type <- match.arg(facet_type)
  df <- df[df$method %in% methods, ]
  if (nrow(df) == 0) stop("No rows to plot after filtering methods.")
  df$method <- factor(df$method, levels = methods)
  df$.y  <- df[[y_var]]
  df$.lo <- df[[y_var]] - 1.96 * df[[y_sd]]
  df$.hi <- df[[y_var]] + 1.96 * df[[y_sd]]

  p <- ggplot(df, aes(x = .data[[x_var]], y = .y,
                      color = method, fill = method, shape = method,
                      group = method))
  if (include_hline) {
    p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey55")
  }
  # Box (fill = light method color; no center line).
  p <- p +
    geom_crossbar(aes(ymin = .lo, ymax = .hi),
                  width = if (identical(dodge, DODGE_C)) 0.5 else 45,
                  fatten = 0, position = dodge, alpha = 0.30,
                  linewidth = 0.45, show.legend = TRUE) +
    geom_point(position = dodge, size = 2.2, stroke = 0.7)

  if (facet_type == "2x2") {
    df$GE_lab <- factor(df$GE_lab, levels = GE_LEVELS)
    p <- p + facet_wrap(~ factor(GE_lab, levels = GE_LEVELS), ncol = 2)
  } else if (facet_type == "1x2") {
    df$G_facet <- factor(df$G_facet, levels = G_LEVELS)
    p <- p + facet_wrap(~ factor(G_facet, levels = G_LEVELS), ncol = 2)
  }

  p <- p +
    scale_color_manual(name = "Method",
                       values = colors_all[methods],
                       labels = labels_all[methods]) +
    scale_fill_manual(name  = "Method",
                       values = colors_all[methods],
                       labels = labels_all[methods]) +
    scale_shape_manual(name = "Method",
                       values = shapes_all[methods],
                       labels = labels_all[methods]) +
    labs(y = y_lab, title = title) +
    theme_biom

  if (x_var == "n") {
    p <- p +
      scale_x_continuous(breaks = sort(unique(df$n))) +
      labs(x = "Sample size n")
  } else {
    lv <- c("n=200, 15%", "n=200, 50%", "n=500, 15%", "n=500, 50%")
    p <- p +
      scale_x_discrete(limits = lv) +
      labs(x = "Sample size n / censoring %")
  }
  p
}

# --- Xc figure -------------------------------------------------------------
y_lab_xc <- expression("Bias of "*hat(X)[c])

pa_xc <- make_plot(panel_a_df(xc), methods_xc, "bias", "ESE", y_lab_xc,
                   facet_type = "2x2", x_var = "n", dodge = DODGE_N,
                   include_hline = TRUE,
                   title = expression("(a) Main grid: "*beta[1]*" = 2, 30% censoring"))
pb_xc <- make_plot(panel_b_df(xc), methods_xc, "bias", "ESE", y_lab_xc,
                   facet_type = "1x2", x_var = "nc_lab", dodge = DODGE_C,
                   include_hline = TRUE,
                   title = expression("(b) Censoring sensitivity ("*epsilon*" = N(0,9), "*beta[1]*" = 2)"))
pc_xc <- make_plot(panel_c_df(xc), methods_xc, "bias", "ESE", y_lab_xc,
                   facet_type = "none", x_var = "n", dodge = DODGE_N,
                   include_hline = TRUE,
                   title = expression("(c) "*beta[1]*" = 1 sensitivity (logistic, min-EV, 30% censoring)"))

fig_xc <- pa_xc / pb_xc / pc_xc +
  plot_layout(heights = c(2.6, 1.1, 1.0), guides = "collect") &
  theme(legend.position = "bottom")

# --- C-index figure -------------------------------------------------------
y_lab_ci <- "Test-set Harrell's C-index"

pa_ci <- make_plot(panel_a_df(ci), methods_ci, "c_test", "c_test_sd", y_lab_ci,
                   facet_type = "2x2", x_var = "n", dodge = DODGE_N,
                   include_hline = FALSE,
                   title = expression("(a) Main grid: "*beta[1]*" = 2, 30% censoring"))
pb_ci <- make_plot(panel_b_df(ci), methods_ci, "c_test", "c_test_sd", y_lab_ci,
                   facet_type = "1x2", x_var = "nc_lab", dodge = DODGE_C,
                   include_hline = FALSE,
                   title = expression("(b) Censoring sensitivity ("*epsilon*" = N(0,9), "*beta[1]*" = 2)"))
pc_ci <- make_plot(panel_c_df(ci), methods_ci, "c_test", "c_test_sd", y_lab_ci,
                   facet_type = "none", x_var = "n", dodge = DODGE_N,
                   include_hline = FALSE,
                   title = expression("(c) "*beta[1]*" = 1 sensitivity (logistic, min-EV, 30% censoring)"))

fig_ci <- pa_ci / pb_ci / pc_ci +
  plot_layout(heights = c(2.6, 1.1, 1.0), guides = "collect") &
  theme(legend.position = "bottom")

# --- Save ------------------------------------------------------------------
ggsave(file.path(out_dir, "fig_sim_xc.pdf"),     fig_xc,
       width = 7.5, height = 8.5, units = "in")
ggsave(file.path(out_dir, "fig_sim_cindex.pdf"), fig_ci,
       width = 7.5, height = 8.5, units = "in")

cat(sprintf("Wrote: %s\n", file.path(out_dir, "fig_sim_xc.pdf")))
cat(sprintf("Wrote: %s\n", file.path(out_dir, "fig_sim_cindex.pdf")))
