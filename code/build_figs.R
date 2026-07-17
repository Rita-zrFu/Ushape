#!/usr/bin/env Rscript
# build_figs.R
# Reproduce the two USR simulation figures from the aggregated CSVs:
#   figures/fig_sim_xc.pdf     -- Xc bias with mean +/- 1.96 * ESE
#   figures/fig_sim_cindex.pdf -- test-set Harrell C-index with mean +/- 1.96 * SD
#
# Layout matches the manuscript captions:
#   Panel (a) main grid: beta1 = 2 and target_cens_rate = 0.30, faceted by
#             (G, epsilon) cells, sample size on the x-axis.
#   Panel (b) censoring sensitivity at epsilon = N(0,9) and beta1 = 2:
#             the 15% and 40% censoring scenarios for each G.
#   Panel (c) milder U-shape at beta1 = 1 (single (G,epsilon) cell).
#
# Method colors and shapes are shared across the two figures.
#
# Inputs (defaults):
#   csv_output/sim_xc.csv, csv_output/sim_cindex.csv
#   SCENARIOS.csv (for panel assignment via beta1, target_cens_rate, ei, G_type, n)
#
# Usage:
#   Rscript code/build_figs.R [csv_dir] [scenarios_csv] [out_dir]

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
csv_dir       <- if (length(args) >= 1) args[1] else "csv_output"
scenarios_csv <- if (length(args) >= 2) args[2] else "SCENARIOS.csv"
out_dir       <- if (length(args) >= 3) args[3] else "figures"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

scen <- read.csv(scenarios_csv, stringsAsFactors = FALSE)
xc   <- read.csv(file.path(csv_dir, "sim_xc.csv"),     stringsAsFactors = FALSE)
ci   <- read.csv(file.path(csv_dir, "sim_cindex.csv"), stringsAsFactors = FALSE)

# --- Column-name compatibility ---
# Older csv_output/ files use xc_mean/CI_cov_*; the current summarizer uses
# bias/ECP_*.  We only need `bias` (present in both) and an SE-like column.
if (!"ESE" %in% names(xc)) stop("sim_xc.csv missing ESE column")
if (!"bias" %in% names(xc)) stop("sim_xc.csv missing bias column")

# C-index: current summarizer uses c_test/c_test_sd; older format used
# c_test_mean/c_test_sd.  Handle both.
if ("c_test_mean" %in% names(ci)) {
  ci$c_test <- ci$c_test_mean
} else if (!"c_test" %in% names(ci)) {
  stop("sim_cindex.csv missing c_test column")
}
if (!"c_test_sd" %in% names(ci)) stop("sim_cindex.csv missing c_test_sd column")

# --- Method display (order, color, shape) ---
methods_xc <- c("MCE", "spline", "quadratic", "quad_interact", "spline_interact")
methods_ci <- c("MCE", "spline", "quadratic", "quad_interact", "spline_interact", "oracle")
labels_all <- c(MCE = "MCE",
                spline = "Cox+spline",
                quadratic = "Cox+quadratic",
                quad_interact = "Cox+quad w/ interactions",
                spline_interact = "Cox+spline w/ interactions",
                oracle = "Oracle")
colors_all <- c(MCE = "#D62728",
                spline = "#1F77B4",
                quadratic = "#2CA02C",
                quad_interact = "#9467BD",
                spline_interact = "#FF7F0E",
                oracle = "#7F7F7F")
shapes_all <- c(MCE = 16,
                spline = 17,
                quadratic = 15,
                quad_interact = 18,
                spline_interact = 25,
                oracle = 4)

# --- Merge CSVs with scenario metadata ---
merge_scen <- function(df) {
  df <- merge(df, scen, by.x = "scenario", by.y = "scenario_id", all.x = TRUE)
  df$G_lab   <- ifelse(df$G_type == "logistic", "G: logistic", "G: exponential")
  df$ei_lab  <- ifelse(df$ei == "norm", "epsilon: N(0,9)", "epsilon: min-EV")
  df$cens_lab_pct <- paste0(round(100 * df$target_cens_rate), "%")
  df
}
xc <- merge_scen(xc)
ci <- merge_scen(ci)

# --- Panel assignment ---
panel_a <- function(df) df[df$beta1 == 2 & df$target_cens_rate == 0.30, ]
panel_b <- function(df) df[df$beta1 == 2 & df$target_cens_rate != 0.30 & df$ei == "norm", ]
panel_c <- function(df) df[df$beta1 == 1, ]

# --- Plotters ---
theme_biom <- theme_bw(base_size = 9) +
  theme(strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 9, face = "bold"))

plot_xc_panel <- function(df, facet_formula, title, x_var = "n",
                          dodge_w = 45) {
  df <- df[df$method %in% methods_xc, ]
  df$method <- factor(df$method, levels = methods_xc)
  df$lo <- df$bias - 1.96 * df$ESE
  df$hi <- df$bias + 1.96 * df$ESE
  pd <- position_dodge(width = dodge_w)
  ggplot(df, aes(x = .data[[x_var]], y = bias,
                 color = method, shape = method, group = method)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_errorbar(aes(ymin = lo, ymax = hi), position = pd, width = 30) +
    geom_point(position = pd, size = 2) +
    facet_grid(facet_formula) +
    scale_color_manual(values = colors_all[methods_xc], labels = labels_all[methods_xc]) +
    scale_shape_manual(values = shapes_all[methods_xc], labels = labels_all[methods_xc]) +
    scale_x_continuous(breaks = sort(unique(df[[x_var]]))) +
    labs(x = "sample size n", y = expression(hat(X)[c]~"bias"), title = title) +
    theme_biom
}

plot_ci_panel <- function(df, facet_formula, title, x_var = "n",
                          dodge_w = 45) {
  df <- df[df$method %in% methods_ci, ]
  df$method <- factor(df$method, levels = methods_ci)
  df$lo <- df$c_test - 1.96 * df$c_test_sd
  df$hi <- df$c_test + 1.96 * df$c_test_sd
  pd <- position_dodge(width = dodge_w)
  ggplot(df, aes(x = .data[[x_var]], y = c_test,
                 color = method, shape = method, group = method)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), position = pd, width = 30) +
    geom_point(position = pd, size = 2) +
    facet_grid(facet_formula) +
    scale_color_manual(values = colors_all[methods_ci], labels = labels_all[methods_ci]) +
    scale_shape_manual(values = shapes_all[methods_ci], labels = labels_all[methods_ci]) +
    scale_x_continuous(breaks = sort(unique(df[[x_var]]))) +
    labs(x = "sample size n", y = "test-set C-index", title = title) +
    theme_biom
}

build_three_panel <- function(df, plotter) {
  pa <- plotter(panel_a(df),
                as.formula("ei_lab ~ G_lab"),
                "(a) Main grid: beta1 = 2, 30% censoring")
  pb <- plotter(panel_b(df),
                as.formula("cens_lab_pct ~ G_lab"),
                "(b) Censoring sensitivity: epsilon = N(0,9), beta1 = 2")
  pc <- plotter(panel_c(df),
                as.formula(". ~ G_lab"),
                "(c) Milder U-shape: beta1 = 1, epsilon = min-EV")
  pa / pb / pc +
    plot_layout(heights = c(2.6, 1.1, 1.0), guides = "collect") &
    theme(legend.position = "bottom")
}

# --- Build both figures ---
fig_xc <- build_three_panel(xc, plot_xc_panel)
fig_ci <- build_three_panel(ci, plot_ci_panel)

ggsave(file.path(out_dir, "fig_sim_xc.pdf"),     fig_xc,
       width = 7.5, height = 8.5, units = "in")
ggsave(file.path(out_dir, "fig_sim_cindex.pdf"), fig_ci,
       width = 7.5, height = 8.5, units = "in")

cat(sprintf("Wrote: %s\n", file.path(out_dir, "fig_sim_xc.pdf")))
cat(sprintf("Wrote: %s\n", file.path(out_dir, "fig_sim_cindex.pdf")))
