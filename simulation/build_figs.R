#!/usr/bin/env Rscript
# build_figs.R
# Build the two USR simulation figures used in the manuscript:
#   fig_sim_xc.pdf     -- X_c bias, mean +/- 1.96 * ESE, filled boxes
#   fig_sim_cindex.pdf -- test-set Harrell's C-index, mean +/- 1.96 * SD
#
# Each figure has three stacked sections (heights ~ 2.6 : 1.1 : 1.0):
#   (a) Main grid:  2x2 by (G, epsilon), x = n, methods dodged.   [12 scenarios]
#   (b) Censoring:  1x2 by G (epsilon = N only), x = "n / cens %"  [ 8 scenarios]
#                   methods dodged.
#   (c) beta_1 = 1: 1x1 (logistic, min-EV, 30% cens),              [ 3 scenarios]
#                   x = n, methods dodged.
#
# Inputs: csv_dir/{sim_xc.csv, sim_cindex.csv, sim_meta.csv}
# Outputs: out_dir/fig_sim_xc.pdf, out_dir/fig_sim_cindex.pdf
#
# Usage:
#   Rscript build_figs.R [csv_dir] [out_dir]
# Defaults:
#   csv_dir = csv_output
#   out_dir = figures

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
csvdir <- if (length(args) >= 1) args[1] else "csv_output"
outdir <- if (length(args) >= 2) args[2] else "figures"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

xc   <- read.csv(file.path(csvdir, "sim_xc.csv"),     stringsAsFactors = FALSE)
ci   <- read.csv(file.path(csvdir, "sim_cindex.csv"), stringsAsFactors = FALSE)
meta <- read.csv(file.path(csvdir, "sim_meta.csv"),   stringsAsFactors = FALSE)

# ---- Method palette / shapes (Dark2) ----
ALL_METHODS <- c("MCE","spline","quadratic","spline_interact","quad_interact","oracle")
M_LABELS    <- c("MCE","Spline","Quadratic","Spline×Z","Quad.×Z","Oracle")
method_cols <- setNames(
  c("#1B9E77","#D95F02","#7570B3","#E7298A","#A6761D","#666666"), M_LABELS)
method_shapes <- setNames(c(16, 17, 15, 18, 8, 4), M_LABELS)

# ---- Scenario sets ----
MAIN_SIDS <- c("S02","S05","S07","S08","S09","S10","S12","S15","S17","S18","S19","S20")
CENS_SIDS <- c("S01","S04","S03","S06","S11","S14","S13","S16")
B1_SIDS   <- c("S21","S22","S23")

# ---- Helpers ----
PANEL_LEVELS_MAIN <- c(
  "logistic, N(0,9)", "logistic, min-EV",
  "exp, N(0,9)",      "exp, min-EV")
panel_lab_main <- function(G, ei) {
  ep <- ifelse(ei == "norm", "N(0,9)", "min-EV")
  paste0(G, ", ", ep)
}
panel_lab_G <- function(G) ifelse(G == "logistic", "G = logistic", "G = exp")

attach_meta <- function(df) {
  df %>% inner_join(meta[, c("scenario","n","G_type","ei","cens_rate","beta1")],
                    by = "scenario")
}

# ---- Build data frames ----
xc_methods_feas <- setdiff(ALL_METHODS, "oracle")
ci_methods_all  <- c(ALL_METHODS[1:5], "oracle_v2")

# X_c data (no oracle)
mk_xc_df <- function(sids) {
  xc %>%
    filter(scenario %in% sids, method %in% xc_methods_feas) %>%
    mutate(method_lab = factor(method,
                               levels = ALL_METHODS[1:5],
                               labels = M_LABELS[1:5])) %>%
    attach_meta()
}
xc_main_df <- mk_xc_df(MAIN_SIDS) %>%
  mutate(panel = factor(panel_lab_main(G_type, ei), levels = PANEL_LEVELS_MAIN),
         n_lab = factor(n, levels = c(200, 500, 1000)))
xc_cens_df <- mk_xc_df(CENS_SIDS) %>%
  mutate(panel = factor(panel_lab_G(G_type),
                        levels = c("G = logistic", "G = exp")),
         cens_pct = round(cens_rate * 100),
         x_lab = factor(sprintf("n=%d, %d%%", n, cens_pct),
                        levels = c("n=200, 15%", "n=200, 50%",
                                   "n=500, 15%", "n=500, 50%")))
xc_b1_df <- mk_xc_df(B1_SIDS) %>%
  mutate(n_lab = factor(n, levels = c(200, 500, 1000)))

# C-index data (all 6 methods; oracle_v2 relabeled to oracle for display)
mk_ci_df <- function(sids) {
  ci %>%
    filter(scenario %in% sids, method %in% ci_methods_all) %>%
    mutate(method = ifelse(method == "oracle_v2", "oracle", method),
           method_lab = factor(method, levels = ALL_METHODS, labels = M_LABELS)) %>%
    attach_meta()
}
ci_main_df <- mk_ci_df(MAIN_SIDS) %>%
  mutate(panel = factor(panel_lab_main(G_type, ei), levels = PANEL_LEVELS_MAIN),
         n_lab = factor(n, levels = c(200, 500, 1000)))
ci_cens_df <- mk_ci_df(CENS_SIDS) %>%
  mutate(panel = factor(panel_lab_G(G_type),
                        levels = c("G = logistic", "G = exp")),
         cens_pct = round(cens_rate * 100),
         x_lab = factor(sprintf("n=%d, %d%%", n, cens_pct),
                        levels = c("n=200, 15%", "n=200, 50%",
                                   "n=500, 15%", "n=500, 50%")))
ci_b1_df <- mk_ci_df(B1_SIDS) %>%
  mutate(n_lab = factor(n, levels = c(200, 500, 1000)))

# ---- Shared theme ----
shared_theme <- theme_bw(base_size = 8.5) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey92", color = "grey80"),
        strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(face = "bold", size = 9,
                                  margin = margin(b = 4)),
        plot.margin = margin(t = 6, r = 6, b = 4, l = 6),
        legend.key.size = unit(0.45, "cm"),
        legend.text = element_text(size = 8),
        axis.text.x = element_text(size = 7.5))

# ============================================================
# Topic builder: X_c bias or C-index in the same three-panel design
# ============================================================
build_topic <- function(main_df, cens_df, b1_df,
                        y_expr_main, y_expr_cens, y_expr_b1,
                        ymin_expr, ymax_expr,
                        y_lab, drop_oracle, is_xc = FALSE) {
  drop_arg <- if (drop_oracle) TRUE else FALSE
  cols_used <- if (drop_oracle) method_cols[1:5] else method_cols

  # ---- (a) main 2x2 grid ----
  pA <- ggplot(main_df, aes(x = n_lab, y = !!y_expr_main,
                            color = method_lab, fill = method_lab,
                            shape = method_lab))
  if (is_xc) pA <- pA + geom_hline(yintercept = 0, linetype = "dashed",
                                   color = "grey50", linewidth = 0.4)
  pA <- pA +
    geom_crossbar(aes(ymin = !!ymin_expr, ymax = !!ymax_expr),
                  width = 0.5, fatten = 0, linewidth = 0.4, alpha = 0.2,
                  position = position_dodge(width = 0.55)) +
    geom_point(size = 1.6, position = position_dodge(width = 0.55), stroke = 0.55) +
    facet_wrap(~ panel, ncol = 2) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    labs(x = "Sample size n", y = y_lab,
         title = "(a) Main grid: β₁ = 2, 30% censoring") +
    shared_theme

  # ---- (b) censoring sensitivity 1x2 (N error only) ----
  pB <- ggplot(cens_df, aes(x = x_lab, y = !!y_expr_cens,
                            color = method_lab, fill = method_lab,
                            shape = method_lab))
  if (is_xc) pB <- pB + geom_hline(yintercept = 0, linetype = "dashed",
                                   color = "grey50", linewidth = 0.4)
  pB <- pB +
    geom_crossbar(aes(ymin = !!ymin_expr, ymax = !!ymax_expr),
                  width = 0.55, fatten = 0, linewidth = 0.4, alpha = 0.2,
                  position = position_dodge(width = 0.6)) +
    geom_point(size = 1.6, position = position_dodge(width = 0.6), stroke = 0.55) +
    facet_wrap(~ panel, ncol = 2) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    labs(x = "Sample size n / censoring %", y = y_lab,
         title = "(b) Censoring sensitivity (ε = N(0,9), β₁ = 2)") +
    shared_theme

  # ---- (c) beta_1 = 1 single panel ----
  pC <- ggplot(b1_df, aes(x = n_lab, y = !!y_expr_b1,
                          color = method_lab, fill = method_lab,
                          shape = method_lab))
  if (is_xc) pC <- pC + geom_hline(yintercept = 0, linetype = "dashed",
                                   color = "grey50", linewidth = 0.4)
  pC <- pC +
    geom_crossbar(aes(ymin = !!ymin_expr, ymax = !!ymax_expr),
                  width = 0.5, fatten = 0, linewidth = 0.4, alpha = 0.2,
                  position = position_dodge(width = 0.55)) +
    geom_point(size = 1.6, position = position_dodge(width = 0.55), stroke = 0.55) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    labs(x = "Sample size n", y = y_lab,
         title = "(c) β₁ = 1 sensitivity (logistic, min-EV, 30% censoring)") +
    shared_theme

  # Apply shared scales to all three panels via patchwork &
  combined <- (pA / pB / pC) +
    plot_layout(heights = c(2.6, 1.1, 1.0), guides = "collect") &
    scale_color_manual(values = method_cols, name = "Method", drop = drop_arg) &
    scale_fill_manual(values = method_cols,  name = "Method", drop = drop_arg) &
    scale_shape_manual(values = method_shapes, name = "Method", drop = drop_arg) &
    guides(fill = "none", color = "none",
           shape = guide_legend(nrow = 1,
                                override.aes = list(size = 2.8,
                                                    color = cols_used))) &
    theme(legend.position = "bottom")
  combined
}

# ============================================================
# Figure 1: X_c bias
# ============================================================
fig_xc <- build_topic(
  xc_main_df, xc_cens_df, xc_b1_df,
  y_expr_main = quote(bias), y_expr_cens = quote(bias), y_expr_b1 = quote(bias),
  ymin_expr = quote(bias - 1.96 * ESE),
  ymax_expr = quote(bias + 1.96 * ESE),
  y_lab = expression("Bias of " * hat(X)[c]),
  drop_oracle = TRUE, is_xc = TRUE)

out_pdf <- file.path(outdir, "fig_sim_xc.pdf")
ggsave(out_pdf, fig_xc, width = 7.5, height = 8.5, units = "in",
       device = cairo_pdf)
cat(sprintf("Wrote %s\n", out_pdf))

# ============================================================
# Figure 2: Test-set C-index
# ============================================================
fig_ci <- build_topic(
  ci_main_df, ci_cens_df, ci_b1_df,
  y_expr_main = quote(c_test_mean), y_expr_cens = quote(c_test_mean),
  y_expr_b1 = quote(c_test_mean),
  ymin_expr = quote(c_test_mean - 1.96 * c_test_sd),
  ymax_expr = quote(c_test_mean + 1.96 * c_test_sd),
  y_lab = expression("Test-set Harrell's " * italic(C) * "-index"),
  drop_oracle = FALSE, is_xc = FALSE)

out_pdf <- file.path(outdir, "fig_sim_cindex.pdf")
ggsave(out_pdf, fig_ci, width = 7.5, height = 8.5, units = "in",
       device = cairo_pdf)
cat(sprintf("Wrote %s\n", out_pdf))
