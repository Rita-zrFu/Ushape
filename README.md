# Ushape

Reproducibility code for the simulation study in the manuscript

> **Rank-Based Estimation of U-Shaped Biomarker Risk Curves and Critical Points for Time-to-Event Outcomes**
> (submitted to *Biometrics*, Biometric Methodology section)

The repository contains the R and C++ code that generated the 23-scenario simulation results reported in the paper (main-text Table 1 and Figures 2--3), together with the SLURM dispatcher used to run the study on an HPC cluster, and an R script that reproduces the two figures from the aggregated CSV outputs.

---

## Repository layout

```
Ushape/
├── README.md
├── SCENARIOS.csv                 source of truth for the 23 simulation scenarios
├── .gitignore
├── code/
│   ├── simulation_one_rep.R      per-replicate simulation runner
│   ├── myCindex.cpp              Rcpp C-index kernel (hard + sigmoid-smoothed)
│   ├── summarize_all_to_csv.R    aggregate per-rep RDS files into 7 CSVs
│   └── build_figs.R              rebuild fig_sim_xc.pdf and fig_sim_cindex.pdf
├── cluster/
│   └── dispatch.sh               SLURM array dispatcher (one job per scenario)
├── real data_lilly.R             legacy real-data prototype (see below)
└── bootstrap_realdata_Lilly.R    legacy stub
```

---

## Simulation design

The study spans **23 scenarios × 6 methods × 1000 replications**. Scenarios vary in

| Factor | Levels |
|---|---|
| sample size `n` | 200, 500, 1000 |
| error `epsilon` | `N(0,9)` (`norm`) or centered min-extreme-value with SD 3 (`ev`) |
| link `G` | logistic (`10 * plogis(s, scale = 5)`) or exponential (`exp((s+1)/5)`) |
| target censoring rate | 15%, 30%, 40% |
| `beta_1` | 1 (milder U-shape) or 2 (steeper U-shape) |

`SCENARIOS.csv` is the single source of truth: each row provides `scenario_id`, `n`, `ei`, `G_type`, `censor_lo`, `censor_hi`, `beta1`, `target_cens_rate`, evaluation times `t0_S` and `t0_CR`, and the corresponding true values (`cp_true`, `s_t0_true`, `crlb_true`, `crub_true`, `sigma_h`).

For each scenario the following six methods are fit:

- **MCE** -- rank-based Maximum C-index Estimator (proposed method)
- **Cox+spline** -- Cox with a natural spline in `X`
- **Cox+quadratic** -- Cox with `X` and `X^2`
- **Cox+quadratic w/ X:Z interactions**
- **Cox+spline w/ X:Z interactions**
- **Oracle** -- Cox with the true change point known

Point estimation for the MCE uses `DEoptim` on the hard C-index; standard errors come from a sigmoid-smoothed bootstrap fitted with L-BFGS-B (500 bootstrap draws, bandwidth `sigma_h = round(n^(1/4))`). Delta-method and percentile CIs for the critical point are both saved.

Full design details are given in Sections 2 and 6 of the manuscript.

---

## How to reproduce

### Dependencies

- R >= 4.0 with packages `DEoptim`, `survival`, `SurvMetrics`, `Rcpp`, `splines`, `dplyr`, `ggplot2`, `patchwork`
- C++ toolchain compatible with `Rcpp::sourceCpp` (gcc 11+ or equivalent)

### Single-replicate local run (smoke test)

```
Rscript code/simulation_one_rep.R S07 1 1000 norm logistic 0.3000 8.3000 TRUE 2
```

Arguments are `scenario_id seed n ei G_type censor_lo censor_hi run_competitors [beta1]`. The script writes one RDS list to `results/S07/rep_0001_j<batch_tag>.rds` containing MCE estimates, bootstrap covariance, all five competitor results, and full provenance (`scenario_id`, `seed`, `batch_tag`, `timestamp`, `n`, `ei`, `G_type`, `beta1`, censoring bounds).

### Full study via SLURM

```
bash cluster/dispatch.sh
```

Environment overrides: `N_REPS` (default 1100), `WALL` (default 4:00:00), `MEM` (default 4G), `RESULTS_DIR`, `SCENARIOS_CSV`, `MYCINDEX_CPP_PATH`, `RUN_COMPETITORS`. The dispatcher submits 23 array jobs of 1100 tasks each (~25,300 total) and keeps the first 1000 successful replicates per scenario at summarization time.

### Aggregate

```
Rscript code/summarize_all_to_csv.R results csv_output SCENARIOS.csv
```

Produces `csv_output/{sim_meta, sim_params, sim_xc, sim_cindex, sim_survival, sim_cr, sim_boot_diag}.csv`.

### Rebuild figures

```
Rscript code/build_figs.R csv_output SCENARIOS.csv figures
```

Produces `figures/fig_sim_xc.pdf` and `figures/fig_sim_cindex.pdf`. Both are three-panel layouts (main grid / censoring sensitivity / milder-U-shape sensitivity) with method colors and shapes shared across the two files.

---

## Output provenance

Every per-replicate RDS embeds `scenario_id`, `seed`, `batch_tag` (SLURM job ID when available, otherwise a UTC timestamp), `timestamp`, and the exact scenario parameters. The summarizer checks provenance consistency across replicates before aggregating; scenarios with mismatched provenance are skipped rather than silently merged.

---

## Real-data analysis (legacy)

`real data_lilly.R` and `bootstrap_realdata_Lilly.R` are earlier working files retained for reference. The UK Biobank analysis reported in Section 7 of the current manuscript is not distributed in this repository because the data are subject to UK Biobank access terms.

---

## Contact

Corresponding author: Yuxin (Daisy) Zhu, Department of Biostatistics, Johns Hopkins Bloomberg School of Public Health (`daisy@jhu.edu`).

Co-authors: Zhirui Fu, Mei-Cheng Wang, Yu Du.

## License

TBD.
