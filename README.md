# Ushape

Code accompanying the manuscript

> **Rank-Based Estimation of U-Shaped Biomarker Risk Curves and Critical Points for Time-to-Event Outcomes**

This repository contains the simulation and real-data analysis code used in the paper. It does **not** contain results (per-replicate RDS files, aggregated CSVs, or figure PDFs); those are regenerated locally by running the scripts described below.

---

## Repository layout

```
Ushape/
├── README.md
├── .gitignore
├── simulation/
│   ├── simulation_one_rep.R          per-replicate simulation runner
│   ├── myCindex.cpp                   Rcpp C-index kernel (hard + smoothed)
│   ├── summarize_all_to_csv.R         aggregate per-rep RDS → 7 CSVs
│   ├── reprocess_oracle_v2.R          add corrected oracle Cox to per-rep RDS
│   ├── build_figs.R                   rebuild fig_sim_xc.pdf and fig_sim_cindex.pdf
│   └── cluster/
│       ├── dispatch_full.sh                    SLURM array dispatcher
│       └── dispatch_reprocess_oracle_v2.sh     SLURM oracle reprocessing dispatcher
└── real_data_analysis/
    ├── real data_lilly.R              UK Biobank / Lilly analysis prototype
    └── bootstrap_realdata_Lilly.R     bootstrap stub
```

Nothing under `simulation/results/`, `simulation/csv_output/`, `simulation/figures/`, or `simulation/cluster/logs/` is tracked in the repository; these are output directories the pipeline populates locally.

---

## Simulation study

The simulation reported in the manuscript spans **23 scenarios × 6 methods × 1000 replications**. The 23 scenarios are enumerated in `simulation/cluster/dispatch_full.sh` (S01–S23) and vary in

| factor | levels |
|---|---|
| sample size `n` | 200, 500, 1000 |
| error distribution `epsilon` | `N(0, 9)` (`norm`) or centered min-extreme-value with SD 3 (`ev`) |
| link `G` | logistic (`10 * plogis(s, scale = 5)`) or exponential (`exp((s+1)/5)`) |
| target censoring rate | 15%, 30%, 40% |
| `beta_1` | 1 (mild U-shape) or 2 (steeper U-shape) |

For each scenario the following six methods are fit on the training set and evaluated on an independent test set of 5,000 observations:

- **MCE** — proposed rank-based Maximum C-index Estimator
- **Cox+spline** — Cox with a natural spline in `X`
- **Cox+quadratic** — Cox with `X` and `X^2`
- **Cox+quadratic with X:Z interactions**
- **Cox+spline with X:Z interactions**
- **Oracle** — Cox piecewise-linear in `X`, with the true change point known (revised specification computed by `reprocess_oracle_v2.R`)

MCE point estimation uses `DEoptim` on the hard C-index objective; SEs come from a multi-start BOBYQA bootstrap (525 draws, first 500 retained) on the hard C-index. Delta-method Wald and bootstrap-percentile CIs for the critical point are both saved.

Full model, DGM, and estimation details are given in Sections 2, 3, and 6 of the manuscript.

---

## Reproducing the simulation results

### Dependencies

- R >= 4.0 with packages `DEoptim`, `nloptr`, `survival`, `SurvMetrics`, `Rcpp`, `splines`, `dplyr`, `ggplot2`, `patchwork`
- A C++ toolchain compatible with `Rcpp::sourceCpp` (gcc 11+ or equivalent)
- SLURM (for the cluster dispatch scripts); a local run is also possible for a single scenario at reduced replicate counts

### Single-replicate smoke test (local)

From `simulation/`:

```
Rscript simulation_one_rep.R S07 1 1000 norm logistic 0.3000 8.3000 TRUE 2
```

Arguments are `scenario_id seed n ei G_type censor_lo censor_hi run_competitors [beta1]`. Output is written to `results/S07/rep_0001_j<batch_tag>.rds` as one R list containing the MCE estimates, bootstrap covariance, all five competitor results, and provenance fields (`scenario_id`, `seed`, `batch_tag`, `timestamp`, `n`, `ei`, `G_type`, `beta1`, censoring bounds).

### Full study on SLURM

From `simulation/`:

```
bash cluster/dispatch_full.sh
```

Submits 23 array jobs of 1050 replicates each (~24,150 tasks total). Environment overrides: `N_REPS`, `WALL`, `MEM`, `RESULTS_DIR`, `MYCINDEX_CPP_PATH`. After completion, run the oracle reprocessing patch that adds the corrected Cox oracle to each per-rep RDS:

```
bash cluster/dispatch_reprocess_oracle_v2.sh
```

### Aggregate to CSV

```
Rscript summarize_all_to_csv.R results csv_output
```

Produces `csv_output/{sim_meta, sim_params, sim_xc, sim_cindex, sim_survival, sim_cr, sim_boot_diag}.csv`. Each row aggregates the first 1000 successful replicates per scenario.

### Rebuild the manuscript figures

```
Rscript build_figs.R csv_output figures
```

Produces `figures/fig_sim_xc.pdf` and `figures/fig_sim_cindex.pdf`. Each is a three-panel layout with method colors and shapes shared across both files. The 23-scenario grouping used for panels (a)/(b)/(c) is embedded in `build_figs.R` and mirrors `cluster/dispatch_full.sh`.

---

## Output provenance

Every per-replicate RDS embeds `scenario_id`, `seed`, `batch_tag` (SLURM job id when available, otherwise a UTC timestamp), and the exact scenario parameters used to generate the record. Cross-check provenance across replicates before pooling or comparing results from different cluster runs.

---

## Real-data analysis

`real_data_analysis/real data_lilly.R` is the Eli Lilly dataset analysis prototype; `bootstrap_realdata_Lilly.R` is a bootstrap stub. The UK Biobank analysis reported in Section 7 of the current manuscript is not distributed in this repository because the data are subject to UK Biobank access terms.

---

## Contact

Corresponding author: Yuxin (Daisy) Zhu, Department of Biostatistics, Johns Hopkins Bloomberg School of Public Health (`daisy@jhu.edu`).

Co-authors: Zhirui Fu, Mei-Cheng Wang, Yu Du.

## License

TBD.
