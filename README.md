# MIDFA: Scalable Bayesian Factor Analysis for Mixed and Incomplete Data

MIDFA is sparse Bayesian factor analysis for data with mixed variable types and missing values. It learns the number of factors from the data.

**The method is `mcem_fa_algorithm.R` plus the C++ code in `src/`. Nothing else is needed to use it.** The rest of the repository reproduces the results of the paper, and its NO.MS analyses need data that is not public.

## Using MIDFA on your own data

You need R (≥ 4.3) with a C++ compiler and the packages `Rcpp`, `RcppArmadillo`, `RcppTN`, `nloptr` and `invgamma`.

```r
source("path/to/mcem_fa_algorithm.R")   # compiles the C++ in src/
Y <- as.matrix(my_data)                 # N x D numeric, NA = missing
set.seed(1)
fit <- dpe_func(v0s = 10^-(1:4), Y = Y, K = 10, px_sigma = FALSE, px_rotate = FALSE,
                varimax_every = -1, conv_param = TRUE, tol = 0.016,
                nburn = 100, n_mcmc = 100, verbose = -1)
final <- fit[[length(fit)]]             # the last, smallest-v0 stage
loadings <- final$lambda_corr * (final$w > 0.5)
loadings <- loadings[, colSums(final$w > 0.5) > 0, drop = FALSE]
```

- **`Y`:** variable types are detected. A column with two distinct values is treated as binary, and any other column enters through its ranks, so ordinal, count and continuous variables need no transformation.
- **`K`:** an upper bound on the number of factors. Factors the data do not support are switched off.
- **`v0s`:** the spike variances, from large to small. Each stage is warm-started from the one before.
- **`final$w`:** the posterior inclusion probabilities. A loading is selected when its probability exceeds 0.5.
- **`final$lambda_corr`:** the loadings on the correlation scale.

The call uses the settings of simulation scenarios 1-3. The other arguments are those of `mcem_algorithm`, which `dpe_func` passes on, and their defaults are in its definition.

## Citation

If you use this code, please cite the manuscript:

> Hutchings, G., Samartsidis, P., Donnay, C., Gaetano, L., Fisher, E., Nichols, T. E., Holmes, C., Häring, D. A., & Ganjgahi, H. (2026). *MIDFA: Scalable Bayesian Factor Analysis for Mixed and Incomplete Data*.

For questions, contact **habib.ganjgahi@bdi.ox.ac.uk**.

---

## Paper Results Guide

Each table and figure in the paper and supplement, with the code that produces it.

| Paper Item | Description | Script(s) | Data Required | Cluster | Output |
|---|---|---|---|---|---|
| **Figure 1** | DPE regularisation path | `simulations/dpe_path.R` | Synthetic | Local (about 1 min) | `simulations/dpe_maintext_pca_r1.pdf` |
| **Figure 2** | True Λ (12×3, scenarios 1-2) | `simulations/create_true_loading_matrix.R` | None | Local | `true_factor_loading_matrix.pdf` |
| **Figure 3** | True Λ (8×3, scenario 3) | `simulations/create_true_loading_matrix_mixed.R` | None | Local | `true_loading_mixed_missing.pdf` |
| **Table 1** | Continuous simulation | `simulations/continuous/run.R` → `simulations/analyse.R` | Synthetic | BMRC | Table printed by `analyse.R` |
| **Table 2** | Binary simulation | `simulations/binary/run.R` → `simulations/analyse.R` | Synthetic | BMRC | Table printed by `analyse.R` |
| **Table 3** | Mixed/missing simulation | `simulations/mixed/run.R` → `simulations/analyse.R` | Synthetic + `missingness_analysis.rds`† | Novartis | Table printed by `analyse.R` |
| **Table 4** | Clinical loading matrix (Λ̃) | `noms-analysis/clinical/run.R` → `noms-analysis/clinical/analyse.R` | NO.MS† | Novartis | Run 1 of the `analyse.R` output |
| **Figure 4** | MoCo spatial maps (100 MoCos) | `noms-analysis/mri/analyse_interactive.R` | NO.MS† | Novartis | `factor_atlas_*.nii.gz` → FSLeyes |
| **Figure 5** | PCA first principal component | `noms-analysis/mri/make_pca_supp_figure.sh` (runs `make_pca_loading_maps.R`) | NO.MS† | Novartis | `pca1.png` |
| **Figure 6** | CDW-significant MoCos | `noms-analysis/mri/clinical_analysis/analyse_cox_survival.R` | NO.MS† | Novartis | `coxph_sig_factors_*.nii.gz` → FSLeyes |
| **Cox p-values** (§3.3) | MoCo 48, MoCo 76, VOLT2 | `noms-analysis/mri/clinical_analysis/analyse_cox_survival.R` | NO.MS† | Novartis | Job log and `cox_survival_all.csv` |
| **Figure S1** | Per-element loading MSE, scenario 1 | `simulations/plot_elementwise_mse.R` on `continuous/elementwise_mse.rds` | Synthetic | Any | `mse_continuous.pdf` |
| **Figure S2** | True Λ (1,956×5, high-dimensional) | `simulations/continuous_highdim/plot_true_loading.R` | None | Local | `true_loading_highdim.pdf` |
| **Table S1** + **Figure S3** | High-dimensional study | `simulations/continuous_highdim/run.R` → `simulations/analyse.R`, then `simulations/continuous_highdim/plot_elementwise_mse.R` | Synthetic | BMRC | Table printed by `analyse.R`, `mse_continuous_highdim.pdf` |
| **Table S2** + **Figure S4** | Sample-size study | `simulations/continuous_samplesize/run.R` → `analyse.R` → `plot_elementwise.R` | Synthetic | BMRC | Table printed by `analyse.R`, `table_rank.csv`, `elementwise_mse_rank.pdf` |
| **Table S3** | E-step sample-size study | `simulations/continuous_estep/run.R` → `analyse.R` | Synthetic | BMRC | Table printed by `analyse.R` |
| **Figure S5** | Per-element loading MSE, scenario 2 | `simulations/plot_elementwise_mse.R` on `binary/elementwise_mse.rds` | Synthetic | Any | `mse_binary.pdf` |
| **Figure S6** | Margins of the NO.MS variables behind scenario 3 | `simulations/mixed/plot_margins.R` | NO.MS† | Novartis | `margins_mixed.pdf` |
| **Figure S7** | Per-element loading MSE, scenario 3 | `simulations/plot_elementwise_mse.R` on `mixed/elementwise_mse.rds` | Synthetic + `missingness_analysis.rds`† | Any | `mse_mixed.pdf` |
| **Figures S8, S9** | PCA second and third principal components | `noms-analysis/mri/make_pca_supp_figure.sh` | NO.MS† | Novartis | `pca2.png`, `pca3.png` |

†Requires NO.MS data (see Data Availability).

### Mapping the `analyse.R` tables to the paper

`analyse.R` prints one row per method (per N or per M in the sample-size and E-step studies, which have their own `analyse.R`). Each entry is a mean with its standard error in brackets, and `n` is the number of realisations.

- **Table 1** (continuous): rows `our method` (MIDFA), `factanal K known`, `pca K known`, `rockova 2016`.
- **Table 2** (binary): rows `our method`, `nmf K known`, `li 2023 mirt`.
- **Table 3** (mixed): rows `full data`, `missing data`, `complete cases`.
- **Table S1** (high-dimensional): rows `our method`, `pca K known`, `rockova 2016`.
- **Table S2** (sample size): one row per N. Also written to `table_rank.csv` by `simulations/continuous_samplesize/analyse.R`.
- **Table S3** (E-step): one row per number of E-step samples M.

Columns: `struct_rec` is the Struct. Rec. Rate, `K_correct` the Correct K Rate, `false_pos` the Avg. FP, `false_neg` the Avg. FN and `lambda_mse` the Λ MSE. `mean_time` is not reported in the paper.

---

## Data Availability

The Novartis-Oxford Multiple Sclerosis (NO.MS) data are not public. Anonymised raw data and related documents can be requested through CSDR (https://www.clinicalstudydatarequest.com) under a data-sharing agreement with Novartis, with requests reviewed and approved by an independent CSDR panel.

Items marked † in the Paper Results Guide need NO.MS data:
- scenario 3 (Table 3, Figure S7), through `missingness_analysis.rds`, which `simulations/mixed/missingness-predictor/missingness-clinical.R` generates from the NO.MS clinical data;
- Figure S6;
- the clinical analysis (Table 4);
- the MRI analyses (Figures 4-6, S8, S9 and the Cox p-values).

Everything else runs on synthetic data alone.

The code reads NO.MS files from these Novartis paths. Edit them for your copy:
- Clinical extract `03_baseline_results.csv`: `/data/ms/processed/mri/MS_Share/4George/` (`noms-analysis/clinical/run.R`, `simulations/mixed/plot_margins.R`) and `/data/ms/processed/clinical/projects/fahmm/interim_tables/` (`missingness-clinical.R`). The clinical analysis also reads `sdmt.sas7bdat` from `/data/ms/unprocessed/clinical/NOMS_Version2_20211022-OXF-ANALYTICS/`.
- MRI: the subjects × voxels lesion matrix `processed_fa_scans.rds` under `/data/users/uu85g9/factor_analysis-oldold/consolodation_may-24/` (rebuilt from the individual lesion maps by `noms-analysis/mri/build_mri_data.sh`, whose header explains its output path) and the brain mask `final_mask.nii.gz` under `/data/ms/processed/mri/MS_Share/4George/`. The Cox analysis also reads `longitudinal.csv` and `T2lesion_ind/` under `/data/ms/processed/mri/MS_Share/4George/`, `longitudinal-table-CORE.csv` under `/data/ms/processed/clinical/projects/longitudinal/`, and `mri-clinical_processed.rds` under `/data/users/uu85g9/factor-analysis-old/mri-clinical/`.

---

## Prerequisites

- **R** (≥ 4.3) with a C++ compiler and the packages listed under *Using MIDFA on your own data*
- **SLURM**, on two clusters, BMRC and Novartis. The submit scripts load `R-bundle-CRAN/2023.12-foss-2023a` on BMRC and `george-bundle/0.1-foss-2023a-R-4.3.2` on Novartis
- Summaries and plots: `dplyr` and `data.table` (`analyse.R`), `ggplot2` (≥ 3.4), `reshape2`, `patchwork`, `scales`, `pheatmap`
- Baseline methods (simulation competitors):
  - NMF: `NMF` (PCA and factanal use base R)
  - Ročková & George: `MCMCpack`, `partitions`, `mvtnorm`, `nloptr`, `glmnet`
  - Li et al.: `pacman`, `tidyverse`, `glmnet`, `TruncatedNormal`, `progress`, `hash`, `nloptr`, `furrr`, `purrr`, `plot.matrix`, `argparser`, `feather`, `yaml`, `here`. `pacman::p_load` installs missing packages at run time, so preinstall them on clusters without internet access
- Clinical analysis: `tidyverse`, `haven`
- MRI analysis: `RNifti`, `survival`, `data.table`, `dplyr`, `tidyr`, `ggplot2`, `gridExtra`
- Figure generation: [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki) (`fslmerge`, [FSLeyes](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLeyes)) and ImageMagick (MRI figures)

Clone the repository to `~/factor-analysis` on each cluster: every `source()` path starts there (`/home/uu85g9/factor-analysis` in the clinical and MRI scripts). An existing directory of that name would be sourced instead, so move it first or edit the paths. Also edit the hard-coded output paths (`/data/users/uu85g9/...`, `/well/nichols-nvs/...`) and the SLURM headers (`--constraint`, `--partition` and, in `simulations/binary/submit-run.sh`, `--account`).

---

## Repository Layout

| Path | Description | Paper Reference |
|---|---|---|
| `mcem_fa_algorithm.R` | Core MCEM algorithm with `dpe_func()` for deterministic path expansion | Algorithms 1-2, §2.1 |
| `src/` | C++ E-step (`mc_e_step.cpp`), truncated normal sampler (`rtn1.*`), RNG (`rng.*`), utilities (`utils.*`) | §2.1 |
| `simulations/continuous/` | Scenario 1: Gaussian data, 12 vars, 3 factors | Table 1, Figure S1 |
| `simulations/binary/` | Scenario 2: binary data, 12 vars, 3 factors | Table 2, Figure S5 |
| `simulations/mixed/` | Scenario 3: mixed types + structured missingness, 8 vars, 3 factors. `missingness-predictor/` fits the missingness models and marginal ECDFs from NO.MS and generates the data, `plot_margins.R` plots the margins | Table 3, Figures S6, S7 |
| `simulations/continuous_highdim/` | Scenario 1 at 1,956 variables and N=100 | Table S1, Figures S2, S3 |
| `simulations/continuous_samplesize/` | Scenario 1 at six sample sizes | Table S2, Figure S4 |
| `simulations/continuous_estep/` | Scenario 1 with the number of E-step samples varied | Table S3 |
| `simulations/methods/` | Baseline implementations: factanal, PCA+varimax, NMF, Li et al. MIRT, Ročková & George 2016 | Tables 1-3, S1 |
| `simulations/shared/` | Evaluation metrics: structure recovery, correct K, false positives and negatives, loading MSE (`evaluate.R`). `evaluate_highdim.R` matches estimated to true factors by support overlap | Tables 1-3, S1-S3 |
| `simulations/analyse.R` | Aggregates SLURM array outputs into summary tables and writes the per-element MSE | Tables 1-3, S1 |
| `simulations/dpe_path.R`, `create_true_loading_matrix*.R`, `plot_elementwise_mse.R` | Figure scripts | Figures 1-3, S1, S5, S7 |
| `noms-analysis/clinical/` | 10-seed clinical factor analysis (9 variables, N=8,023) | Table 4, §3.2 |
| `noms-analysis/mri/` | MRI lesion mask analysis → Modes of lesion Co-occurrence (MoCos), and the PCA maps | §3.3, Figures 4, 5, S8, S9 |
| `noms-analysis/mri/clinical_analysis/` | Cox PH survival models | Figure 6, §3.3 |
| `paper-results/` | Pre-computed NIfTI atlases, Cox CSV, rendered figures (not tracked, see retrieval instructions) | Figures 4-6 |

---

## Simulation Settings

### Scenarios 1-3

Settings shared by the three scenarios:

| Parameter | Value |
|-----------|-------|
| v0 (DPE schedule) | 0.1, 0.01, 0.001, 0.0001 |
| v1 | 10 |
| IBP alpha | 2 |
| K (truncation) | 10 |
| N | 200 |
| Realisations | 500 |
| Parameter expansion | None (`px_sigma = FALSE`, `px_rotate = FALSE`) |
| Varimax | Applied to loading matrix between DPE stages (not within) |
| Convergence | Max absolute loading change < 0.016 (at most 2,000 iterations per stage) |
| E-step | Burn-in of 100 sweeps, run once at the start of each DPE stage. Each EM iteration then averages M = 100 samples |
| Seeds | Data: `12345 + r` for realisation r (scenario 3: `1234 + r`). Fit: 12345 |

- **Scenario 1:** Λ is 12×3 (variables 1-6, 7-10 and 11-12 load on one factor each, with loading 1), and Y = Z has Gaussian margins.
- **Scenario 2:** the same Λ, with all 12 variables made binary through the Gaussian copula (success probabilities from 0.25 to 0.75).
- **Scenario 3:** Λ is 8×3 (variables 1-4, 4-6 and 7-8, so variable 4 loads on two factors). Margins follow the NO.MS empirical CDFs, and missingness is drawn from models fitted to NO.MS. The conditions are `full data`, `missing data` and `complete cases`.
- **Comparators:** factanal, PCA+varimax and NMF are given the true K=3. MIDFA, Ročková & George and Li et al. use K=10. The Ročková & George comparator (`simulations/methods/veronika2016_wrapper.R`) runs lambda0 = 5, 10, 20, 30, 40 (lambda1 = 0.001, epsilon = 0.05, alpha = 1/D), each stage warm-started from the previous one as in the authors' `EXAMPLE_FACTOR_ROTATE.R`. In scenario 1 it is not varimax-rotated (MIDFA is not rotated within stages either). In the high-dimensional study it is rotated after every EM iteration (`varimax_all = TRUE`), since the lambda0 = 20 stage does not converge without it.

### Higher-dimensional, sample-size and E-step studies

These repeat scenario 1 with the changes below. Everything else, including the seeds, is as above.

| Setting | High-dimensional (Table S1) | Sample size (Table S2) | E-step (Table S3) |
|---|---|---|---|
| Data | D=1,956, N=100, true K=5 (blocks of 500 variables, adjacent blocks overlap in 136), Y = Z | D=12, true K=3 as scenario 1, N = 50, 250, 500, 1000, 2500, 5000 | As scenario 1 |
| K (truncation) | 20 (PCA given K=5) | 10 | 10 |
| v0 (DPE schedule) | 0.1, 0.075, 0.05 | 0.05, 0.03, 0.015, 0.005 | As scenario 1 |
| Iterations per stage | At most 100 | 100 | As scenario 1 |
| Convergence tolerance | 0.045 | 1e-6 | 0.016 |
| Parameter expansion | `px_sigma` and `px_rotate` | `px_sigma` only | None |
| E-step samples M | 100 | 100 | 50, 100, 200, 300, 400 |
| Burn-in | 100 | 100 | 100 |
| Methods | MIDFA, Ročková & George, PCA | MIDFA | MIDFA |

---

## Simulations

### Before the first run

- Submit each job from its study directory, for example `cd ~/factor-analysis/simulations/continuous`. Logs go to `logs/` there, which SLURM needs to exist already, so each simulation directory (and `noms-analysis/clinical/`) ships an empty `logs/.gitkeep`.
- Compile the C++ code once before submitting an array, so the tasks do not all build `src/*.o` at once. On a compute node, with the cluster's module loaded:
  ```bash
  Rscript -e 'source("~/factor-analysis/mcem_fa_algorithm.R")'
  ```
  The algorithm finds `src/` from its own path, so source it at the top level of a script, not inside a function.
- Create the output directories of `continuous`, `binary` and `mixed` first (the other studies create theirs):
  ```bash
  mkdir -p /well/nichols-nvs/users/peo100/factor-analysis/paper-code/simulations/{continuous,binary}      # BMRC
  mkdir -p /data/users/uu85g9/factor-analysis/paper-code/simulations/mixed                                  # Novartis
  ```
- Start with empty output directories: `analyse.R` reads every `<number>.rds` file in its directory.

### Studies

| Study | Directory (under `simulations/`) | Cluster | `--array` | Tasks | Time, memory | Task id → (method or setting, realisation r) |
|---|---|---|---|---|---|---|
| Scenario 1 (Table 1) | `continuous/` | BMRC | `0-1999` | 2,000 | 1 h, 2 GB | `id %% 4` = 0 `our method`, 1 `rockova 2016`, 2 `factanal K known`, 3 `pca K known`. `r = id %/% 4 + 1` |
| Scenario 2 (Table 2) | `binary/` | BMRC | `0-1499` | 1,500 | 1 h, 8 GB | `id %% 3` = 0 `our method`, 1 `nmf K known`, 2 `li 2023 mirt`. `r = id %/% 3 + 1` |
| Scenario 3 (Table 3) | `mixed/` | Novartis | `0-1499` | 1,500 | 2 h, 2 GB | `id %% 3` = 0 `full data`, 1 `missing data`, 2 `complete cases`. `r = id %/% 3 + 1` |
| High-dimensional (Table S1) | `continuous_highdim/` | BMRC | `0-1499` | 1,500 | 12 h, 8 GB | `id %% 3` = 0 `our method`, 1 `rockova 2016`, 2 `pca K known`. `r = id %/% 3 + 1` |
| Sample size (Table S2) | `continuous_samplesize/` | BMRC | `0-2999` | 3,000 | 3 h, 4 GB | `id %% 6` = 0-5 gives N = 50, 250, 500, 1000, 2500, 5000. `r = id %/% 6 + 1` |
| E-step (Table S3) | `continuous_estep/` | BMRC | `0-2499` | 2,500 | 10 h, 2 GB | `id %% 5` = 0-4 gives M = 50, 100, 200, 300, 400. `r = id %/% 5 + 1` |

The high-dimensional (about 190 core-hours) and sample-size (about 300 core-hours) studies are the most expensive.

### Running a study

Chain the analysis job with `afterany`, not `afterok`, which a single failed task would leave pending forever. Missing tasks show up in `n` and can be resubmitted with `sbatch --array=<ids> submit-run.sh`.

```bash
# BMRC
for study in continuous binary continuous_highdim continuous_samplesize continuous_estep
do
  cd ~/factor-analysis/simulations/$study
  jid=$(sbatch --parsable submit-run.sh)
  sbatch --dependency=afterany:$jid submit-analyse.sh
done
```

Each `submit-analyse.sh` hard-codes the directory it reads (`BASE` in `continuous_samplesize`), so edit it if you change where the run writes.

Scenario 3 (NO.MS data, Novartis) first fits the missingness model:

```bash
cd ~/factor-analysis/simulations/mixed
Rscript missingness-predictor/missingness-clinical.R     # writes missingness_analysis.rds

jid=$(sbatch --parsable submit-run.sh)
sbatch --dependency=afterany:$jid submit-analyse.sh
```

If you change the path of `missingness_analysis.rds`, edit `save_path` in `missingness-clinical.R` and `missingness_rds_path` in `simulate_missing_data_functions.R` to match (`run.R` reads the latter).

#### Outputs

`<BMRC>` is `/well/nichols-nvs/users/peo100/factor-analysis/paper-code/simulations`. The summary table is printed to the log of the analysis job in `logs/`.

| Study | Output directory | Written by the analysis step |
|---|---|---|
| Scenario 1, scenario 2, high-dimensional | `<BMRC>/continuous/`, `<BMRC>/binary/`, `<BMRC>/continuous_highdim/` | `elementwise_mse.rds` and `all_results.rds` |
| Scenario 3 | `/data/users/uu85g9/factor-analysis/paper-code/simulations/mixed/` (also holds `missingness_analysis.rds`) | `elementwise_mse.rds` and `all_results.rds` |
| Sample size | `<BMRC>/continuous_samplesize/n<N>/` for each N | in `<BMRC>/continuous_samplesize/`: `table_rank.csv`, `failures.csv`, `all_results.rds`, `all_runs.rds`, `elementwise_mse.rds`, `elementwise_mse_rank.pdf` |
| E-step | `<BMRC>/continuous_estep/` | `all_results.rds` |

#### Running a single task

To run one task outside SLURM:

```bash
cd ~/factor-analysis/simulations/continuous_estep
SLURM_ARRAY_TASK_ID=1 SIM_OUT=$HOME/test_out/ Rscript run.R
```

Output goes to `save_dir` in `run.R`. `SIM_OUT` (ending in `/`) overrides it in `continuous_highdim`, `continuous_samplesize` and `continuous_estep`. In the others, edit `save_dir` and create the directory. Use a scratch directory for a test, as `analyse.R` reads every result file in its directory. Use one R process per task: forked processes (`mclapply`) collide, because the algorithm compiles inline C++ on each call.

#### Reproducibility

Each task is reproducible on its own: the data seed depends only on the realisation, and the fit seed is constant. MIDFA fits repeat exactly on the same machine but not across machines, even on the BMRC skylake nodes and module used for the paper. The C++ code is compiled with the local, possibly CPU-specific, compiler flags, and a last-digit floating-point difference makes the Monte Carlo draws diverge. Individual realisations and iteration counts can therefore differ, and summary-table entries agree to within about two standard errors. Li et al. fits also vary between runs. PCA, factanal, NMF and Ročková & George reproduce exactly.

### Figures 1-3 and the simulation figures

```bash
cd ~/factor-analysis/simulations
Rscript dpe_path.R                           # Figure 1: fits once (cached in dpe_maintext_pca_r1.rds), then plots
Rscript create_true_loading_matrix.R         # Figure 2: true_factor_loading_matrix.pdf
Rscript create_true_loading_matrix_mixed.R   # Figure 3: true_loading_mixed_missing.pdf
```

`dpe_path.R` writes next to the script, and the other two write to the current directory. The per-element MSE figures read the analysis step's `elementwise_mse.rds` (copy it off the cluster to plot locally). Always pass the methods list: it sets the panel order and the panels that fix the shared colour scale.

```bash
cd ~/factor-analysis/simulations
BMRC=/well/nichols-nvs/users/peo100/factor-analysis/paper-code/simulations
MIXED=/data/users/uu85g9/factor-analysis/paper-code/simulations/mixed

# Figure S1, S5, S7
Rscript plot_elementwise_mse.R $BMRC/continuous/elementwise_mse.rds mse_continuous.pdf 'our method,factanal K known,pca K known,rockova 2016'
Rscript plot_elementwise_mse.R $BMRC/binary/elementwise_mse.rds mse_binary.pdf 'our method,nmf K known,li 2023 mirt'
Rscript plot_elementwise_mse.R $MIXED/elementwise_mse.rds mse_mixed.pdf 'full data,missing data,complete cases'

# Figure S2, S3 (high-dimensional)
Rscript continuous_highdim/plot_true_loading.R true_loading_highdim.pdf
Rscript continuous_highdim/plot_elementwise_mse.R $BMRC/continuous_highdim/elementwise_mse.rds mse_continuous_highdim.pdf 'our method,rockova 2016,pca K known'

# Figure S6 (NO.MS, Novartis): writes margins_mixed.pdf and prints N and the missing fraction of each variable
Rscript mixed/plot_margins.R                 # optional arguments: [csv] [out_dir]
```

Figure S4 (`elementwise_mse_rank.pdf`) is written by `continuous_samplesize/submit-analyse.sh`, which runs `plot_elementwise.R` after `analyse.R`, and can be rerun on its own with `Rscript continuous_samplesize/plot_elementwise.R <output directory>/`.

---

## Clinical Analysis (Table 4)

Uses the NO.MS clinical data. Recovers K=5 latent dimensions of MS disease status from 9 clinical variables (N=8,023, with PASAT missing for 1,711 and SDMT for 4,105 participants, which the model handles directly).

**Variables:** EDSS, T25FWM, HPT9M, PASAT, SDMT, VOLT2, NBV2, NUMGDT1, RELAPSE.

**Settings:** DPE v0s = {1, 0.8, 0.4, 0.2, 0.1, 0.05}, K=10, v1=10, 101 EM iterations per stage (the tolerance of 1e-7 is never reached, so every stage runs all 101), 100 burn-in sweeps and 100 MC samples, PX rotation + PX sigma, ibp_alpha=2. The array index is the random initialisation (seed 1234 + index).

```bash
cd ~/factor-analysis/noms-analysis/clinical
sbatch submit.sh       # array 1-10: 10 random initialisations, about 13 minutes each

# After jobs complete:
sbatch --dependency=afterany:<JOBID> --partition=general,legacy \
  --wrap="module load george-bundle/0.1-foss-2023a-R-4.3.2 && Rscript ~/factor-analysis/noms-analysis/clinical/analyse.R"
```

Compile the C++ E-step once before submitting, as described under Simulations. The job log prints N and the number of missing values per variable. A single task can be run with `SLURM_ARRAY_TASK_ID=1 Rscript run.R`.

Results are saved to `/data/users/uu85g9/factor-analysis/paper-code/noms/clinical/results_{r}.rds`.

`analyse.R` prints, for each of the 10 runs, the number of factors found, the total iterations and the final Q function value, followed by the run's loading matrix. The loadings are those of the last DPE stage, with entries of posterior inclusion probability below 0.5 set to zero, all-zero columns dropped, and normalised as in the paper's Eq. (6), with rows and columns ordered as in Table 4. **Table 4 is the matrix of Run 1.**

All 10 runs find the same five dimensions:

1. **Physical disability** - EDSS, T25FWM, HPT9M
2. **Cognitive function** - PASAT, SDMT
3. **Asymptomatic activity** - VOLT2, NUMGDT1
4. **Brain damage** - VOLT2, NBV2 (opposite signs)
5. **Relapse** - RELAPSE

The main difference between runs is a negative SDMT loading on physical disability, which is present in Run 1 and in 7 of the 10 runs behind the paper.

---

## MRI Analysis (Figures 4-6, S8-S9; Cox p-values)

NO.MS binary lesion masks, D=51,595 voxels, N=2,093.

**Settings:** K=100, v0 ∈ {0.1, 0.01, 0.001}, v1=10, 50 samples per E-step, convergence tolerance 0.01 (first stage), PX sigma.

The paper's MRI fits ran on an earlier copy of the algorithm. Its E-step, M-step and parameter expansion are the same as in `mcem_fa_algorithm.R`, but it stored the log prior instead of the Q function, so a fresh run reports different convergence diagnostics.

### Step 1: Primary run

```bash
cd ~/factor-analysis/noms-analysis/mri/dpe_random_pxsigma
mkdir -p logs
sbatch submit.sh
```

### Step 2: Restart from DPE stage 3 (paper model)

The primary run is restarted at DPE stage 3 (v0=0.001) from iteration 70, where the parameter differences had plateaued.

```bash
cd ~/factor-analysis/noms-analysis/mri/dpe_random_pxsigma_restart_dpe3
mkdir -p logs
sbatch submit.sh
```

### Step 3: Generate MoCo atlas (Figure 4)

```bash
# Set condition and base_path in the script, then run interactively:
Rscript ~/factor-analysis/noms-analysis/mri/analyse_interactive.R
```

This creates `factor_atlas_*.nii.gz`, visualised in FSLeyes (see Figure Generation below). The paper's atlas is the restart run's DPE stage 1 at iteration 70 (`factor_atlas_dpe_random_pxsigma_restart_dpe3_iter70_dpe1.nii.gz`), so set `iter <- 70` and `dpe_stage <- 1` at the top of the script. The defaults of -1 read the final stage.

### Step 4: Cox PH survival models (Figure 6, Cox p-values)

Compresses subjects into MoCo space and fits Cox models predicting 3-month confirmed disability worsening (M3CDW).

```bash
cd ~/factor-analysis/noms-analysis/mri/clinical_analysis
mkdir -p logs compressed_data results

Rscript build_clinical_analysis_inputs.R               # once: writes inputs/subjects.rds and inputs/Y.rds

jid=$(sbatch --parsable --array=2 submit_compress.sh)  # compress subjects into MoCo space
sbatch --dependency=afterany:$jid submit_cox.sh        # Cox PH models: p-values in logs/cox_surv_<jobid>.out
```

Only compression task 2 (`restart_dpe3_dpe1_iter70`) feeds the paper. Tasks 1 and 3 are other checkpoints, and task 4 needs a `dpe_pca_nopx` run that is not in this repository.

Paper reports: MoCo 48 (p=0.0058), MoCo 76 (p=0.0083). T2 lesion volume is not significant (p=0.5554).

### Step 5: PCA loading maps (Figure 5, Figures S8-S9)

Maps the first three principal components of the lesion data, for comparison with the MoCos. Run on a compute node, with FSL on the `PATH` and ImageMagick available (`module load ImageMagick` if needed):

```bash
bash ~/factor-analysis/noms-analysis/mri/make_pca_supp_figure.sh
```

`make_pca_loading_maps.R` computes the PCA of the standardised lesion matrix through its 2,093 × 2,093 Gram matrix (no random numbers, so reruns are identical), signs each component so that its largest loading is positive, and writes `pca_loadings.nii.gz`. The shell script then renders each component over the MNI152 1 mm template with `fsleyes` (unthresholded axial lightbox, red-yellow for positive and blue for negative loadings) and autocrops the images. Outputs are written to `/data/users/uu85g9/factor-analysis/noms/mri/pca_supp/`: `pca1.png` (Figure 5), `pca2.png` (Figure S8), `pca3.png` (Figure S9), `pca_loadings.nii.gz` and `pca_variance_explained.txt`.

---

## Retrieving Paper Results from HPC

`paper-results/` is not tracked. To populate it:

```bash
# MoCo atlas NIfTI (Figure 4)
scp uu85g9@ms-login-00:/data/users/uu85g9/factor-analysis/noms/mri/runs_jan18_v0_0.1/dpe_random_pxsigma_restart_dpe3/factor_atlas_dpe_random_pxsigma_restart_dpe3_iter70_dpe1.nii.gz ~/factor-analysis/paper-results/

# Significant Cox factors NIfTI (Figure 6)
scp uu85g9@ms-login-00:/data/users/uu85g9/factor-analysis/noms/mri/runs_jan18_v0_0.1/clinical_analysis/results/coxph_sig_factors_restart_dpe3_dpe1_iter70.nii.gz ~/factor-analysis/paper-results/

# Cox results table
scp uu85g9@ms-login-00:/data/users/uu85g9/factor-analysis/noms/mri/runs_jan18_v0_0.1/clinical_analysis/results/cox_survival_all.csv ~/factor-analysis/paper-results/

# PCA maps (Figures 5, S8, S9)
scp "uu85g9@ms-login-00:/data/users/uu85g9/factor-analysis/noms/mri/pca_supp/pca*.png" ~/factor-analysis/paper-results/
```

These are the locations of the results used in the paper. The MRI scripts write to the same paths without the `runs_jan18_v0_0.1/` level (for example `/data/users/uu85g9/factor-analysis/noms/mri/dpe_random_pxsigma_restart_dpe3/`).

### Extracting Cox p-values from CSV

```r
library(dplyr)
cox <- read.csv("paper-results/cox_survival_all.csv")

cox %>%
  filter(dataset == "restart_dpe3_dpe1_iter70", grepl("^V", variable), p_value < 0.05) %>%
  select(variable, coef, exp_coef, se, p_value) %>%
  arrange(p_value)
```

---

## Figure Generation

Figures 4 and 6 are generated from NIfTI files using [FSLeyes](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLeyes). After retrieving the NIfTI files (see above), set `PAPER_RESULTS` to your local `paper-results/` directory, for example `PAPER_RESULTS=~/factor-analysis/paper-results`, and check the MNI template path (`/usr/local/fsl/data/standard/`) matches your FSL install:

**Figure 4** - MoCo atlas (`moco_atlas.png`):
```bash
fsleyes --scene lightbox --worldLoc 16.863185978095004 -18.000099182128906 30.609235989629497 --displaySpace $PAPER_RESULTS/coxph_sig_factors_restart_dpe3_dpe1_iter70.nii.gz --zaxis 2 --zrange 0.3199999999068677 0.6800000000931322 --sliceSpacing 0.029945054973519132 --sampleSlices centre --labelSpace none --sliceOverlap 0.0 --hideCursor --cursorWidth 1.0 --bgColour 0.0 0.0 0.0 --fgColour 1.0 1.0 1.0 --cursorColour 0.0 1.0 0.0 --colourBarLocation top --colourBarLabelSide top-left --colourBarSize 100.0 --labelSize 12 --performance 3 --movieSync /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz --name "MNI152_T1_1mm_brain" --overlayType volume --alpha 100.0 --brightness 49.75000000000001 --contrast 49.90029860765409 --cmap greyscale --negativeCmap greyscale --displayRange 0.0 8447.64 --clippingRange 0.0 8447.64 --modulateRange 0.0 8364.0 --gamma 0.0 --cmapResolution 256 --interpolation none --numSteps 150 --blendFactor 0.1 --smoothing 0 --resolution 100 --numInnerSteps 10 --clipMode intersection --volume 0 $PAPER_RESULTS/coxph_sig_factors_restart_dpe3_dpe1_iter70.nii.gz --name "coxph_sig_factors_restart_dpe3_dpe1_iter70" --disabled --overlayType volume --alpha 100.0 --brightness 49.74999999999999 --contrast 49.90029860765409 --cmap greyscale --negativeCmap greyscale --displayRange 0.0 76.76 --clippingRange 0.0 76.76 --modulateRange 0.0 76.0 --gamma 0.0 --cmapResolution 256 --interpolation none --numSteps 150 --blendFactor 0.1 --smoothing 0 --resolution 100 --numInnerSteps 10 --clipMode intersection --volume 0 $PAPER_RESULTS/factor_atlas_dpe_random_pxsigma_restart_dpe3_iter70_dpe1.nii.gz --name "factor_atlas_dpe_random_pxsigma_restart_dpe3_iter70_dpe1" --overlayType volume --alpha 100.0 --brightness 49.75000000000001 --contrast 49.90029860765409 --cmap random --negativeCmap greyscale --displayRange 0.0 101.0 --clippingRange 0.0 101.0 --modulateRange 0.0 100.0 --gamma 0.0 --cmapResolution 256 --interpolation none --numSteps 150 --blendFactor 0.1 --smoothing 0 --resolution 100 --numInnerSteps 10 --clipMode intersection --volume 0
```

**Figure 6** - CDW-significant MoCos (`moco_cdw_significant.png`):
```bash
fsleyes --scene lightbox --worldLoc 16.863185978095004 -18.000099182128906 30.609235989629497 --displaySpace $PAPER_RESULTS/coxph_sig_factors_restart_dpe3_dpe1_iter70.nii.gz --zaxis 2 --zrange 0.44000000004656614 0.5933333334031825 --sliceSpacing 0.016904761899069067 --sampleSlices centre --labelSpace none --sliceOverlap 0.0 --hideCursor --cursorWidth 1.0 --bgColour 0.0 0.0 0.0 --fgColour 1.0 1.0 1.0 --cursorColour 0.0 1.0 0.0 --colourBarLocation top --colourBarLabelSide top-left --colourBarSize 100.0 --labelSize 12 --performance 3 --movieSync /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz --name "MNI152_T1_1mm_brain" --overlayType volume --alpha 100.0 --cmap greyscale --negativeCmap greyscale --useNegativeCmap --displayRange 0.0 8447.64 --clippingRange 0.0 8447.64 --modulateRange 0.0 8364.0 --gamma 0.0 --cmapResolution 256 --interpolation none --numSteps 150 --blendFactor 0.1 --smoothing 0 --resolution 100 --numInnerSteps 10 --clipMode intersection --volume 0 $PAPER_RESULTS/coxph_sig_factors_restart_dpe3_dpe1_iter70.nii.gz --name "coxph_sig_factors_restart_dpe3_dpe1_iter70" --overlayType volume --alpha 100.0 --cmap fsleyes_hsv --negativeCmap greyscale --useNegativeCmap --displayRange 0.0 76.76 --clippingRange 0.0 76.76 --modulateRange 0.0 76.0 --gamma 0.0 --cmapResolution 256 --interpolation none --numSteps 150 --blendFactor 0.1 --smoothing 0 --resolution 100 --numInnerSteps 10 --clipMode intersection --volume 0 $PAPER_RESULTS/factor_atlas_dpe_random_pxsigma_restart_dpe3_iter70_dpe1.nii.gz --name "factor_atlas_dpe_random_pxsigma_restart_dpe3_iter70_dpe1" --disabled --overlayType volume --alpha 100.0 --cmap random --negativeCmap greyscale --useNegativeCmap --displayRange 0.0 101.0 --clippingRange 0.0 101.0 --modulateRange 0.0 100.0 --gamma 0.0 --cmapResolution 256 --interpolation none --numSteps 150 --blendFactor 0.1 --smoothing 0 --resolution 100 --numInnerSteps 10 --clipMode intersection --volume 0
```

After saving screenshots from FSLeyes, crop whitespace with ImageMagick:

```bash
mogrify -trim +repage ~/factor-analysis/paper-results/moco_atlas.png ~/factor-analysis/paper-results/moco_cdw_significant.png
```
