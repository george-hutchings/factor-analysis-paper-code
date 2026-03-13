# Scalable Bayesian Factor Analysis for Multi-Modal MS Data

This repository accompanies *Scalable Bayesian Statistical Machine Learning For Multi-Modal Data With Applications to Multiple Sclerosis* (Hutchings et al., 2026). The code implements a scalable Bayesian factor analysis framework combining a semiparametric Gaussian copula, an Indian buffet process (IBP) prior on latent dimensionality, and a continuous spike-and-slab prior for sparse, interpretable loadings. Monte Carlo Expectation-Maximisation (MCEM) with a C++-accelerated E-step drives inference on large, mixed-modality biomedical datasets such as the NO.MS cohort.

## Citation

If you use this code, please cite the manuscript:

> Hutchings, G., Samartsidis, P., Gaetano, L., Fisher, E., Nichols, T. E., Holmes, C., Häring, D. A., & Ganjgahi, H. (2026). *Scalable Bayesian Statistical Machine Learning For Multi-Modal Data With Applications to Multiple Sclerosis*.

For questions, contact **habib.ganjgahi@bdi.ox.ac.uk**.

---

## Paper Results Guide

This table maps every table and figure in the paper to the code that produces it. Pre-computed simulation results are included in `simulations/*/results.txt` so readers can verify paper numbers without re-running 7,500 SLURM jobs.

| Paper Item | Description | Script(s) | Data Required | Output |
|---|---|---|---|---|
| **Table 1** | Continuous simulation | `simulations/continuous/run.R` → `simulations/analyse.R` | Synthetic | `simulations/continuous/results.txt` |
| **Table 2** | Binary simulation | `simulations/binary/run.R` → `simulations/analyse.R` | Synthetic | `simulations/binary/results.txt` |
| **Table 3** | Mixed/missing simulation | `simulations/mixed/run.R` → `simulations/analyse.R` | Synthetic + `missingness_analysis.rds`† | `simulations/mixed/results.txt` |
| **Table 4** | Clinical loading matrix (Λ) | `noms-analysis/clinical/run.R` → `noms-analysis/clinical/analyse.R` | NO.MS† | Run 1 output |
| **Figure 1** | True Λ (12×3, scenarios 1–2) | Defined in `simulations/continuous/run.R` lines 25–28 | None | — |
| **Figure 2** | True Λ (8×K, scenario 3) | Defined in `simulations/mixed/run.R` lines 34–37 | None | — |
| **Figure 3** | MoCo spatial maps (100 MoCos) | `noms-analysis/mri/analyse_interactive.R` | NO.MS† | `factor_atlas_*.nii.gz` → FSLeyes |
| **Figure 4** | PCA first principal component | Standard `prcomp()` on binary lesion masks | NO.MS† | Not in repo |
| **Figure 5** | CDW-significant MoCos | `noms-analysis/mri/clinical_analysis/analyse_cox_survival.R` | NO.MS† | `coxph_sig_factors_*.nii.gz` → FSLeyes |
| **Cox p-values** (§3.3) | MoCo 48, MoCo 76, VOLT2 | `noms-analysis/mri/clinical_analysis/review_results.R` | NO.MS† | Console output |

†Requires NO.MS data access (see Data Availability below).

### Mapping results.txt rows to paper tables

The `results.txt` files also contain "long v0" variant results — these are additional experiments **not reported in the paper**.

- **Table 1** (continuous): rows `our method`, `factanal K known`, `pca K known`, `rockova 2016`.
- **Table 2** (binary): rows `our method`, `nmf K known`, `li 2023 mirt`.
- **Table 3** (mixed): rows `full data`, `missing data`, `complete cases`.

---

## Data Availability

The NO.MS dataset is proprietary and cannot be distributed. Analyses are divided into two tiers:

**Fully reproducible** (no proprietary data):
- Simulation scenario 1 — continuous data (`simulations/continuous/`)
- Simulation scenario 2 — binary data (`simulations/binary/`)

**Requires NO.MS data access:**
- Simulation scenario 3 — mixed/missing data (needs `missingness_analysis.rds` generated from NO.MS clinical data via `simulations/mixed/missingness-predictor/missingness-clinical.R`)
- Clinical factor analysis (`noms-analysis/clinical/`)
- MRI factor analysis (`noms-analysis/mri/`)
- Cox PH survival models (`noms-analysis/mri/clinical_analysis/`)

For data access enquiries, contact **habib.ganjgahi@bdi.ox.ac.uk**.

---

## Prerequisites

- **R** (≥ 4.3) with `Rcpp` and `RcppArmadillo` (C++ E-step compilation)
- **SLURM** workload manager (for cluster job submission)
- Baseline methods: `NMF`, `mirt` (simulation competitors)
- MRI analysis: `RNifti`, `survival`
- Mixed missingness model: `pROC`, `caret`
- Figure generation: [FSLeyes](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLeyes) and ImageMagick

Hard-coded cluster output paths (e.g. `/data/users/uu85g9/...`, `/well/nichols-nvs/...`) must be updated for your environment. Two SLURM clusters are used: BMRC (continuous/binary simulations) and Novartis (mixed simulations, NO.MS analyses).

---

## Repository Layout

| Path | Description | Paper Reference |
|---|---|---|
| `mcem_fa_algorithm.R` | Core MCEM algorithm with `dpe_func()` for deterministic path expansion | Algorithm 1, §2.1 |
| `src/` | C++ E-step (`mc_e_step.cpp`), truncated normal sampler (`rtn1.*`), RNG (`rng.*`), utilities (`utils.*`) | §2.1 |
| `simulations/continuous/` | Scenario 1: Gaussian data, 12 vars, 3 factors | Table 1, Figure 1 |
| `simulations/binary/` | Scenario 2: binary data, 12 vars, 3 factors | Table 2, Figure 1 |
| `simulations/mixed/` | Scenario 3: mixed types + structured missingness, 8 vars, 3 factors | Table 3, Figure 2 |
| `simulations/methods/` | Baseline implementations: factanal, PCA+varimax, NMF, Li et al. MIRT, Ročková & George 2016 | Tables 1–3 |
| `simulations/shared/evaluate.R` | Evaluation metrics: structure recovery, correct K, loading MSE | Tables 1–3 |
| `simulations/analyse.R` | Aggregates SLURM array outputs into summary tables | Tables 1–3 |
| `noms-analysis/clinical/` | 10-seed clinical factor analysis (9 variables, N=8,023) | Table 4, §3.2 |
| `noms-analysis/mri/` | MRI lesion mask analysis → Modes of lesion Co-occurrence (MoCos) | §3.3, Figure 3 |
| `noms-analysis/mri/clinical_analysis/` | Cox PH survival models and EDSS correlations | Figure 5, §3.3 |
| `paper-results/` | Pre-computed NIfTI atlases, Cox CSV, rendered figures (not tracked — see retrieval instructions) | Figures 3, 5 |

The EDSS correlation scripts (`analyse_edss_correlation.R`, `submit_edss_corr.sh`) are exploratory analyses not reported in the paper.

---

## Simulation Settings

All three scenarios use the same core algorithm (`mcem_fa_algorithm.R`) with these shared settings:

| Parameter | Value |
|-----------|-------|
| v0 (DPE schedule) | 0.1, 0.01, 0.001, 0.0001 |
| v1 | 10 |
| IBP alpha | 2 |
| K (truncation) | 10 |
| N | 200 |
| Realisations | 500 |
| Parameter expansion | None (within DPE stages) |
| Varimax | Applied to loading matrix between DPE stages (not within) |
| Convergence | Max absolute loading change < 0.016 |
| E-step | 100 burn-in + 100 retained samples per iteration |

---

## Simulations

### Scenarios 1–2: Continuous and Binary (Tables 1–2)

Fully reproducible — requires only R + SLURM.

```bash
# First-time setup: create output directories on BMRC cluster
mkdir -p /well/nichols-nvs/users/peo100/factor-analysis/simulations/continuous/logs
mkdir -p /well/nichols-nvs/users/peo100/factor-analysis/simulations/binary/logs

# Continuous (Table 1): 5 methods × 500 realisations = 2,500 jobs
cd ~/factor-analysis/simulations/continuous && sbatch submit-run.sh

# Binary (Table 2): 4 methods × 500 realisations = 2,000 jobs
cd ~/factor-analysis/simulations/binary && sbatch submit-run.sh

# After all jobs complete, aggregate results:
cd ~/factor-analysis/simulations/continuous && sbatch submit-analyse.sh
cd ~/factor-analysis/simulations/binary    && sbatch submit-analyse.sh
```

### Scenario 3: Mixed/Missing (Table 3)

Requires `missingness_analysis.rds` generated from NO.MS data.

```bash
# First-time setup on Novartis cluster
mkdir -p /data/users/uu85g9/factor-analysis/simulations/mixed/logs

# Generate missingness model (requires NO.MS clinical data)
cd ~/factor-analysis/simulations/mixed
Rscript missingness-predictor/missingness-clinical.R

# Run: 6 conditions × 500 realisations = 3,000 jobs
sbatch submit-run.sh

# Aggregate:
sbatch submit-analyse.sh
```

For interactive testing without SLURM, use `simulations/*/interactive.R` — set `slurm_id` to select the (method, realisation) pair, then source the file.

---

## Clinical Analysis (Table 4)

Requires NO.MS clinical data. Recovers K=5 latent dimensions of MS disease status from 9 clinical variables (N=8,023).

**Variables:** EDSS, T25FWM, HPT9M, PASAT, SDMT, VOLT2, NBV2, NUMGDT1, RELAPSE.

**Settings:** DPE v0s = {1, 0.8, 0.4, 0.2, 0.1, 0.05}, K=10, v1=10, 101 EM iterations per stage, 100 MC samples, PX rotation + PX sigma, ibp_alpha=2.

```bash
cd ~/factor-analysis/noms-analysis/clinical
sbatch submit.sh       # 10 random seeds

# After jobs complete:
sbatch --dependency=afterok:<JOBID> --partition=general,legacy \
  --wrap="module load george-bundle/0.1-foss-2023a-R-4.3.2; Rscript ~/factor-analysis/noms-analysis/clinical/analyse.R"
```

Results are saved to `/data/users/uu85g9/factor-analysis/noms/clinical/results_{r}.rds`.

Run 1 is selected for the paper (Table 4). All 10 runs find K=5 factors recovering the same five dimensions:

1. **Physical disability** — EDSS, T25FWM, HPT9M
2. **Cognitive function** — PASAT, SDMT
3. **Asymptomatic activity** — VOLT2, NUMGDT1
4. **Brain damage** — VOLT2, NBV2 (opposite signs)
5. **Relapse** — RELAPSE

---

## MRI Analysis (Figures 3, 5; Cox p-values)

Requires NO.MS MRI binary lesion masks (D=51,595 voxels, N=2,093).

**Settings:** K=100, v0 ∈ {0.1, 0.01, 0.001}, v1=10, convergence tolerance 0.01 (first stage), PX sigma.

### Step 1: Primary run

```bash
cd ~/factor-analysis/noms-analysis/mri/dpe_random_pxsigma
sbatch submit.sh
```

### Step 2: Restart from DPE stage 3 (paper model)

The primary run is restarted at DPE stage 3 (v0=0.001) from iteration 70, selected by inspecting convergence plots (`plot_convergence_stats.R`) which showed parameter differences had plateaued.

```bash
cd ~/factor-analysis/noms-analysis/mri/dpe_random_pxsigma_restart_dpe3
sbatch submit.sh
```

### Step 3: Generate MoCo atlas (Figure 3)

```bash
# Set condition and base_path in the script, then run interactively:
Rscript ~/factor-analysis/noms-analysis/mri/analyse_interactive.R
```

This creates `factor_atlas_*.nii.gz`, visualised in FSLeyes (see Figure Generation below).

### Step 4: Cox PH survival models (Figure 5, Cox p-values)

Compresses subjects into MoCo space and fits Cox models predicting 3-month confirmed disability worsening (M3CDW).

```bash
cd ~/factor-analysis/noms-analysis/mri/clinical_analysis
mkdir -p logs compressed_data results

sbatch submit_compress.sh     # compress subjects into MoCo space
sbatch submit_cox.sh          # Cox PH models

# Print formatted results (HRs, CIs, p-values):
Rscript ~/factor-analysis/noms-analysis/mri/clinical_analysis/review_results.R
```

Paper reports: MoCo 48 (p=0.0058), MoCo 76 (p=0.0083); T2 lesion volume not significant (p=0.5554).

---

## Retrieving Paper Results from HPC

The `paper-results/` directory is not tracked in git. To populate it after running analyses:

```bash
# MoCo atlas NIfTI (Figure 3)
scp uu85g9@ms-login-00:/data/users/uu85g9/factor-analysis/noms/mri/dpe_random_pxsigma_restart_dpe3/factor_atlas_dpe_random_pxsigma_restart_dpe3_iter70_dpe1.nii.gz ~/factor-analysis/paper-results/

# Significant Cox factors NIfTI (Figure 5)
scp uu85g9@ms-login-00:/data/users/uu85g9/factor-analysis/noms/mri/clinical_analysis/results/coxph_sig_factors_restart_dpe3_dpe1_iter70.nii.gz ~/factor-analysis/paper-results/

# Cox results table
scp uu85g9@ms-login-00:/data/users/uu85g9/factor-analysis/noms/mri/clinical_analysis/results/cox_survival_all.csv ~/factor-analysis/paper-results/
```

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

Figures 3 and 5 are generated from NIfTI files using [FSLeyes](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLeyes). After retrieving the NIfTI files (see above), update the file paths in the commands below to match your local `paper-results/` directory:

**Figure 3** — MoCo atlas (`moco_atlas.png`):
```bash
fsleyes --scene lightbox --worldLoc 16.863185978095004 -18.000099182128906 30.609235989629497 --displaySpace /home/hutchings/OneDrive/Documents/academic/24-25/factor-analysis/paper-results/coxph_sig_factors_restart_dpe3_dpe1_iter70.nii.gz --zaxis 2 --zrange 0.3199999999068677 0.6800000000931322 --sliceSpacing 0.029945054973519132 --sampleSlices centre --labelSpace none --sliceOverlap 0.0 --hideCursor --cursorWidth 1.0 --bgColour 0.0 0.0 0.0 --fgColour 1.0 1.0 1.0 --cursorColour 0.0 1.0 0.0 --colourBarLocation top --colourBarLabelSide top-left --colourBarSize 100.0 --labelSize 12 --performance 3 --movieSync /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz --name "MNI152_T1_1mm_brain" --overlayType volume --alpha 100.0 --brightness 49.75000000000001 --contrast 49.90029860765409 --cmap greyscale --negativeCmap greyscale --displayRange 0.0 8447.64 --clippingRange 0.0 8447.64 --modulateRange 0.0 8364.0 --gamma 0.0 --cmapResolution 256 --interpolation none --numSteps 150 --blendFactor 0.1 --smoothing 0 --resolution 100 --numInnerSteps 10 --clipMode intersection --volume 0 /home/hutchings/OneDrive/Documents/academic/24-25/factor-analysis/paper-results/coxph_sig_factors_restart_dpe3_dpe1_iter70.nii.gz --name "coxph_sig_factors_restart_dpe3_dpe1_iter70" --disabled --overlayType volume --alpha 100.0 --brightness 49.74999999999999 --contrast 49.90029860765409 --cmap greyscale --negativeCmap greyscale --displayRange 0.0 76.76 --clippingRange 0.0 76.76 --modulateRange 0.0 76.0 --gamma 0.0 --cmapResolution 256 --interpolation none --numSteps 150 --blendFactor 0.1 --smoothing 0 --resolution 100 --numInnerSteps 10 --clipMode intersection --volume 0 /home/hutchings/OneDrive/Documents/academic/24-25/factor-analysis/paper-results/factor_atlas_dpe_random_pxsigma_restart_dpe3_iter70_dpe1.nii.gz --name "factor_atlas_dpe_random_pxsigma_restart_dpe3_iter70_dpe1" --overlayType volume --alpha 100.0 --brightness 49.75000000000001 --contrast 49.90029860765409 --cmap random --negativeCmap greyscale --displayRange 0.0 101.0 --clippingRange 0.0 101.0 --modulateRange 0.0 100.0 --gamma 0.0 --cmapResolution 256 --interpolation none --numSteps 150 --blendFactor 0.1 --smoothing 0 --resolution 100 --numInnerSteps 10 --clipMode intersection --volume 0
```

**Figure 5** — CDW-significant MoCos (`moco_cdw_significant.png`):
```bash
fsleyes --scene lightbox --worldLoc 16.863185978095004 -18.000099182128906 30.609235989629497 --displaySpace /home/hutchings/OneDrive/Documents/academic/24-25/factor-analysis/paper-results/coxph_sig_factors_restart_dpe3_dpe1_iter70.nii.gz --zaxis 2 --zrange 0.44000000004656614 0.6333333334497486 --sliceSpacing 0.018534798512027186 --sampleSlices centre --labelSpace none --sliceOverlap 0.0 --hideCursor --cursorWidth 1.0 --bgColour 0.0 0.0 0.0 --fgColour 1.0 1.0 1.0 --cursorColour 0.0 1.0 0.0 --colourBarLocation top --colourBarLabelSide top-left --colourBarSize 100.0 --labelSize 12 --performance 3 --movieSync /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz --name "MNI152_T1_1mm_brain" --overlayType volume --alpha 100.0 --brightness 49.75000000000001 --contrast 49.90029860765409 --cmap greyscale --negativeCmap greyscale --displayRange 0.0 8447.64 --clippingRange 0.0 8447.64 --modulateRange 0.0 8364.0 --gamma 0.0 --cmapResolution 256 --interpolation none --numSteps 150 --blendFactor 0.1 --smoothing 0 --resolution 100 --numInnerSteps 10 --clipMode intersection --volume 0 /home/hutchings/OneDrive/Documents/academic/24-25/factor-analysis/paper-results/coxph_sig_factors_restart_dpe3_dpe1_iter70.nii.gz --name "coxph_sig_factors_restart_dpe3_dpe1_iter70" --overlayType volume --alpha 100.0 --brightness 49.74999999999999 --contrast 49.90029860765409 --cmap random --negativeCmap greyscale --displayRange 0.0 76.76 --clippingRange 0.0 76.76 --modulateRange 0.0 76.0 --gamma 0.0 --cmapResolution 256 --interpolation none --numSteps 150 --blendFactor 0.1 --smoothing 0 --resolution 100 --numInnerSteps 10 --clipMode intersection --volume 0 /home/hutchings/OneDrive/Documents/academic/24-25/factor-analysis/paper-results/factor_atlas_dpe_random_pxsigma_restart_dpe3_iter70_dpe1.nii.gz --name "factor_atlas_dpe_random_pxsigma_restart_dpe3_iter70_dpe1" --disabled --overlayType volume --alpha 100.0 --brightness 49.75000000000001 --contrast 49.90029860765409 --cmap random --negativeCmap greyscale --displayRange 0.0 101.0 --clippingRange 0.0 101.0 --modulateRange 0.0 100.0 --gamma 0.0 --cmapResolution 256 --interpolation none --numSteps 150 --blendFactor 0.1 --smoothing 0 --resolution 100 --numInnerSteps 10 --clipMode intersection --volume 0
```

After saving screenshots from FSLeyes, crop whitespace with ImageMagick:

```bash
mogrify -trim +repage ~/factor-analysis/paper-results/moco_atlas.png ~/factor-analysis/paper-results/moco_cdw_significant.png
```
