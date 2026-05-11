# MCI-to-AD Conversion Prediction via Multimodal Joint Modeling

A thesis project predicting conversion from Mild Cognitive Impairment (MCI) to Alzheimer's Disease (AD) using a multimodal deep learning autoencoder whose latent features are fed into a Bayesian joint model (JMBayes2) for survival analysis.

---

## Overview

This project combines neuroimaging (MRI), tabular clinical data, and longitudinal cognitive assessments (MMSE, ADAS13) to predict time-to-AD-conversion in MCI patients from the ADNI dataset. A supervised Tensor Fusion autoencoder with multitask learning (MTL) extracts a 128-dimensional latent representation per visit, which is then used as a covariate in a joint longitudinal-survival model fitted in R.

The pipeline is evaluated against six ablation baselines:

| Method | Description |
|---|---|
| Clinical Cox | Age, sex, education, ADAS13 only |
| Image-Only CNN | ResNet18 features, no tabular data |
| Tabular-Only NN | Clinical features, no imaging |
| Concatenation Fusion | Simple feature concatenation |
| No Autoencoder | Raw features, no learned latent space |
| AE-Only (No MTL) | Autoencoder without survival supervision |
| **MAE-JM (Thesis)** | **Full model — multimodal AE + joint model** |

---

## Architecture

### Deep Learning (Python)

The core model is a `SupervisedTemporalFusionAutoencoder` with three components:

- **Image encoder** — attention-gated ResNet18 projecting to 256-dim
- **Tabular encoder** — MLP projecting 16 clinical features to 64-dim
- **Time encoder** — linear projection of visit time to 16-dim

These are fused via **Tensor Fusion** (outer product), producing a `(256+1)×(64+1)×(16+1) = 283,985`-dim vector projected down to a 128-dim latent space `z`. A BiLSTM then processes the sequence of per-visit `z` vectors, and the final hidden state is passed to:

- A **survival head** trained with Efron-corrected Cox partial likelihood
- An **MMSE head** trained with MSE on 12-month future MMSE (the MTL signal)

### Joint Model (R)

The exported latent features (`z_final`, `z_slope`) are joined with longitudinal MMSE trajectories and fitted using `JMBayes2`:

```r
lme(MMSE ~ Years_bl + z_feature, random = ~Years_bl | PTID)
jm(coxFit, lmeFit, functional_forms = ~value(MMSE) + slope(MMSE))
```

This jointly estimates how MMSE trajectory and latent imaging features associate with conversion hazard.

---

## Repository Structure

```
.
├── README.md
│
├── Data/
│   ├── baseline_clinical_features.csv       # Method 1: age, sex, ADAS13, MMSE
│   ├── image_only_features.csv              # Method 2: ResNet18 image features
│   ├── tabular_only_features.csv            # Method 3: tabular-only NN features
│   ├── tabular_patient_level.csv
│   ├── latent_concat.csv                    # Method 4: concatenation fusion (per-visit)
│   ├── concat_patient_level.csv             # Method 4: patient-level summary
│   ├── latent_no_ae.csv                     # Method 5: no autoencoder (per-visit)
│   ├── no_ae_patient_level.csv
│   ├── latent_ae_only.csv                   # Method 6: AE without MTL (per-visit)
│   ├── ae_only_patient_level.csv
│   ├── latent_improved_autoencoder_v6.csv   # Method 7: thesis model (per-visit)
│   └── latent_patient_level_v6.csv          # Method 7: patient-level summary
│
├── Dataset_Creation_Preprocessing/
│   ├── Identify_Valid_Patients.ipynb        # Filter ADNI patients with MCI label + imaging
│   ├── AD_Baseline.ipynb                    # Build baseline_clinical_features.csv
│   ├── AD_ImageOnly.ipynb                   # Extract ResNet18 image features
│   ├── AD_Structured.ipynb                  # Build tabular feature CSVs
│   ├── Autoencoder_Model.ipynb              # AE-only ablation (no MTL)
│   ├── Concatenation_Model.ipynb            # Concatenation fusion ablation
│   ├── Tensor_Fusion_Only.ipynb             # No-autoencoder ablation
│   ├── TensorFusion.py                      # Standalone tensor fusion module
│   └── Thesis_Model.ipynb                   # Full MAE-JM model (main thesis model)
│
├── R/
│   ├── 01_thesis_model.Rmd                  # Method 7 (MAE-JM): the thesis model
│   └── 02_comparison_methods.Rmd            # Methods 1–6: baseline and ablation models
│
└── thesis_text.md
```

---

## Data

All patient data is sourced from the **Alzheimer's Disease Neuroimaging Initiative (ADNI)**. Access requires registration at [adni.loni.usc.edu](https://adni.loni.usc.edu).

The `Data/` directory contains exported latent features and derived clinical summaries — not raw imaging or clinical records. Raw MRI scans and ADNI tabular data must be downloaded separately and placed at the paths expected by the preprocessing notebooks.

### Required Data Files

All data files must be placed in the `Data/` directory before running the R or Python steps.

**For `01_thesis_model.Rmd` (Method 7 — MAE-JM)**

| File | Description |
|------|-------------|
| `latent_patient_level_v6.csv` | Patient-level features (latent z-vectors + survival labels). Must include: `PTID`, `split`, `time_to_event`, `event`, `z_final_*`, `z_slope_*` |
| `latent_improved_autoencoder_v6.csv` | Longitudinal visit-level data. Must include: `PTID`, `split`, `Years_bl`, `MMSE` |

**For `02_comparison_methods.Rmd` (Methods 1–6)**

| File | Description | Used By |
|------|-------------|---------|
| `baseline_clinical_features.csv` | Clinical baseline features + longitudinal MMSE | Method 1 (Clinical Cox) |
| `image_only_features.csv` | CNN image features per patient | Method 2 (Image-Only) |
| `tabular_patient_level.csv` | Tabular deep NN patient-level features | Method 3 (Tabular-Only) |
| `concat_patient_level.csv` | Concatenation fusion patient-level features | Method 4 (Concatenation) |
| `latent_concat.csv` | Longitudinal data for concatenation model | Method 4 (Concatenation) |
| `no_ae_patient_level.csv` | No-autoencoder patient-level features | Method 5 (No AE) |
| `latent_no_ae.csv` | Longitudinal data for no-AE model | Method 5 (No AE) |
| `ae_only_patient_level.csv` | AE-only (no MTL) patient-level features | Method 6 (AE-Only) |
| `latent_ae_only.csv` | Longitudinal data for AE-only model | Method 6 (AE-Only) |

> **Note:** All CSVs must contain a `split` column with values `"train"` and `"val"`. If missing, re-export from the Python pipeline.

### Feature Schema

**Per-visit CSVs** (`latent_improved_autoencoder_v6.csv`, etc.) contain:

| Column | Description |
|---|---|
| `PTID` | Patient ID |
| `Years_bl` | Time from baseline (years) |
| `MMSE` | Mini-Mental State Examination score |
| `time_to_event` | Time from MCI diagnosis to AD conversion or censoring |
| `event` | 1 = converted to AD, 0 = censored |
| `split` | `train` or `val` |
| `z_0` … `z_127` | 128-dim latent visit embedding |

**Patient-level CSVs** (`latent_patient_level_v6.csv`, etc.) additionally contain:

| Column | Description |
|---|---|
| `z_final_0` … `z_final_127` | Final BiLSTM hidden state (sequence summary) |
| `z_slope_0` … `z_slope_127` | Linear slope of each latent dim over visit time |
| `risk_score` | Cox linear predictor from survival head |
| `n_visits` | Number of valid visits used |

---

## Requirements

### Python

```
torch >= 2.0
torchvision
pandas
numpy
scikit-learn
lifelines
Pillow
```

Install with:

```bash
pip install torch torchvision pandas numpy scikit-learn lifelines Pillow
```

### R

```r
install.packages(c(
  "JMbayes2", "nlme", "survival", "survRM2", "boot",
  "pec", "survminer", "glmnet", "ggplot2", "gridExtra",
  "dplyr", "tidyr", "e1071", "scales", "knitr", "rmarkdown"
))
```

---

## Usage

### Step 1 — Identify valid patients

Run `Identify_Valid_Patients.ipynb` to produce `VALID_PATIENTS.pkl`, which filters ADNI patients who have at least one MCI diagnosis, a valid imaging scan, and required tabular covariates.

### Step 2 — Build ablation CSVs

Run the following notebooks in order to produce input CSVs for each comparison method:

```
AD_Baseline.ipynb          → baseline_clinical_features.csv
AD_ImageOnly.ipynb         → image_only_features.csv
AD_Structured.ipynb        → tabular CSVs
Autoencoder_Model.ipynb    → latent_ae_only.csv, ae_only_patient_level.csv
Concatenation_Model.ipynb  → latent_concat.csv, concat_patient_level.csv
Tensor_Fusion_Only.ipynb   → latent_no_ae.csv, no_ae_patient_level.csv
```

### Step 3 — Train the thesis model

Run `Thesis_Model.ipynb` (or `TensorFusion.py` directly). This trains the full MAE-JM model and exports:

```
Data/latent_improved_autoencoder_v6.csv   ← per-visit latent features (train + val)
Data/latent_patient_level_v6.csv          ← patient-level features (train + val)
```

Key config options at the top of the script:

```python
CONFIG = {
    'latent_dim':     128,
    'lstm_hidden':    128,
    'epochs':         100,
    'batch_size':     16,
    'alpha_survival': 0.85,   # Cox loss weight
    'alpha_mmse':     0.15,   # MMSE future prediction weight (MTL)
}
```

> **Note on MTL loss scaling:** MMSE MSE loss is numerically much larger than Cox loss (~200 vs ~1.5 at convergence). Consider normalizing MMSE targets to zero mean / unit variance or reducing `alpha_mmse` to ~0.05 to better balance the two signals.

### Step 4 — Run the thesis joint model (Method 7)

Open `R/01_thesis_model.Rmd` in RStudio and click **Knit**, or run:

```r
rmarkdown::render("R/01_thesis_model.Rmd")
```

This produces:
- An HTML report with embedded output
- All figures saved to `thesis_figures/MSTE-JM/`
- A results CSV at `thesis_figures/ALL_METHODS_RESULTS.csv`

> **Runtime:** Method 7 runs MCMC via `JMbayes2` (2 chains × 10,000 iterations). Expect ~20–40 minutes depending on your machine. All random seeds are set to `42`; results may vary slightly across machines due to floating-point differences in MCMC.

### Step 5 — Run comparison methods in R *(optional)*

```r
rmarkdown::render("R/02_comparison_methods.Rmd")
```

---

## Figures (Joint Modeling Chapter)

All figures are saved to `thesis_figures/MSTE-JM/` after running `01_thesis_model.Rmd`.

| Figure File | Thesis Reference | Generated By |
|-------------|-----------------|--------------|
| `MSTE-JM_jm_survival_curves.png` | Joint model predicted survival curves (converter vs. non-converter) | Section 4 — `plot_joint_model_survival_curves()` |
| `MSTE-JM_baseline_survival.png` | Baseline survival function from Cox sub-model | Section 3 — `basehaz()` plot |
| `MSTE-JM_KM_stratified.png` | Kaplan-Meier curves stratified by predicted risk group | Section 5 — `ggsurvplot()` |
| `MSTE-JM_brier_scores.png` | Brier score at 1, 3, 5 years | Section 4 — `pec()` |
| `MSTE-JM_predicted_vs_observed.png` | Predicted vs. observed conversion time (converters only) | Section 7 — `get_median_survival()` |

---

## Results

Best model (MAE-JM) achieved a validation C-index of **0.8247** at epoch 29 (early stopping at epoch 44). Event rate: ~46% converters in training set, ~48% in validation.

---

## Key Design Decisions

**Why Tensor Fusion instead of simple concatenation?** Tensor fusion captures all pairwise and three-way interactions between image, tabular, and time modalities via the outer product, at the cost of a large intermediate representation (283,985-dim) that is immediately projected down. The ablation in Method 4 (concatenation) tests whether this complexity is justified.

**Why a BiLSTM?** Patient visits are longitudinal sequences of variable length. The BiLSTM reads the sequence of per-visit latent vectors and produces a single patient-level summary that the joint model uses. The backward pass is particularly useful for patients with irregular visit spacing.

**Why JMBayes2?** The joint model correctly accounts for the fact that MMSE is measured with error and is subject to informative dropout (sicker patients drop out), which a naive Cox model on last-observed MMSE would ignore. The association parameters `alpha(value)` and `alpha(slope)` quantify how current MMSE level and rate of decline each independently relate to conversion hazard.

---

## Citation

If you use this code, please cite the ADNI dataset:

> Data used in preparation of this article were obtained from the Alzheimer's Disease Neuroimaging Initiative (ADNI) database (adni.loni.usc.edu). The ADNI was launched in 2003 as a public-private partnership, led by Principal Investigator Michael W. Weiner, MD.

---

## License

MIT
