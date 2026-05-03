# TCGA Multi-Modal SHAP Analysis

> **PhD Thesis Research** | Exposing the *Modality Gap* in standard SHAP-based explainability across five TCGA cancer cohorts using GDC Portal data.

---

## Overview

This repository contains five Jupyter notebooks that execute a rigorous multi-phase audit of SHAP (SHapley Additive exPlanations) applied to multi-modal clinical-genomic survival prediction. The central argument is that standard SHAP suffers from a systematic **Modality Gap** — it over-attributes feature importance to high-dimensional genomic modalities while suppressing the genuine predictive signal held in lower-dimensional clinical features.

The same audit framework is applied independently across five TCGA cancer cohorts, providing cross-cancer evidence that this failure is a **systematic methodological flaw**, not dataset-specific noise.

---

## Cancer Cohorts

| Notebook | Cancer | Abbrev. | Key Clinical Driver | Modalities |
|---|---|---|---|---|
| `TCGA_LUAD_SHAP_Analysis.ipynb` | Lung Adenocarcinoma | LUAD | Stage + EGFR/KRAS + Smoking history | Clinical, Mutation, CNV, RNAseq, miRNA, RPPA |
| `TCGA_LUSC_SHAP_Analysis.ipynb` | Lung Squamous Cell Carcinoma | LUSC | Stage + Pack-years smoked | Clinical, Mutation, CNV, RNAseq, miRNA, RPPA, Methylation |
| `TCGA_OV_SHAP_Analysis.ipynb` | High-Grade Serous Ovarian Carcinoma | OV | Platinum response + Residual disease | Clinical, Mutation, CNV, RNAseq, miRNA, RPPA, Methylation |
| `TCGA_STAD_SHAP_Analysis.ipynb` | Gastric Adenocarcinoma | STAD | Pathologic stage + Lymph nodes | Clinical, Mutation, CNV, RNAseq, miRNA, RPPA |
| `TGCA_COAD_SHAP_Analysis.ipynb` | Colorectal Adenocarcinoma | COAD | Pathologic stage + APC mutation | Clinical, Mutation, CNV, RNAseq, miRNA, RPPA |

---

## The Modality Gap — Core Thesis Argument

Standard SHAP assigns feature importance via a coalition game over all feature subsets. In multi-modal genomic datasets this produces three systematic failure modes:

**1. The Dimensionality Illusion** — High-dimensional modalities (RNA: ~500–1000 features, miRNA: ~500, RPPA: ~200) collectively dominate the SHAP attribution budget regardless of their actual predictive contribution. Clinical features with genuine prognostic value (e.g., `pathologic_stage`, `platinum_response`) are crowded out.

**2. The Modality Attribution Contradiction** — Hierarchical Marginal Gain Analysis (Phase 4) repeatedly shows that modalities ranked highest by SHAP add negligible or even negative Δ AUC when stacked onto clinical features. SHAP says genomic features matter; AUC says they don't.

**3. Attribution Instability** — Bootstrap stability analysis (Phase 3) shows that SHAP importance scores for genomic features have a high Coefficient of Variation (CV) across resamples, confirming they capture sample-specific noise rather than robust biological signal.

### Proposed Solution

The domain-aware SHAP reformulation developed in the thesis addresses these failure modes:

$$\varphi_i = \sum_{S \subseteq F \setminus \{i\}} w(S,i) \cdot n! \cdot |S|! \cdot (n-|S|-1)! \cdot [f(S \cup \{i\} \mid D) - f(S \mid D)]$$

where:
- **w(S,i)** down-weights coalitions containing modality-redundant features
- **D** encodes the dependency structure, including cross-modal biological relationships and sample-level coverage asymmetry
- **f(S|D)** conditions marginal contributions on clinical workflow ordering (Clinical → Mutation → Genomics)

---

## Audit Framework — 4-Phase Pipeline (+ Cancer-Specific Analyses)

Every notebook implements the same four-phase audit, enabling direct cross-cancer comparison.

### Phase 1 — Baseline SHAP
Trains an XGBoost classifier on the full multi-modal feature matrix and computes SHAP values using `TreeExplainer`. Produces a ranked importance plot and identifies which modalities SHAP promotes.

### Phase 2 — ROAR Audit (Remove And Retrain)
Systematically removes the top-ranked SHAP features one by one and measures the resulting AUC degradation. Compares the guided (SHAP-ranked) removal curve against a random removal baseline. A flat curve — where guided removal tracks random removal — is direct evidence that SHAP-nominated features do not hold independent predictive signal.

### Phase 3 — Stability Analysis (Bootstrap)
Runs 10 bootstrap iterations (90% sampling) and records XGBoost feature importances at each iteration. Computes mean importance and Coefficient of Variation (CV) per feature. Clinical features consistently show low CV; genomic features show high CV, confirming they reflect noise rather than stable biology.

### Phase 4 — Hierarchical Marginal Gain Analysis
Builds a cumulative modality stack in clinical workflow order:

`Clinical → +Mutation → +CNV → +RNAseq → +miRNA → +RPPA [→ +Methylation for OV/LUSC]`

Measures the Δ AUC at each step using 5-fold stratified cross-validation. The gap between what SHAP ranks and what actually improves AUC is the quantified modality gap.

### Cancer-Specific Novel Analyses

| Cancer | Novel Analysis | Key Finding |
|---|---|---|
| **LUSC** | miRNA Dual-Quantification Bias | 1046 UUID folders = 523 isoform + 523 mature-miRNA files; every patient has both, creating structural duplication bias |
| **LUSC** | NRF2/KEAP1 Cross-Modality Dependency | NQO1 protein p<0.0001 in NRF2-active vs WT, yet SHAP CV=1.19 — biologically confirmed signal receives unstable attribution |
| **LUSC** | Smoking-TMB Independence | r=0.042, p=0.391 — pack-years and TMB are NOT correlated but capture different carcinogenesis signals; SHAP still over-ranks TMB |
| **LUSC** | Modality Coverage Asymmetry | Missing GDC barcode manifest collapses 4/7 modalities to 0 patients; 100% patient loss at full intersection |
| **OV** | BRCA1 Tri-Modal Redundancy | BRCA1 silencing encoded simultaneously via germline mutation, somatic mutation, and promoter methylation — standard SHAP sums these rather than treating them as one signal |
| **OV** | Coverage Asymmetry | miRNA: ~507 patients vs RNAseq: ~213 patients — zero-imputation for missing modalities creates systematic SHAP bias |
| **LUAD** | Mutual Exclusivity Blindness | KRAS and EGFR are mutually exclusive drivers; SHAP assumes feature independence and cannot model this conditional structure |

---

## Data Source & Format

All data is sourced from the **GDC (Genomic Data Commons) Portal** in GDC format. Each cohort directory follows the structure:

```
TGCA_<COHORT>/
├── barcode/
│   └── TCGA_<COHORT>_file_mapping.csv   # UUID → patient barcode crosswalk
├── Clinical/
│   └── **/*.xml                          # One XML per patient
├── Mutation/
│   └── **/*.maf.gz                       # Ensemble masked MAF (GDC)
├── CNV/
│   └── **/*.gene_level_copy_number.v36.tsv
├── RNAseq/
│   └── **/                              # STAR TPM counts
├── miRNA/
│   └── **/                              # RPM quantification
├── RPPA/
│   └── **/*_RPPA_data.tsv
└── Methylation/  (OV and LUSC only)
    └── **/                              # β-values
```

> **Note for LUSC:** No pre-built barcode manifest is shipped with the GDC download. The notebooks recover the UUID → patient crosswalk via a four-step fallback: pre-built CSV → MAF `Tumor_Sample_Barcode` → RPPA filename pattern → Clinical XML `bcr_patient_barcode`.

### Path Configuration

At the top of each notebook, update the `BASE` path to point to your local or Google Drive data directory:

```python
# Google Colab
BASE = "/content/drive/MyDrive/PHD_dataset_shap/TGCA_LUAD"

# Local
BASE = os.path.expanduser("~/TGCA_LUAD")
```

---

## Requirements

```
python >= 3.9
pandas
numpy
scikit-learn
xgboost
shap
matplotlib
seaborn
lxml
```

Install all dependencies:

```bash
pip install pandas numpy scikit-learn xgboost shap matplotlib seaborn lxml
```

### Library Version Audit

Each notebook includes a version checkpoint cell at startup to ensure reproducibility:

```python
import sys, sklearn, xgboost, shap, pandas, numpy, matplotlib, seaborn
print(f"Python      : {sys.version}")
print(f"scikit-learn: {sklearn.__version__}")
print(f"xgboost     : {xgboost.__version__}")
print(f"shap        : {shap.__version__}")
```

---

## Running the Notebooks

The notebooks are designed to run end-to-end in **Google Colab** with data stored on Google Drive. They can also be run locally by updating the `BASE` path.

```bash
# Clone the repository
git clone https://github.com/indrajeetapache/tgca-luad-lusc-ov-stad-analysis.git
cd tgca-luad-lusc-ov-stad-analysis

# Launch Jupyter
jupyter notebook
```

Open the notebook for your target cohort and run all cells sequentially. Each phase prints diagnostic output and generates inline plots.

---

## Cross-Cancer Evidence Summary

The table below summarises the confirmed SHAP failure patterns across all five cohorts. The consistency across independent cancer types is the core PhD thesis argument.

| Phase | LUAD | LUSC | OV | STAD | COAD |
|---|---|---|---|---|---|
| **SHAP top-ranked modality** | RNA/miRNA | TMB (Mutation) | Methylation/RNA | RNA/RPPA | RNA/RPPA |
| **Clinical AUC alone** | Strong | 0.473 | Strong (residual disease + Pt response) | Strong | Strong |
| **Genomic Δ AUC vs Clinical** | Minimal | +0.067 (RPPA only) | Diminishing returns | Minimal | Minimal |
| **Genomic feature CV (stability)** | High | High (rppa_nqo1 CV=1.19) | High (Methylation) | High | High |
| **Clinical feature CV (stability)** | Low | Low | Low | Low | Low |
| **Mutation modality AUC effect** | Negligible | **Negative (−0.009)** | Negligible | Negligible | Negligible |
| **SHAP vs Marginal Gain contradiction** | ✓ | ✓ (mutation ranked 1st, negative gain) | ✓ | ✓ | ✓ |

---

## Repository Structure

```
tgca-luad-lusc-ov-stad-analysis/
├── README.md
└── notebook/
    ├── coad/
    │   └── TGCA_COAD_SHAP_Analysis.ipynb      # Colorectal Adenocarcinoma
    ├── luad/
    │   └── TCGA_LUAD_SHAP_Analysis.ipynb      # Lung Adenocarcinoma
    ├── lusc/
    │   └── TCGA_LUSC_SHAP_Analysis.ipynb      # Lung Squamous Cell Carcinoma
    ├── ov/
    │   └── TCGA_OV_SHAP_Analysis.ipynb        # Ovarian Carcinoma
    └── stad/
        └── TCGA_STAD_SHAP_Analysis.ipynb      # Gastric Adenocarcinoma
```

