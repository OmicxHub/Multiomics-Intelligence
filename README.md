# Clinical–Genomics Harmonisation and AI Baseline Module

> A reproducible R pipeline that harmonises clinical metadata with somatic mutation calls from TCGA-BRCA,
> constructs a machine-learning-ready feature table with a correctly defined five-year survival outcome, and
> benchmarks five baseline classifiers plus a Cox proportional-hazards model.

![R](https://img.shields.io/badge/R-tidyverse%20%7C%20maftools%20%7C%20glmnet-276DC3?style=flat-square&logo=r)
![Cohort](https://img.shields.io/badge/TCGA--BRCA-964%20matched%20patients-185FA5?style=flat-square)
![ML cohort](https://img.shields.io/badge/Confirmed%20outcomes-305-0F3F75?style=flat-square)
![Best AUC](https://img.shields.io/badge/Elastic%20Net%20AUC-0.815-22c55e?style=flat-square)
![Scripts](https://img.shields.io/badge/Scripts-01%20%E2%86%92%2009b-6A3D9A?style=flat-square)

---

## 1. Overview and objective

Clinical annotation and somatic mutation data for the same patients arrive from different files, with different
identifiers, different granularity (one row per patient vs one row per mutation) and different conventions for
missing values. Before any model can be trained, those layers have to be reconciled into a single patient-level
table in which every row is one person and every column means the same thing for everybody.

This module performs that harmonisation for TCGA-BRCA and then asks a concrete downstream question:

> **How well can clinical variables and somatic mutation features predict five-year survival, and which model
> family handles the resulting class imbalance?**

The pipeline runs end to end in R, from raw cBioPortal downloads to trained models, saved figures and result
tables.

### Relationship to the earlier version of this pipeline

This repository documents the **revised pipeline (scripts 01–09b)**, which supersedes an initial version.
Three defects were identified in review and corrected here:

| # | Issue in the initial version | Correction in this revision |
|---|---|---|
| 1 | The five-year survival target counted right-censored patients (alive, followed < 60 months) as "survived" | Only patients with a **confirmed** five-year outcome are labelled; censored patients are excluded (964 → 305 patients) |
| 2 | Three models trained (Logistic Regression, Random Forest, XGBoost) | **Five classifiers** plus a Cox proportional-hazards model |
| 3 | Feature scaling was fitted on all patients **before** the train/test split — data leakage | The split happens first; the scaler is fitted on training data only |

A fourth change is editorial: statements about mutation co-occurrence were rewritten to describe what the data
shows rather than asserting mechanism, since molecular subtypes were not assigned in this module.

---

## 2. Scientific motivation

Breast cancer outcome is driven by factors recorded at different levels of resolution. Stage and age are
recorded per patient; somatic mutations are recorded per variant and must be aggregated into patient-level
quantities such as tumour mutation burden or per-gene flags. Integrating the two is a prerequisite for any
multi-omics prognostic work, and the way the integration is done — which patients are matched, how the outcome
is defined, when scaling is fitted — determines whether downstream model performance is real or an artefact.

The module therefore treats harmonisation and outcome definition as first-class analytical steps, documented and
validated, rather than as plumbing.

---

## 3. Datasets and sources

Three TCGA cohorts were reviewed in script 01 before one was selected.

| Metric | **TCGA-BRCA** | TCGA-LUAD | TCGA-LAML |
|---|---|---|---|
| Role | **Primary dataset** | Reviewed, set aside | Reviewed, set aside |
| Clinical patients | 1,098 | 585 | 200 |
| Patients with mutation calls | 964 | 557 | 199 |
| Total mutations (MAF) | 89,012 | 190,868 | 10,841 |
| Age available | Yes (16 missing) | Yes (94 missing) | **No — entirely missing** |
| AJCC stage available | Yes (12 missing) | Yes (71 missing) | **Not applicable** |
| Overall survival | 945 alive / 152 deceased | 10.8% missing OS_STATUS | 67 alive / 133 deceased |
| Top gene by raw count | PIK3CA | TTN (passenger gene) | DNMT3A |

**TCGA-BRCA was selected** as the largest cohort with the most complete clinical annotation (age, stage, overall
and disease-free survival) and standard AJCC staging. TCGA-LUAD was set aside for this task because of its
missing outcome data and because its top raw-count gene is a known passenger. TCGA-LAML remains a candidate for
cross-cancer testing but lacks age and staging entirely.

An independent cohort, **CPTAC-BRCA 2020** (122 patients), was evaluated as an external test set — see
Section 8.

```
TCGA-BRCA  : https://www.cbioportal.org/study/summary?id=brca_tcga_gdc
TCGA-LUAD  : https://www.cbioportal.org/study/summary?id=luad_tcga_gdc
TCGA-LAML  : https://www.cbioportal.org/study/summary?id=laml_tcga_pan_can_atlas_2018
CPTAC-BRCA : https://www.cbioportal.org/study/summary?id=brca_cptac_2020
```

Each study supplies three files: `data_clinical_patient.txt`, `data_clinical_sample.txt`, `data_mutations.txt`.
All are open-access, de-identified releases; no restricted or patient-identifiable data is used anywhere in this
pipeline.

---

## 4. Harmonisation workflow

### 4.1 Clinical cleaning (script 02)

The clinical patient file (1,098 × 31) and sample file (1,102 × 10) were standardised into **21 analysis-ready
columns**, including `os_event`, `os_months`, `survival_days`, `stage_clean`, `stage_numeric`, pathological
T/N/M stage, year of diagnosis and primary diagnosis.

- `os_event` was derived as 0 = alive, 1 = deceased from `OS_STATUS`; `survival_days` takes `DAYS_TO_DEATH` where
  the patient is deceased and `DAYS_LAST_FOLLOWUP` otherwise.
- **Stage parsing order matters.** Stage strings are matched in the order IV → III → II → I. Checking "Stage II"
  first would misclassify Stage III patients, because `"Stage IIIA"` contains the substring `"Stage II"` — a known
  TCGA cleaning pitfall that silently corrupts the stage variable if missed.

### 4.2 Mutation processing (script 03)

The raw MAF (89,012 records × 140 columns) was filtered to protein-altering variant classes.

| Category | Variant classes | Decision |
|---|---|---|
| Coding — retained | Missense, Nonsense, Frame_Shift_Del, Frame_Shift_Ins, Splice_Site, In_Frame_Del, In_Frame_Ins, Nonstop, Translation_Start_Site | **Kept** — alter the protein product |
| Non-coding — removed | Silent, RNA, Intron, 3'/5'Flank, 3'/5'UTR, Splice_Region, IGR | Removed — synonymous or outside coding regions |

**Result:** 66,687 coding mutations retained (74.9%), 22,325 removed. Across 964 patients this gives a median of
33 coding mutations per patient (mean 69.2), a median tumour mutation burden of 0.87 mut/Mb, 15,344 mutated genes
and 17 hypermutated patients (TMB > 10 mut/Mb).

### 4.3 Patient identifier harmonisation (script 04)

Patient identifiers were extracted from the first 12 characters of the MAF tumour-sample barcode
(`TCGA-AN-A046-01A-11D-A12B-09` → `TCGA-AN-A046`) and joined to the cleaned clinical table.

| Category | N | % of clinical | Note |
|---|---|---|---|
| Total clinical patients | 1,098 | 100% | All patients with clinical metadata |
| Patients with ≥ 1 coding mutation | 964 | 87.8% | After the coding filter |
| **Matched — in both datasets** | **964** | **87.8%** | Enter the merged dataset |
| Clinical only (not sequenced) | 134 | 12.2% | Insufficient tumour material or no exome sequencing |
| Genomic only (no clinical record) | 0 | 0% | Expected: TCGA requires clinical enrolment |

The merged clinical–genomic dataset contains **964 patients × 34 columns**. A match rate of 87.8% is typical for
a TCGA cohort; the 134 unmatched patients are retained in the cleaned clinical file and excluded only from
modelling.

---

## 5. Outcome definition and the cohort correction

No survival-status column exists in the raw data, so `five_year_survival` was engineered from `os_event` and
`os_months` (60 months = 5 years).

| Case | Condition | `five_year_survival` | Interpretation |
|---|---|---|---|
| 1 | `os_event == 0` and `os_months ≥ 60` | 1 | Confirmed survivor — alive beyond five years |
| 2 | `os_event == 1` and `os_months ≥ 60` | 1 | Confirmed survivor — died after five years |
| 3 | `os_event == 1` and `os_months < 60` | 0 | Confirmed non-survivor |
| 4 | `os_event == 0` and `os_months < 60` | `NA` | **Censored — excluded.** Outcome not yet resolved |

Case 4 is the correction. A patient alive after 20 months of follow-up has not been shown to survive five years,
so labelling them a survivor inflates the majority class and trains the model on an assumption rather than an
observation.

| Outcome group | Patients | % of 964 | In ML cohort |
|---|---|---|---|
| Confirmed survived (alive ≥ 60 months) | 173 | 17.9% | Yes (label 1) |
| Confirmed survived (died after 5 years) | 45 | 4.7% | Yes (label 1) |
| Confirmed died within 5 years | 87 | 9.0% | Yes (label 0) |
| **Censored — excluded** | **659** | **68.4%** | No |

**Modelling cohort: 305 patients — 218 survived (71.5%) and 87 did not (28.5%), a class ratio of about 2.5 : 1.**
The correction reduces cohort size by two thirds but changes the class balance from 91%/9% to 71.5%/28.5%, which
makes threshold-dependent metrics such as specificity and MCC informative instead of degenerate.

| | |
|---|---|
| ![Censoring overview](outputs/TCGA-BRCA/figures/fig07_censoring_overview.png) | ![Class balance](outputs/TCGA-BRCA/figures/fig07_class_balance.png) |
| Which patients have a confirmed five-year outcome | Corrected target distribution, n = 305 |

Within the confirmed-outcome cohort, five-year survival falls steeply with stage: **Stage I 86.2%** (n = 58),
Stage II 77.1% (n = 153), Stage III 55.9% (n = 68), **Stage IV 21.4%** (n = 14).

![Survival by stage](outputs/TCGA-BRCA/figures/fig07_survival_by_stage.png)

---

## 6. Feature table

| Group | Features | Type | N |
|---|---|---|---|
| Clinical | `age`, `stage_numeric` (I = 1 … IV = 4), `primary_diagnosis` (one-hot encoded) | Continuous / ordinal / binary | ~6 |
| Genomic summary | `tmb_nonsynonymous`, `n_truncating`, `n_nonsense`, `n_frameshift_del`, `n_frameshift_ins`, `n_splice_site` | Continuous | 6 |
| Gene flags | TP53, PIK3CA, TTN, CDH1, GATA3, MUC16, KMT2C, MAP3K1, HMCN1, FLG, RYR2, SYNE1, PTEN, USH2A, NCOR1, ZFHX4, SPTA1, DMD, MUC17, NEB | Binary (1 = mutated) | 20 |
| Target | `five_year_survival` | Binary | 1 |
| Tracking only | `patient_id` — removed before model fitting | Identifier | 1 |

The 20 gene flags are the most frequently mutated genes in the cohort, included **without** pre-filtering for
driver versus passenger status so that the models themselves weigh their relevance. TP53 (34.3%) and PIK3CA
(34.0%) are co-dominant; large genes such as TTN and MUC16 are also present and are expected to carry passenger
mutations.

![Top mutated genes](outputs/TCGA-BRCA/figures/fig06_top_mutated_genes.png)

---

## 7. Modelling and evaluation

### 7.1 Design

- **Split before scaling.** The 305-patient cohort is split 80/20, stratified on the outcome, seed 42. The scaler
  (mean and standard deviation) is fitted on the training set only and then applied to the test set; binary gene
  flags and one-hot columns are left unscaled.
- **Five classifiers:** Logistic Regression (`glm`), Elastic Net (`glmnet`, α = 0.5), Random Forest, XGBoost and a
  radial-kernel SVM (`e1071`).
- **Survival model:** a Cox proportional-hazards model (`survival`) fitted on the full matched cohort, which uses
  censored follow-up and therefore complements the binary classifiers.
- **Held-out test set:** 62 patients, of whom 18 are non-survivors.

### 7.2 Classifier results (test set, n = 62)

| Model | AUC | Accuracy | Recall | Precision | Specificity | F1 | MCC |
|---|---|---|---|---|---|---|---|
| **Elastic Net** | **0.815** | 0.790 | 1.000 | 0.772 | 0.278 | 0.871 | **0.463** |
| XGBoost | 0.792 | 0.710 | 1.000 | 0.710 | 0.000 | 0.830 | 0.000 |
| SVM | 0.780 | 0.790 | 0.977 | 0.782 | **0.333** | 0.869 | 0.446 |
| Logistic Regression | 0.726 | 0.661 | 0.864 | 0.717 | 0.167 | 0.784 | 0.039 |
| Random Forest | 0.725 | 0.710 | 0.932 | 0.732 | 0.167 | 0.820 | 0.151 |

| | |
|---|---|
| ![ROC curves](outputs/TCGA-BRCA/figures/fig08_roc_overlay.png) | ![Model comparison](outputs/TCGA-BRCA/figures/fig08_model_comparison.png) |
| ROC curves, all five models | Metrics side by side |

**Interpretation.** Elastic Net gives the best discrimination and, unlike XGBoost, retains some ability to
identify non-survivors. XGBoost reaches a comparable AUC (0.792) but at the default 0.5 threshold labels every
test patient a survivor, identifying none of the 18 true non-survivors — a class-imbalance failure mode shared
in weaker form by Logistic Regression and Random Forest. SVM achieves the highest specificity (0.333) while
keeping a strong AUC. Elastic Net's advantage is plausibly its combination of regularisation, which suits a
moderate sample size, with retained class separation. **AUC alone does not capture threshold-dependent clinical
usability**, which is why specificity and MCC are reported alongside it.

Elastic Net selected **11 non-zero features**; stage and age are the strongest predictors of reduced five-year
survival, with several gene flags contributing smaller effects.

![Elastic Net coefficients](outputs/TCGA-BRCA/figures/fig08_en_importance.png)

### 7.3 Cox proportional-hazards model

| Term | Hazard ratio | 95% CI | p |
|---|---|---|---|
| `stage_numeric` | **2.383** | 1.895 – 2.997 | < 0.001 |
| `age` | **1.034** | 1.020 – 1.049 | < 0.001 |
| `diag_mucinous_carcinoma` | 4.855 | 1.505 – 15.659 | 0.008 |
| `DMD` mutated | 2.721 | 1.194 – 6.198 | 0.017 |
| `RYR2` mutated | 0.106 | 0.025 – 0.455 | 0.003 |
| `tmb_nonsynonymous` | 1.239 | 0.955 – 1.608 | 0.107 |
| `KMT2C`, `TTN`, mutation-type counts | 0.978 – 1.390 | overlapping 1 | n.s. |

![Cox forest plot](outputs/TCGA-BRCA/figures/fig08_cox_forest_plot.png)

Stage and age dominate here too, agreeing with the classifier feature importances. The individual gene terms —
mucinous histology, DMD, RYR2 — rest on small numbers of patients and wide confidence intervals and should be
read as **exploratory signals, not established prognostic effects**.

---

## 8. External test on CPTAC-BRCA 2020 — prediction only

The trained models were applied to an independent cohort to test transferability.

| Check | Result |
|---|---|
| Patients (clinical / with mutations) | 122 / 122 |
| Coding mutations | 19,548 (median TMB 1.78 mut/Mb vs 0.87 in TCGA) |
| Overlap of the 20 training genes | **20 / 20 present** |
| `AGE` | Available, 86.1% complete |
| `OS_STATUS`, `OS_MONTHS`, AJCC stage, primary diagnosis | **Not present in the CPTAC release** |

Because the outcome columns are absent, **the five-year outcome cannot be computed in CPTAC and no AUC,
sensitivity or specificity can be calculated there.** Script 09b therefore runs in prediction-only mode: the five
trained models plus an unweighted ensemble generate predicted survival probabilities, with missing stage imputed
as Stage II (the training-set median) and missing diagnosis one-hot columns set to 0.

- 72 of 122 patients (59.0%) fall in the "high survival predicted" band (ensemble probability ≥ 0.75); the
  remaining 50 (41.0%) fall in 0.50–0.75. No patient is predicted below the 0.5 threshold.
- Model agreement is high: all five models agree on the predicted class for 96 of 122 patients (78.7%).
- Predicted probability declines gradually with age, consistent with the age effect learned on TCGA-BRCA.

![CPTAC predicted risk](outputs/TCGA-BRCA/figures/fig09b_risk_distribution.png)

> **This is a plausibility check, not external validation.** Without ground-truth outcomes in CPTAC the results
> describe predicted risk distributions only. The stage imputation further flattens real heterogeneity in that
> cohort. Formal external validation requires a cohort with complete survival annotation.

---

## 9. Technologies and tools

| Layer | Tools |
|---|---|
| Language | R |
| Data handling | tidyverse (dplyr, tidyr, readr, stringr, ggplot2), scales |
| Genomics | maftools (MAF summary, oncoplot, Ti/Tv), cBioPortal flat files |
| Survival analysis | survival, survminer |
| Machine learning | glmnet (Elastic Net), randomForest, xgboost, e1071 (SVM), pROC |
| Visualisation | ggplot2, pheatmap, maftools built-in plots |

Reproducibility: the random seed is fixed at 42 in every script that samples, and the train/test assignment is
written to `outputs/TCGA-BRCA/models/train_test_split_ids.csv` so the split can be reconstructed exactly.

---

## 10. Repository structure

```
01R_clinical_genomics/
├── scripts/
│   ├── 01_data_overview.R              review three TCGA cohorts, select the primary dataset
│   ├── 02_clean_clinical_data.R        31 → 21 standardised clinical columns, stage parsing
│   ├── 03_mutation_processing.R        coding-variant filter, TMB, per-patient mutation counts
│   ├── 04_match_patient_ids.R          barcode harmonisation, clinical ∩ genomic merge
│   ├── 05_generate_summary_tables.R    cohort summaries, gene frequency by stage, TMB tables
│   ├── 06_generate_figures.R           mutation-landscape figures (ggplot2 + maftools)
│   ├── 07_build_ml_feature_table.R     target construction, gene flags, one-hot diagnosis
│   ├── 08_train_ml_models.R            five classifiers + Cox model, leakage-safe scaling
│   ├── 09a_eda_cptac_brca.R            CPTAC variable-availability and overlap audit
│   └── 09b_test_model_on_cptac.R       apply trained models to CPTAC (prediction only)
│
├── data/
│   ├── raw/                            cBioPortal study files (see Section 3)
│   └── processed/TCGA-BRCA/            clinical_cleaned · genomics_cleaned ·
│                                       clinical_genomics_merged · feature_table_ml · feature_table_cox
│
└── outputs/TCGA-BRCA/
    ├── figures/                        dataset review, mutation landscape, target correction,
    │                                   per-model ROC / confusion / importance, Cox forest, CPTAC
    ├── tables/                         numbered result tables, 01_* … 09b_*
    └── models/                         ml1_logistic_regression · ml2_elastic_net · ml3_random_forest ·
                                        ml4_xgboost · ml5_svm · cox_model · feature_scaler · split IDs
```

Key result files: `tables/08_model_performance_all.csv` (the metrics table above),
`tables/08_cox_hazard_ratios.csv`, `tables/07_five_year_outcome_summary.csv` (the cohort correction) and
`tables/09a_cptac_outcome_feasibility.csv` (why CPTAC cannot be scored).

---

## 11. Reproducing the pipeline

```r
install.packages(c("tidyverse", "maftools", "survival", "survminer", "scales",
                   "glmnet", "randomForest", "xgboost", "e1071", "pROC", "pheatmap"))

# 1. Download the study folders from the cBioPortal links in Section 3
#    into data/raw/TCGA-BRCA, data/raw/TCGA-LAML, data/raw/TCGA-LUAD
#    and data/raw/TEST_brca_cptac_2020
# 2. Set the working directory to this folder, then run in order:

source("scripts/01_data_overview.R")
source("scripts/02_clean_clinical_data.R")
source("scripts/03_mutation_processing.R")
source("scripts/04_match_patient_ids.R")
source("scripts/05_generate_summary_tables.R")
source("scripts/06_generate_figures.R")
source("scripts/07_build_ml_feature_table.R")
source("scripts/08_train_ml_models.R")
source("scripts/09a_eda_cptac_brca.R")
source("scripts/09b_test_model_on_cptac.R")
```

Each script prints its progress and writes its tables and figures to `outputs/TCGA-BRCA/`. The scripts set an
absolute working directory in their header; adjust that line to your own path before running.

---

## 12. Limitations

- **Cohort size relative to feature count.** 305 confirmed outcomes support 32 features, 20 of them gene flags.
  Individual gene coefficients (SPTA1, DMD, TTN, RYR2) are exploratory and should not be read as validated
  prognostic markers.
- **659 censored patients are unused** by the classification task. Their five-year outcome is genuinely
  unresolved; the Cox model is the part of this module that can use their follow-up.
- **Class imbalance persists** at 71.5% / 28.5%. At the default 0.5 threshold, XGBoost identifies no
  non-survivors, and Logistic Regression and Random Forest identify few.
- **CPTAC cannot validate the model.** It lacks survival, staging and diagnosis data; stage was imputed as
  Stage II for every patient there, which understates real heterogeneity.
- **Molecular subtypes were not assigned**, so statements about TP53/PIK3CA marking distinct subtypes are
  framed as consistency with published patterns rather than as findings of this module.
- **Mutation calls are taken as provided** by the cBioPortal GDC-harmonised release; no re-calling, VAF filtering
  or tumour-purity adjustment was performed.

---

## 13. Future directions

1. **Address class imbalance directly** — class weighting, threshold tuning on training data, or resampling —
   and report specificity-oriented metrics rather than accuracy.
2. **Obtain a cohort with complete survival annotation** for a formal, AUC-based external validation.
3. **Extend to a second cancer type** (TCGA-LUAD, or TCGA-LAML with its missing age and stage acknowledged) to
   test cross-cancer generalisability.
4. **Use time-to-event modelling more widely**, so that censored patients contribute to the primary analysis
   rather than only to a secondary model.
5. **Add further omics layers** — this pipeline is the clinical–genomic foundation that the DNA methylation
   module builds upon.

---

## 14. Contribution

Pipeline design, implementation, analysis and documentation by **Muteeba Azhar** as the Month 1 deliverable of
the Bioinformatics and AI-Ready Multi-Omics Internship at MultiOmics Intelligence: cohort review and dataset
selection, clinical and mutation data cleaning, identifier harmonisation, outcome engineering and its
correction, feature-table construction, model training and evaluation, the CPTAC feasibility audit, and all
figures, tables and written documentation in this repository.

The revision documented here was prepared in response to a code review that identified the target-definition,
model-coverage and data-leakage issues listed in Section 1.

---

## 15. Acknowledgements and data attribution

- Data generated by **The Cancer Genome Atlas (TCGA) Research Network** (https://www.cancer.gov/tcga) and
  accessed through **cBioPortal for Cancer Genomics**. External test data from the **Clinical Proteomic Tumor
  Analysis Consortium (CPTAC)** BRCA 2020 release, also via cBioPortal.
- Analysis conducted during the Bioinformatics and AI-Ready Multi-Omics Internship at **MultiOmics
  Intelligence**.
- Software: the R Core Team and the authors of tidyverse, maftools, survival, survminer, glmnet, randomForest,
  xgboost, e1071 and pROC.

This repository contains analysis code and derived results from open-access, de-identified public datasets. No
restricted or patient-identifiable data is included.

### References

- Cerami E *et al.* *The cBio Cancer Genomics Portal: an open platform for exploring multidimensional cancer
  genomics data.* Cancer Discovery, 2012.
- Mayakonda A *et al.* *Maftools: efficient and comprehensive analysis of somatic variants in cancer.* Genome
  Research, 2018.
- The Cancer Genome Atlas Network. *Comprehensive molecular portraits of human breast tumours.* Nature, 2012.
- Friedman J, Hastie T, Tibshirani R. *Regularization paths for generalized linear models via coordinate
  descent.* Journal of Statistical Software, 2010.

---

*R · TCGA-BRCA · cBioPortal · scripts 01–09b · five classifiers + Cox model · leakage-safe preprocessing*
