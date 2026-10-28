# ============================================================
# Script 01: Data Overview and Dataset Selection
# Project : Clinical-Genomics Harmonisation Pipeline
# Dataset : TCGA-BRCA (primary), TCGA-LAML, TCGA-LUAD (reviewed)
# Author  : Muteeba Azhar | MOI Internship | July 2026
#
# PURPOSE:
#   Load all three TCGA datasets, inspect their structure,
#   and select the primary dataset for the pipeline.
#   No data modification is performed in this script.
#
# KEY CORRECTION vs previous pipeline:
#   five_year_survival will be defined ONLY for patients with
#   CONFIRMED 5-year outcomes (os_months >= 60 OR os_event==1).
#   Right-censored patients (alive, os_months < 60) will be
#   EXCLUDED — their 5-year outcome is unknown.
#
# OUTPUTS:
#   outputs/TCGA-BRCA/tables/01_dataset_comparison.csv
# ============================================================


# ═══════════════════════════════════════════════════════════════
# SECTION 0 — Install and load all required packages
# ═══════════════════════════════════════════════════════════════

cat("====================================================\n")
cat(" SECTION 0: Package installation and loading\n")
cat("====================================================\n\n")

# ── CRAN packages ─────────────────────────────────────────────
cran_packages <- c(
  "tidyverse",    # dplyr, ggplot2, readr, stringr, tidyr, purrr
  "survival",     # Cox regression, Surv() objects
  "pROC",         # AUC-ROC calculation and ROC curves
  "randomForest", # Random Forest classifier
  "xgboost",      # XGBoost classifier
  "broom",        # tidy model outputs (glm, coxph)
  "scales",       # axis formatting in ggplot2
  "viridis",      # colour palettes
  "ggrepel",      # non-overlapping text labels in ggplot2
  "pheatmap",     # heatmaps
  "yaml"          # read config.yaml
)

cat("Checking and installing CRAN packages...\n")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing : %s\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    cat(sprintf("  OK         : %s\n", pkg))
  }
}

# ── Bioconductor packages ─────────────────────────────────────
cat("\nChecking Bioconductor packages...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_packages <- c("maftools")   # MAF file analysis + oncoplot

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing (Bioconductor): %s\n", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    cat(sprintf("  OK         : %s\n", pkg))
  }
}

# ── Load all libraries ────────────────────────────────────────
cat("\nLoading all libraries...\n")
suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(pROC)
  library(randomForest)
  library(xgboost)
  library(maftools)
  library(broom)
  library(scales)
  library(viridis)
  library(ggrepel)
  library(pheatmap)
  library(yaml)
})

cat("All packages loaded successfully.\n")
cat(sprintf("R version: %s\n\n", R.version.string))


# ═══════════════════════════════════════════════════════════════
# SECTION 1 — Working directory and paths
# ═══════════════════════════════════════════════════════════════

cat("====================================================\n")
cat(" SECTION 1: Set working directory and paths\n")
cat("====================================================\n\n")

setwd("E:/R_Studio/Projects/MOI/01R_clinical_genomics")
cat(sprintf("Working directory set to:\n  %s\n\n", getwd()))

# ── All project paths (relative to working directory) ─────────
PATHS <- list(
  raw_brca         = "data/raw/TCGA-BRCA",
  raw_laml         = "data/raw/TCGA-LAML",
  raw_luad         = "data/raw/TCGA-LUAD",
  proc_brca        = "data/processed/TCGA-BRCA",
  out_brca_tables  = "outputs/TCGA-BRCA/tables",
  out_brca_figures = "outputs/TCGA-BRCA/figures",
  out_brca_models  = "outputs/TCGA-BRCA/models",
  out_brca_reports = "outputs/TCGA-BRCA/reports",
  out_laml_tables  = "outputs/TCGA-LAML/tables",
  out_laml_figures = "outputs/TCGA-LAML/figures",
  docs             = "docs"
)

cat("Verifying output directories exist:\n")
for (d in unlist(PATHS)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    cat(sprintf("  Created : %s\n", d))
  } else {
    cat(sprintf("  Exists  : %s\n", d))
  }
}

# ── NA encodings used across cBioPortal files ─────────────────
NA_STRINGS <- c("", "NA", "N/A", "n/a",
                "[Not Available]", "[Not Applicable]",
                "[Unknown]", "[Discrepancy]")

cat(sprintf("\nNA strings: %s\n", paste(NA_STRINGS, collapse=" | ")))

# ── Save PATHS to global environment for use in all scripts ───
# (Each script also defines PATHS locally — this is belt-and-braces)
assign("PATHS", PATHS, envir = .GlobalEnv)
assign("NA_STRINGS", NA_STRINGS, envir = .GlobalEnv)


# ═══════════════════════════════════════════════════════════════
# SECTION 2 — Helper function: read cBioPortal files
# ═══════════════════════════════════════════════════════════════

read_cbio <- function(path, show_dims = TRUE) {
  # Reads a cBioPortal tab-separated file safely.
  # Skips # metadata/comment rows at top of file.
  # Applies standard NA string recognition.
  # Args:
  #   path      : character — full file path
  #   show_dims : logical — print dimensions after loading
  # Returns: tibble

  if (!file.exists(path)) {
    stop(sprintf("File not found: %s\nCheck that raw data is in data/raw/", path))
  }

  dat <- read_tsv(
    path,
    comment        = "#",
    na             = NA_STRINGS,
    show_col_types = FALSE,
    progress       = FALSE
  )

  if (show_dims) {
    cat(sprintf("  Loaded : %s\n  Size   : %d rows x %d columns\n",
                basename(path), nrow(dat), ncol(dat)))
  }
  return(dat)
}


# ═══════════════════════════════════════════════════════════════
# SECTION 3 — Load TCGA-BRCA
# ═══════════════════════════════════════════════════════════════

cat("\n====================================================\n")
cat(" SECTION 3: Load TCGA-BRCA\n")
cat("====================================================\n\n")

brca_clin <- read_cbio(file.path(PATHS$raw_brca, "data_clinical_patient.txt"))
brca_samp <- read_cbio(file.path(PATHS$raw_brca, "data_clinical_sample.txt"))
brca_maf  <- read_cbio(file.path(PATHS$raw_brca, "data_mutations.txt"))

cat("\n--- BRCA clinical column names ---\n")
print(colnames(brca_clin))

cat("\n--- BRCA clinical: first 3 rows ---\n")
print(head(brca_clin, 3))

cat("\n--- BRCA clinical: structure ---\n")
str(brca_clin)

cat("\n--- BRCA MAF: first 20 column names ---\n")
print(head(colnames(brca_maf), 20))

# ── Check key survival columns ────────────────────────────────
cat("\n--- BRCA OS_STATUS values ---\n")
print(table(brca_clin$OS_STATUS, useNA = "always"))

cat("\n--- BRCA OS_MONTHS summary ---\n")
print(summary(brca_clin$OS_MONTHS))


# ═══════════════════════════════════════════════════════════════
# SECTION 4 — Load TCGA-LAML
# ═══════════════════════════════════════════════════════════════

cat("\n====================================================\n")
cat(" SECTION 4: Load TCGA-LAML\n")
cat("====================================================\n\n")

laml_clin <- read_cbio(file.path(PATHS$raw_laml, "data_clinical_patient.txt"))
laml_samp <- read_cbio(file.path(PATHS$raw_laml, "data_clinical_sample.txt"))
laml_maf  <- read_cbio(file.path(PATHS$raw_laml, "data_mutations.txt"))

cat("\n--- LAML clinical column names ---\n")
print(colnames(laml_clin))

cat("\n--- LAML: AGE available? ---\n")
cat(sprintf("  AGE column present: %s\n",
            "AGE" %in% colnames(laml_clin)))
if ("AGE" %in% colnames(laml_clin)) {
  cat(sprintf("  Non-missing AGE   : %d / %d\n",
              sum(!is.na(laml_clin$AGE)), nrow(laml_clin)))
}

cat("\n--- LAML: AJCC Stage available? ---\n")
cat(sprintf("  PATH_STAGE present: %s\n",
            "PATH_STAGE" %in% colnames(laml_clin)))


# ═══════════════════════════════════════════════════════════════
# SECTION 5 — Load TCGA-LUAD
# ═══════════════════════════════════════════════════════════════

cat("\n====================================================\n")
cat(" SECTION 5: Load TCGA-LUAD\n")
cat("====================================================\n\n")

luad_clin <- read_cbio(file.path(PATHS$raw_luad, "data_clinical_patient.txt"))
luad_samp <- read_cbio(file.path(PATHS$raw_luad, "data_clinical_sample.txt"))
luad_maf  <- read_cbio(file.path(PATHS$raw_luad, "data_mutations.txt"))

cat("\n--- LUAD clinical column names ---\n")
print(colnames(luad_clin))


# ═══════════════════════════════════════════════════════════════
# SECTION 6 — Dataset comparison table
# ═══════════════════════════════════════════════════════════════

cat("\n====================================================\n")
cat(" SECTION 6: Dataset comparison\n")
cat("====================================================\n\n")

n_avail <- function(dat, col) {
  if (!col %in% colnames(dat)) return("NOT IN FILE")
  n <- sum(!is.na(dat[[col]]))
  sprintf("%d (%.1f%%)", n, n/nrow(dat)*100)
}

os_summ <- function(dat) {
  if (!"OS_STATUS" %in% colnames(dat)) return("NOT IN FILE")
  alive    <- sum(str_detect(coalesce(dat$OS_STATUS,""),"^0"), na.rm=TRUE)
  deceased <- sum(str_detect(coalesce(dat$OS_STATUS,""),"^1"), na.rm=TRUE)
  na_n     <- sum(is.na(dat$OS_STATUS))
  sprintf("%d living / %d deceased / %d NA", alive, deceased, na_n)
}

top_gene <- function(maf_df) {
  if (!"Hugo_Symbol" %in% colnames(maf_df)) return("N/A")
  top <- sort(table(maf_df$Hugo_Symbol), decreasing=TRUE)[1]
  sprintf("%s (%d muts)", names(top), top)
}

comparison_tbl <- tibble(
  Metric = c(
    "Clinical patients (rows)",
    "Clinical columns",
    "Sample records",
    "Sample columns",
    "Total mutations (MAF rows)",
    "MAF columns",
    "Unique patients in MAF",
    "AGE available",
    "AJCC Stage available",
    "OS_STATUS",
    "OS_MONTHS available",
    "DFS_STATUS available",
    "Top mutated gene"
  ),
  TCGA_BRCA = c(
    nrow(brca_clin), ncol(brca_clin),
    nrow(brca_samp), ncol(brca_samp),
    nrow(brca_maf),  ncol(brca_maf),
    n_distinct(str_sub(brca_maf$Tumor_Sample_Barcode, 1, 12)),
    n_avail(brca_clin, "AGE"),
    n_avail(brca_clin, "PATH_STAGE"),
    os_summ(brca_clin),
    n_avail(brca_clin, "OS_MONTHS"),
    n_avail(brca_clin, "DFS_STATUS"),
    top_gene(brca_maf)
  ),
  TCGA_LUAD = c(
    nrow(luad_clin), ncol(luad_clin),
    nrow(luad_samp), ncol(luad_samp),
    nrow(luad_maf),  ncol(luad_maf),
    n_distinct(str_sub(luad_maf$Tumor_Sample_Barcode, 1, 12)),
    n_avail(luad_clin, "AGE"),
    n_avail(luad_clin, "PATH_STAGE"),
    os_summ(luad_clin),
    n_avail(luad_clin, "OS_MONTHS"),
    n_avail(luad_clin, "DFS_STATUS"),
    top_gene(luad_maf)
  ),
  TCGA_LAML = c(
    nrow(laml_clin), ncol(laml_clin),
    nrow(laml_samp), ncol(laml_samp),
    nrow(laml_maf),  ncol(laml_maf),
    n_distinct(str_sub(laml_maf$Tumor_Sample_Barcode, 1, 12)),
    n_avail(laml_clin, "AGE"),
    n_avail(laml_clin, "PATH_STAGE"),
    os_summ(laml_clin),
    n_avail(laml_clin, "OS_MONTHS"),
    n_avail(laml_clin, "DFS_STATUS"),
    top_gene(laml_maf)
  )
)

cat("--- Full dataset comparison ---\n")
print(comparison_tbl, n = Inf, width = Inf)

write_csv(comparison_tbl,
          file.path(PATHS$out_brca_tables, "01_dataset_comparison.csv"))
cat(sprintf("\nSaved: %s/01_dataset_comparison.csv\n",
            PATHS$out_brca_tables))


# ═══════════════════════════════════════════════════════════════
# SECTION 7 — Dataset selection rationale
# ═══════════════════════════════════════════════════════════════

cat("\n====================================================\n")
cat(" SECTION 7: Dataset selection rationale\n")
cat("====================================================\n\n")

cat("PRIMARY SELECTED: TCGA-BRCA\n")
cat("  Reason 1: Largest cohort — most statistical power for ML\n")
cat("  Reason 2: Most complete clinical data (Age >98%, Stage >98%)\n")
cat("  Reason 3: Standard AJCC staging applicable as a feature\n")
cat("  Reason 4: Both OS and DFS available\n\n")

cat("SECONDARY: TCGA-LAML\n")
cat("  Purpose: cross-cancer test — does pipeline structure transfer?\n")
cat("  Known limitation: AGE all missing, AJCC not applicable\n\n")

cat("NOT SELECTED: TCGA-LUAD\n")
cat("  Reason: smaller cohort, more missing OS data than BRCA\n\n")

cat("CRITICAL NOTE — five_year_survival definition:\n")
cat("  Right-censored patients (alive but <60 months follow-up) will be\n")
cat("  EXCLUDED from the binary ML outcome. Their 5-year status is unknown.\n")
cat("  This reduces the analysis cohort but ensures scientific validity.\n")
cat("  Expected cohort sizes:\n")
cat("    Confirmed survived (os_months >= 60) : ~218 patients\n")
cat("    Confirmed died within 5yr            :  ~87 patients\n")
cat("    EXCLUDED (censored)                  : ~659 patients\n")
cat("    TOTAL for ML analysis                : ~305 patients\n")
cat("    New class balance                    : ~2.5:1 (vs old 10:1)\n")


# ═══════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════

cat("\n====================================================\n")
cat(" SCRIPT 01 COMPLETE\n")
cat("====================================================\n")
cat(sprintf("  TCGA-BRCA : %d patients | %d coding muts | %d cols\n",
            nrow(brca_clin), nrow(brca_maf), ncol(brca_maf)))
cat(sprintf("  TCGA-LAML : %d patients | %d mutations\n",
            nrow(laml_clin), nrow(laml_maf)))
cat(sprintf("  TCGA-LUAD : %d patients | %d mutations\n",
            nrow(luad_clin), nrow(luad_maf)))
cat(sprintf("\n  Primary dataset : TCGA-BRCA\n"))
cat(sprintf("  Saved           : %s/01_dataset_comparison.csv\n",
            PATHS$out_brca_tables))
cat("\n  Objects in memory for next script:\n")
cat("    brca_clin, brca_samp, brca_maf\n")
cat("    PATHS, NA_STRINGS\n")
cat("\nNext: source('scripts/02_clean_clinical_data.R')\n")
cat("====================================================\n")


# Dataset comparison bar plot — run in console after Script 01
p_compare <- comparison_tbl |>
  select(Metric, TCGA_BRCA, TCGA_LUAD, TCGA_LAML) |>
  filter(Metric %in% c("Clinical patients (rows)",
                       "Total mutations (MAF rows)",
                       "Unique patients in MAF")) |>
  pivot_longer(-Metric, names_to="Dataset", values_to="Value") |>
  mutate(Value = as.numeric(Value)) |>
  filter(!is.na(Value)) |>
  ggplot(aes(x=Dataset, y=Value, fill=Dataset)) +
  geom_col(width=0.6, colour="white") +
  geom_text(aes(label=scales::comma(Value)),
            vjust=-0.4, size=4, fontface="bold") +
  facet_wrap(~Metric, scales="free_y") +
  scale_fill_manual(values=c("TCGA_BRCA"="#185FA5",
                             "TCGA_LAML"="#EF4444",
                             "TCGA_LUAD"="#F59E0B")) +
  scale_y_continuous(labels=scales::comma,
                     expand=expansion(mult=c(0,0.2))) +
  labs(title="Dataset Comparison — TCGA-BRCA, LAML, LUAD",
       x=NULL, y="Count", fill="Dataset") +
  theme_minimal(base_size=12) +
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="bottom")

print(p_compare)
ggsave("outputs/TCGA-BRCA/figures/fig01_dataset_comparison.png",
       p_compare, width=10, height=5, dpi=300)

