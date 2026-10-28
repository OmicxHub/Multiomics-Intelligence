# ============================================================
# Script 04: Patient ID Matching and Dataset Merge
# Project : Clinical-Genomics Harmonisation Pipeline
# Author  : Muteeba Azhar | MOI Internship | July 2026
#
# PURPOSE:
#   Match clinical and genomic datasets by patient_id.
#   Produce the merged dataset used for all downstream
#   analysis. Document all patient counts with correct
#   denominators to resolve Stage II inconsistency.
#
# KEY FIXES vs previous pipeline:
#   1. All patient counts reported with explicit denominator
#      (e.g. "552/964 matched" not "552 Stage II patients")
#   2. Stage II count discrepancy explained:
#      552 = of 964 matched patients
#      574 = of 1,098 total clinical patients
#   3. Censoring categories printed for downstream awareness
#
# INPUTS:
#   data/processed/TCGA-BRCA/clinical_cleaned.csv
#   data/processed/TCGA-BRCA/sample_cleaned.csv
#   data/processed/TCGA-BRCA/genomics_cleaned.csv
#   outputs/TCGA-BRCA/tables/03_mutation_counts_per_patient.csv
#
# OUTPUTS:
#   data/processed/TCGA-BRCA/clinical_genomics_merged.csv
#   outputs/TCGA-BRCA/tables/04_id_matching_summary.csv
#   outputs/TCGA-BRCA/tables/04_stage_counts_both_denominators.csv
#   outputs/TCGA-BRCA/figures/fig04_matching_overview.png
#   outputs/TCGA-BRCA/figures/fig04_stage_by_denominator.png
# ============================================================

library(tidyverse)
library(scales)

setwd("E:/R_Studio/Projects/MOI/01R_clinical_genomics")

PATHS <- list(
  proc_brca        = "data/processed/TCGA-BRCA",
  out_brca_tables  = "outputs/TCGA-BRCA/tables",
  out_brca_figures = "outputs/TCGA-BRCA/figures"
)

cat("====================================================\n")
cat(" SCRIPT 04: Patient ID Matching — TCGA-BRCA\n")
cat("====================================================\n\n")


# ═══════════════════════════════════════════════════════════════
# SECTION 1 — Load processed files
# ═══════════════════════════════════════════════════════════════

cat("SECTION 1: Loading processed files\n")
cat("────────────────────────────────────\n")

clinical_cleaned <- read_csv(
  file.path(PATHS$proc_brca, "clinical_cleaned.csv"),
  show_col_types=FALSE)

sample_cleaned <- read_csv(
  file.path(PATHS$proc_brca, "sample_cleaned.csv"),
  show_col_types=FALSE)

mutation_counts <- read_csv(
  file.path(PATHS$out_brca_tables,
            "03_mutation_counts_per_patient.csv"),
  show_col_types=FALSE)

genomics_cleaned <- read_csv(
  file.path(PATHS$proc_brca, "genomics_cleaned.csv"),
  show_col_types=FALSE)

cat(sprintf("clinical_cleaned  : %d patients\n",
            nrow(clinical_cleaned)))
cat(sprintf("sample_cleaned    : %d records\n",
            nrow(sample_cleaned)))
cat(sprintf("mutation_counts   : %d patients\n",
            nrow(mutation_counts)))
cat(sprintf("genomics_cleaned  : %d mutations\n",
            nrow(genomics_cleaned)))

cat("\n--- clinical_cleaned: str() ---\n")
str(clinical_cleaned)

cat("\n--- mutation_counts: str() ---\n")
str(mutation_counts)


# ═══════════════════════════════════════════════════════════════
# SECTION 2 — ID matching audit
# Check which patients appear in each dataset
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 2: ID matching audit\n")
cat("──────────────────────────────\n")

ids_clinical <- clinical_cleaned$patient_id
ids_genomic  <- mutation_counts$patient_id
ids_sample   <- sample_cleaned$patient_id

# Patients in both datasets (will form merged)
ids_matched   <- intersect(ids_clinical, ids_genomic)

# Clinical only (no genomic data — not sequenced)
ids_clin_only <- setdiff(ids_clinical, ids_genomic)

# Genomic only (no clinical record — should be 0)
ids_gen_only  <- setdiff(ids_genomic, ids_clinical)

cat(sprintf("Clinical patients total          : %d\n",
            length(ids_clinical)))
cat(sprintf("Genomic patients total           : %d\n",
            length(ids_genomic)))
cat(sprintf("Matched (in BOTH datasets)       : %d (%.1f%% of clinical)\n",
            length(ids_matched),
            length(ids_matched)/length(ids_clinical)*100))
cat(sprintf("Clinical only (no genomic data)  : %d\n",
            length(ids_clin_only)))
cat(sprintf("Genomic only (no clinical record): %d\n",
            length(ids_gen_only)))

if (length(ids_gen_only) > 0) {
  cat("\n--- Genomic-only patient IDs (should be 0) ---\n")
  print(ids_gen_only)
} else {
  cat("All genomic patients have clinical records. Good.\n")
}

# Save matching summary
matching_summary <- tibble(
  Category      = c(
    "Total clinical patients",
    "Total genomic patients (with >=1 coding mutation)",
    "Matched — in both datasets",
    "Clinical only (not sequenced)",
    "Genomic only (no clinical record)"
  ),
  N             = c(
    length(ids_clinical),
    length(ids_genomic),
    length(ids_matched),
    length(ids_clin_only),
    length(ids_gen_only)
  ),
  Pct_of_clinical = round(c(
    100,
    length(ids_genomic)/length(ids_clinical)*100,
    length(ids_matched)/length(ids_clinical)*100,
    length(ids_clin_only)/length(ids_clinical)*100,
    length(ids_gen_only)/length(ids_clinical)*100
  ), 1),
  Notes = c(
    "All patients with clinical metadata",
    "12-char patient_id extraction from MAF barcode",
    "Enter the merged analysis dataset",
    "Excluded: insufficient tumour material or no WES",
    "Should be 0 — TCGA requires clinical enrolment"
  )
)

cat("\n--- ID matching summary ---\n")
print(matching_summary)

write_csv(matching_summary,
          file.path(PATHS$out_brca_tables,
                    "04_id_matching_summary.csv"))
cat(sprintf("Saved: %s/04_id_matching_summary.csv\n",
            PATHS$out_brca_tables))


# ═══════════════════════════════════════════════════════════════
# SECTION 3 — Merge datasets
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 3: Merging clinical + sample + mutation counts\n")
cat("────────────────────────────────────────────────────────\n")

# Step 1: inner join clinical + mutation counts
# Only matched patients are retained
merged <- clinical_cleaned |>
  inner_join(mutation_counts, by="patient_id")

cat(sprintf("After clinical + mutation_counts join: %d patients\n",
            nrow(merged)))

# Step 2: left join sample data (TMB from sample file)
# Left join: keep all matched patients,
# add sample data where available
if ("tmb_nonsynonymous" %in% colnames(sample_cleaned)) {
  merged <- merged |>
    left_join(
      sample_cleaned |>
        select(patient_id, tmb_nonsynonymous) |>
        rename(tmb_sample_file = tmb_nonsynonymous),
      by = "patient_id"
    )
  cat("TMB from sample file added as 'tmb_sample_file'\n")
  cat("TMB from mutation_counts retained as 'tmb_nonsynonymous'\n")
}

cat(sprintf("\nFinal merged dataset: %d patients, %d columns\n",
            nrow(merged), ncol(merged)))

cat("\n--- merged: column names ---\n")
print(colnames(merged))

cat("\n--- merged: head(5) ---\n")
print(head(merged, 5))

cat("\n--- merged: str() ---\n")
str(merged)


# ═══════════════════════════════════════════════════════════════
# SECTION 4 — Stage counts with CORRECT denominators
#
# REVIEWER ISSUE: Stage II count inconsistency
# 552 vs 574 are BOTH correct — different denominators
#   552 = Stage II patients in MATCHED set (964 patients)
#   574 = Stage II patients in ALL clinical (1,098 patients)
# Fix: always state denominator clearly
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 4: Stage counts — both denominators explained\n")
cat("────────────────────────────────────────────────────────\n")

# Stage counts from ALL 1,098 clinical patients
stage_all_clin <- clinical_cleaned |>
  count(stage_clean) |>
  rename(n_all_clinical = n) |>
  mutate(
    pct_of_1098 = round(n_all_clinical/nrow(clinical_cleaned)*100, 1)
  )

# Stage counts from 964 MATCHED patients only
stage_matched <- merged |>
  count(stage_clean) |>
  rename(n_matched = n) |>
  mutate(
    pct_of_964 = round(n_matched/nrow(merged)*100, 1)
  )

# Combined table showing both denominators
stage_both <- stage_all_clin |>
  left_join(stage_matched, by="stage_clean") |>
  arrange(match(stage_clean,
    c("Stage I","Stage II","Stage III","Stage IV","Other/Unknown")))

cat("--- Stage II inconsistency RESOLVED ---\n")
cat("Both values are correct — different denominators:\n\n")
print(stage_both)

cat("\n  Stage II / 1,098 (all clinical)  =",
    stage_both$n_all_clinical[stage_both$stage_clean=="Stage II"],
    "(this was the '574' figure)\n")
cat("  Stage II / 964  (matched)         =",
    stage_both$n_matched[stage_both$stage_clean=="Stage II"],
    "(this was the '552' figure)\n")
cat("\n  All downstream analyses use the MATCHED set (n=964).\n")
cat("  Stage II = 552/964 = %.1f%%\n" |>
      sprintf(
        stage_both$pct_of_964[stage_both$stage_clean=="Stage II"]))

write_csv(stage_both,
          file.path(PATHS$out_brca_tables,
                    "04_stage_counts_both_denominators.csv"))
cat(sprintf("\nSaved: %s/04_stage_counts_both_denominators.csv\n",
            PATHS$out_brca_tables))


# ═══════════════════════════════════════════════════════════════
# SECTION 5 — Five-year outcome preview in merged cohort
# Show exact counts for the corrected five_year_survival
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 5: Five-year outcome preview (matched cohort)\n")
cat("────────────────────────────────────────────────────────\n")

outcome_preview <- merged |>
  mutate(
    outcome_5yr = case_when(
      os_event == 0 & os_months >= 60 ~
        "Confirmed survived (alive at 5yr+)",
      os_event == 1 & os_months >= 60 ~
        "Confirmed survived (died after 5yr)",
      os_event == 1 & os_months <  60 ~
        "Confirmed died within 5yr",
      os_event == 0 & os_months <  60 ~
        "CENSORED — 5yr outcome unknown"
    )
  ) |>
  count(outcome_5yr) |>
  mutate(
    pct_of_matched = round(n/nrow(merged)*100, 1),
    in_ML_analysis = str_detect(outcome_5yr, "Confirmed")
  ) |>
  arrange(desc(in_ML_analysis), desc(n))

cat("--- five_year_survival categories in matched cohort ---\n")
print(outcome_preview)

n_ml <- sum(outcome_preview$n[outcome_preview$in_ML_analysis])
n_survived_ml <- merged |>
  filter((os_event==0 & os_months>=60) |
         (os_event==1 & os_months>=60)) |>
  nrow()
n_died_ml <- merged |>
  filter(os_event==1 & os_months<60) |>
  nrow()

cat(sprintf("\n  ML analysis cohort (confirmed outcomes only): %d patients\n",
            n_ml))
cat(sprintf("  Confirmed survived : %d (%.1f%% of ML cohort)\n",
            n_survived_ml, n_survived_ml/n_ml*100))
cat(sprintf("  Confirmed died     : %d (%.1f%% of ML cohort)\n",
            n_died_ml, n_died_ml/n_ml*100))
cat(sprintf("  Class balance      : %.1f:1\n",
            n_survived_ml/n_died_ml))
cat(sprintf("  Censored (excluded): %d (%.1f%% of matched)\n",
            nrow(merged)-n_ml,
            (nrow(merged)-n_ml)/nrow(merged)*100))


# ═══════════════════════════════════════════════════════════════
# SECTION 6 — Save merged dataset
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 6: Saving merged dataset\n")
cat("──────────────────────────────────\n")

write_csv(merged,
          file.path(PATHS$proc_brca,
                    "clinical_genomics_merged.csv"))
cat(sprintf("Saved: %s/clinical_genomics_merged.csv\n",
            PATHS$proc_brca))
cat(sprintf("  Dimensions: %d patients x %d columns\n",
            nrow(merged), ncol(merged)))


# ═══════════════════════════════════════════════════════════════
# SECTION 7 — Plots (print to screen THEN save)
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 7: Generating and printing plots\n")
cat("─────────────────────────────────────────\n")

# ── PLOT 1: Matching overview waterfall ───────────────────────
p1_df <- tibble(
  Category = c("All clinical\n(1,098)",
               "Not sequenced\n(excluded)",
               "Matched cohort\n(964)"),
  N        = c(nrow(clinical_cleaned),
               length(ids_clin_only),
               length(ids_matched)),
  fill_col = c("#185FA5","#EF4444","#0F6E56")
)

p1 <- ggplot(p1_df, aes(x=Category, y=N, fill=fill_col)) +
  geom_col(width=0.6, colour="white") +
  geom_text(aes(label=scales::comma(N)),
            vjust=-0.4, size=5, fontface="bold") +
  scale_fill_identity() +
  scale_y_continuous(
    labels=scales::comma,
    expand=expansion(mult=c(0,0.2))) +
  labs(
    title   = "Patient ID Matching — TCGA-BRCA",
    subtitle= sprintf(
      "%.1f%% match rate | %d clinical-only excluded",
      length(ids_matched)/length(ids_clinical)*100,
      length(ids_clin_only)),
    x       = NULL,
    y       = "Number of Patients",
    caption = "Clinical-only: not sequenced (insufficient tumour material or no WES)"
  ) +
  theme_minimal(base_size=13) +
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="none")

print(p1)
ggsave(file.path(PATHS$out_brca_figures,
                 "fig04_matching_overview.png"),
       p1, width=8, height=6, dpi=300)
cat("Printed and saved: fig04_matching_overview.png\n")

# ── PLOT 2: Stage counts — both denominators ──────────────────
p2_df <- stage_both |>
  filter(stage_clean != "Other/Unknown") |>
  mutate(stage_clean = factor(stage_clean,
    levels=c("Stage I","Stage II","Stage III","Stage IV"))) |>
  pivot_longer(
    cols      = c(n_all_clinical, n_matched),
    names_to  = "denominator",
    values_to = "n"
  ) |>
  mutate(
    denominator = recode(denominator,
      "n_all_clinical" = "All clinical patients (n=1,098)",
      "n_matched"      = "Matched patients (n=964)")
  )

p2 <- ggplot(p2_df,
             aes(x=stage_clean, y=n,
                 fill=denominator)) +
  geom_col(position="dodge", width=0.7, colour="white") +
  geom_text(aes(label=n),
            position=position_dodge(0.7),
            vjust=-0.4, size=3.5, fontface="bold") +
  scale_fill_manual(
    values=c("All clinical patients (n=1,098)"="#9CA3AF",
             "Matched patients (n=964)"        ="#185FA5")) +
  scale_y_continuous(
    labels=scales::comma,
    expand=expansion(mult=c(0,0.2))) +
  labs(
    title   =
      "Stage Counts — All Clinical vs Matched Patients",
    subtitle=
      "Stage II: 574/1,098 (all clinical) vs 552/964 (matched) — both are correct",
    x       = "Cancer Stage",
    y       = "Number of Patients",
    fill    = "Denominator",
    caption =
      "Discrepancy resolved: different denominators. All downstream analyses use matched set (n=964)."
  ) +
  theme_minimal(base_size=12) +
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="bottom")

print(p2)
ggsave(file.path(PATHS$out_brca_figures,
                 "fig04_stage_by_denominator.png"),
       p2, width=10, height=6, dpi=300)
cat("Printed and saved: fig04_stage_by_denominator.png\n")


# ═══════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════

cat("\n====================================================\n")
cat(" SCRIPT 04 COMPLETE\n")
cat("====================================================\n")
cat(sprintf("  Total clinical patients  : %d\n",
            length(ids_clinical)))
cat(sprintf("  Total genomic patients   : %d\n",
            length(ids_genomic)))
cat(sprintf("  Matched (inner join)     : %d (%.1f%%)\n",
            length(ids_matched),
            length(ids_matched)/length(ids_clinical)*100))
cat(sprintf("  Excluded (no genomics)   : %d\n",
            length(ids_clin_only)))
cat(sprintf("  Merged dimensions        : %d x %d\n",
            nrow(merged), ncol(merged)))
cat("\n  Stage II counts (RESOLVED):\n")
cat(sprintf("    Of ALL clinical (1,098): %d (%.1f%%)\n",
            stage_both$n_all_clinical[stage_both$stage_clean=="Stage II"],
            stage_both$pct_of_1098[stage_both$stage_clean=="Stage II"]))
cat(sprintf("    Of MATCHED  (964)      : %d (%.1f%%)\n",
            stage_both$n_matched[stage_both$stage_clean=="Stage II"],
            stage_both$pct_of_964[stage_both$stage_clean=="Stage II"]))
cat("\n  Five-year ML cohort:\n")
cat(sprintf("    Confirmed survived : %d\n", n_survived_ml))
cat(sprintf("    Confirmed died     : %d\n", n_died_ml))
cat(sprintf("    Censored (excl.)   : %d\n",
            nrow(merged)-n_ml))
cat(sprintf("    ML cohort total    : %d\n", n_ml))
cat(sprintf("    Clas
            s balance      : %.1f:1\n",
            n_survived_ml/n_died_ml))
cat("\n  Plots printed + saved:\n")
cat("    fig04_matching_overview.png\n")
cat("    fig04_stage_by_denominator.png\n")
cat("\nNext: source('scripts/05_generate_summary_tables.R')\n")
cat("====================================================\n")

