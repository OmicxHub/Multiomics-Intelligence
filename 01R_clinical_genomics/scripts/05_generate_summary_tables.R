# ============================================================
# Script 05: Summary Statistics and Tables
# Project : Clinical-Genomics Harmonisation Pipeline
# Author  : Muteeba Azhar | MOI Internship | July 2026
#
# PURPOSE:
#   Generate comprehensive summary statistics for the
#   matched cohort (n=964). All counts reported with
#   explicit denominator. Language is cautious throughout —
#   no formal statistical tests are performed here;
#   associations are described as "consistent with" rather
#   than "confirms" or "proves".
#
# INPUTS:
#   data/processed/TCGA-BRCA/clinical_genomics_merged.csv
#   outputs/TCGA-BRCA/tables/03_top_mutated_genes.csv
#
# OUTPUTS:
#   outputs/TCGA-BRCA/tables/05_cohort_summary.csv
#   outputs/TCGA-BRCA/tables/05_tmb_by_stage.csv
#   outputs/TCGA-BRCA/tables/05_tmb_by_vital_status.csv
#   outputs/TCGA-BRCA/tables/05_gene_frequency_by_stage.csv
#   outputs/TCGA-BRCA/tables/05_outcome_categories.csv
# ============================================================

library(tidyverse)
library(scales)

setwd("E:/R_Studio/Projects/MOI/01R_clinical_genomics")

PATHS <- list(
  proc_brca       = "data/processed/TCGA-BRCA",
  out_brca_tables = "outputs/TCGA-BRCA/tables"
)

cat("====================================================\n")
cat(" SCRIPT 05: Summary Statistics — TCGA-BRCA\n")
cat(" All counts reported with explicit denominator.\n")
cat(" Language is cautious — no formal tests performed.\n")
cat("====================================================\n\n")


# ═══════════════════════════════════════════════════════════════
# SECTION 1 — Load merged dataset
# ═══════════════════════════════════════════════════════════════

cat("SECTION 1: Loading merged dataset\n")
cat("────────────────────────────────────\n")

merged <- read_csv(
  file.path(PATHS$proc_brca, "clinical_genomics_merged.csv"),
  show_col_types=FALSE)

top_genes <- read_csv(
  file.path(PATHS$out_brca_tables, "03_top_mutated_genes.csv"),
  show_col_types=FALSE)

N_TOTAL <- nrow(merged)

cat(sprintf("Merged dataset loaded: %d patients, %d columns\n\n",
            N_TOTAL, ncol(merged)))

cat("--- merged: str() ---\n")
str(merged)

cat("\n--- merged: head(3) ---\n")
print(head(merged, 3))


# ═══════════════════════════════════════════════════════════════
# SECTION 2 — Overall cohort summary
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 2: Overall cohort summary\n")
cat("────────────────────────────────────\n")

cohort_summary <- tibble(
  Metric  = c(
    "Total matched patients (denominator for all below)",
    "Alive at last contact",
    "Deceased",
    "Missing vital status",
    "Mean age at diagnosis (years)",
    "Median age at diagnosis (years)",
    "Age range (years)",
    "Mean OS months",
    "Median OS months",
    "Stage I — of 964 matched",
    "Stage II — of 964 matched",
    "Stage III — of 964 matched",
    "Stage IV — of 964 matched",
    "Stage Other/Unknown — of 964 matched",
    "Mean TMB (mut/Mb)",
    "Median TMB (mut/Mb)",
    "Hypermutated patients (TMB > 10 mut/Mb)",
    "Total coding mutations",
    "Unique genes mutated"
  ),
  Value = c(
    as.character(N_TOTAL),
    sprintf("%d (%.1f%% of %d)",
            sum(merged$vital_status_clean=="Alive",na.rm=TRUE),
            sum(merged$vital_status_clean=="Alive",na.rm=TRUE)/N_TOTAL*100,
            N_TOTAL),
    sprintf("%d (%.1f%% of %d)",
            sum(merged$vital_status_clean=="Deceased",na.rm=TRUE),
            sum(merged$vital_status_clean=="Deceased",na.rm=TRUE)/N_TOTAL*100,
            N_TOTAL),
    as.character(sum(is.na(merged$vital_status_clean))),
    sprintf("%.1f", mean(merged$age, na.rm=TRUE)),
    sprintf("%.1f", median(merged$age, na.rm=TRUE)),
    sprintf("%.0f – %.0f",
            min(merged$age, na.rm=TRUE),
            max(merged$age, na.rm=TRUE)),
    sprintf("%.1f", mean(merged$os_months, na.rm=TRUE)),
    sprintf("%.1f", median(merged$os_months, na.rm=TRUE)),
    sprintf("%d (%.1f%% of %d)",
            sum(merged$stage_clean=="Stage I",na.rm=TRUE),
            sum(merged$stage_clean=="Stage I",na.rm=TRUE)/N_TOTAL*100,
            N_TOTAL),
    sprintf("%d (%.1f%% of %d)",
            sum(merged$stage_clean=="Stage II",na.rm=TRUE),
            sum(merged$stage_clean=="Stage II",na.rm=TRUE)/N_TOTAL*100,
            N_TOTAL),
    sprintf("%d (%.1f%% of %d)",
            sum(merged$stage_clean=="Stage III",na.rm=TRUE),
            sum(merged$stage_clean=="Stage III",na.rm=TRUE)/N_TOTAL*100,
            N_TOTAL),
    sprintf("%d (%.1f%% of %d)",
            sum(merged$stage_clean=="Stage IV",na.rm=TRUE),
            sum(merged$stage_clean=="Stage IV",na.rm=TRUE)/N_TOTAL*100,
            N_TOTAL),
    sprintf("%d (%.1f%% of %d)",
            sum(merged$stage_clean=="Other/Unknown",na.rm=TRUE),
            sum(merged$stage_clean=="Other/Unknown",na.rm=TRUE)/N_TOTAL*100,
            N_TOTAL),
    sprintf("%.2f", mean(merged$tmb_nonsynonymous, na.rm=TRUE)),
    sprintf("%.2f", median(merged$tmb_nonsynonymous, na.rm=TRUE)),
    sprintf("%d (%.1f%% of %d)",
            sum(merged$tmb_nonsynonymous>10, na.rm=TRUE),
            sum(merged$tmb_nonsynonymous>10, na.rm=TRUE)/N_TOTAL*100,
            N_TOTAL),
    scales::comma(sum(merged$n_total_coding, na.rm=TRUE)),
    "15,344 (from Script 03)"
  )
)

cat("--- Full cohort summary (n=964 matched) ---\n")
print(cohort_summary, n=Inf)

write_csv(cohort_summary,
          file.path(PATHS$out_brca_tables, "05_cohort_summary.csv"))
cat(sprintf("\nSaved: %s/05_cohort_summary.csv\n",
            PATHS$out_brca_tables))


# ═══════════════════════════════════════════════════════════════
# SECTION 3 — Five-year outcome categories
# Corrected definition — censored patients excluded
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 3: Five-year outcome categories (corrected)\n")
cat("─────────────────────────────────────────────────────\n")

outcome_cats <- merged |>
  mutate(
    outcome_5yr = case_when(
      os_event==0 & os_months >= 60 ~
        "Confirmed survived (alive >=60 months)",
      os_event==1 & os_months >= 60 ~
        "Confirmed survived (died after 5yr)",
      os_event==1 & os_months <  60 ~
        "Confirmed died within 5yr",
      os_event==0 & os_months <  60 ~
        "CENSORED — 5yr outcome unknown (EXCLUDED from ML)"
    )
  ) |>
  count(outcome_5yr) |>
  mutate(
    pct_of_964   = round(n/N_TOTAL*100, 1),
    used_in_ML   = !str_detect(outcome_5yr, "CENSORED")
  ) |>
  arrange(desc(used_in_ML), desc(n))

cat("--- five_year_survival categories ---\n")
print(outcome_cats)

n_ml      <- sum(outcome_cats$n[outcome_cats$used_in_ML])
n_survived <- sum(outcome_cats$n[
  outcome_cats$used_in_ML &
  str_detect(outcome_cats$outcome_5yr, "survived")])
n_died     <- sum(outcome_cats$n[
  str_detect(outcome_cats$outcome_5yr, "died within")])

cat(sprintf("\n  ML cohort total    : %d patients\n", n_ml))
cat(sprintf("  Confirmed survived : %d (%.1f%%)\n",
            n_survived, n_survived/n_ml*100))
cat(sprintf("  Confirmed died     : %d (%.1f%%)\n",
            n_died, n_died/n_ml*100))
cat(sprintf("  Class balance      : %.1f:1\n",
            n_survived/n_died))
cat(sprintf("  Censored (excluded): %d (%.1f%% of 964)\n",
            N_TOTAL-n_ml, (N_TOTAL-n_ml)/N_TOTAL*100))

write_csv(outcome_cats,
          file.path(PATHS$out_brca_tables,
                    "05_outcome_categories.csv"))
cat(sprintf("Saved: %s/05_outcome_categories.csv\n",
            PATHS$out_brca_tables))


# ═══════════════════════════════════════════════════════════════
# SECTION 4 — TMB by stage
# CAUTIOUS LANGUAGE: "consistent with" not "proves"
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 4: TMB by cancer stage\n")
cat("────────────────────────────────\n")
cat("Note: descriptive statistics only.\n")
cat("No formal hypothesis test performed.\n\n")

tmb_by_stage <- merged |>
  filter(stage_clean != "Other/Unknown") |>
  group_by(stage_clean) |>
  summarise(
    N            = n(),
    pct_of_964   = round(n()/N_TOTAL*100, 1),
    mean_tmb     = round(mean(tmb_nonsynonymous, na.rm=TRUE), 3),
    median_tmb   = round(median(tmb_nonsynonymous, na.rm=TRUE), 3),
    sd_tmb       = round(sd(tmb_nonsynonymous, na.rm=TRUE), 3),
    n_hypermutated = sum(tmb_nonsynonymous>10, na.rm=TRUE),
    pct_hyper    = round(sum(tmb_nonsynonymous>10,
                             na.rm=TRUE)/n()*100, 1),
    .groups      = "drop"
  ) |>
  arrange(match(stage_clean,
    c("Stage I","Stage II","Stage III","Stage IV")))

cat("--- TMB by stage (descriptive) ---\n")
print(tmb_by_stage)
cat("\nInterpretation: No clear linear relationship between stage\n")
cat("and TMB is observed in this cohort, which is consistent with\n")
cat("published breast cancer literature. Formal statistical testing\n")
cat("(e.g. Kruskal-Wallis) not performed here.\n")

write_csv(tmb_by_stage,
          file.path(PATHS$out_brca_tables, "05_tmb_by_stage.csv"))
cat(sprintf("Saved: %s/05_tmb_by_stage.csv\n",
            PATHS$out_brca_tables))


# ═══════════════════════════════════════════════════════════════
# SECTION 5 — TMB by vital status
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 5: TMB by vital status\n")
cat("────────────────────────────────\n")

tmb_by_vital <- merged |>
  filter(!is.na(vital_status_clean)) |>
  group_by(vital_status_clean) |>
  summarise(
    N          = n(),
    pct_of_964 = round(n()/N_TOTAL*100, 1),
    mean_tmb   = round(mean(tmb_nonsynonymous, na.rm=TRUE), 3),
    median_tmb = round(median(tmb_nonsynonymous, na.rm=TRUE), 3),
    sd_tmb     = round(sd(tmb_nonsynonymous, na.rm=TRUE), 3),
    .groups    = "drop"
  )

cat("--- TMB by vital status (descriptive) ---\n")
print(tmb_by_vital)

write_csv(tmb_by_vital,
          file.path(PATHS$out_brca_tables,
                    "05_tmb_by_vital_status.csv"))
cat(sprintf("Saved: %s/05_tmb_by_vital_status.csv\n",
            PATHS$out_brca_tables))


# ═══════════════════════════════════════════════════════════════
# SECTION 6 — Driver gene frequency by stage
# CAUTIOUS LANGUAGE throughout
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 6: Driver gene frequency by stage\n")
cat("───────────────────────────────────────────\n")
cat("Descriptive frequencies only.\n")
cat("Patterns described as 'consistent with' published literature.\n")
cat("Formal co-mutation or exclusivity tests not performed.\n\n")

# Top 10 driver genes for cross-stage comparison
driver_genes <- top_genes |>
  filter(!Hugo_Symbol %in% c("TTN","MUC16","HMCN1",
                              "FLG","RYR2","SYNE1")) |>
  head(10) |>
  pull(Hugo_Symbol)

cat(sprintf("Driver genes analysed: %s\n\n",
            paste(driver_genes, collapse=", ")))

# Compute frequency per stage for each gene
# We use genomics_cleaned to get gene-level data per patient
genomics_cleaned <- read_csv(
  file.path(PATHS$proc_brca, "genomics_cleaned.csv"),
  show_col_types=FALSE)

# Binary flag per patient per gene
gene_flags <- genomics_cleaned |>
  filter(Hugo_Symbol %in% driver_genes) |>
  group_by(patient_id, Hugo_Symbol) |>
  summarise(mutated=1L, .groups="drop")

# Join with merged to get stage
gene_by_stage <- merged |>
  select(patient_id, stage_clean) |>
  filter(stage_clean != "Other/Unknown") |>
  cross_join(tibble(Hugo_Symbol=driver_genes)) |>
  left_join(gene_flags, by=c("patient_id","Hugo_Symbol")) |>
  mutate(mutated=coalesce(mutated, 0L)) |>
  group_by(stage_clean, Hugo_Symbol) |>
  summarise(
    n_stage      = n_distinct(patient_id),
    n_mutated    = sum(mutated, na.rm=TRUE),
    pct_mutated  = round(sum(mutated,na.rm=TRUE)/
                         n_distinct(patient_id)*100, 1),
    .groups      = "drop"
  ) |>
  arrange(Hugo_Symbol, match(stage_clean,
    c("Stage I","Stage II","Stage III","Stage IV")))

cat("--- Gene frequency by stage (% patients mutated) ---\n")
print(gene_by_stage |>
  pivot_wider(id_cols=Hugo_Symbol,
              names_from=stage_clean,
              values_from=pct_mutated),
  n=Inf)

write_csv(gene_by_stage,
          file.path(PATHS$out_brca_tables,
                    "05_gene_frequency_by_stage.csv"))
cat(sprintf("Saved: %s/05_gene_frequency_by_stage.csv\n",
            PATHS$out_brca_tables))

cat("\nObservation (cautious language):\n")
cat("TP53 mutation frequency appears to increase from Stage I\n")
cat("to Stage III, which is consistent with published breast\n")
cat("cancer literature associating TP53 with more aggressive\n")
cat("disease (TCGA Network, Nature 2012). Formal statistical\n")
cat("testing not performed.\n")
cat("PIK3CA frequency appears to decrease across stages, which\n")
cat("is consistent with its enrichment in Luminal A subtype\n")
cat("(generally lower stage at diagnosis).\n")


# ═══════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════

cat("\n====================================================\n")
cat(" SCRIPT 05 COMPLETE\n")
cat("====================================================\n")
cat(sprintf("  Cohort size             : %d matched patients\n",
            N_TOTAL))
cat(sprintf("  ML analysis cohort      : %d patients\n", n_ml))
cat(sprintf("  Confirmed survived      : %d (%.1f%%)\n",
            n_survived, n_survived/n_ml*100))
cat(sprintf("  Confirmed died          : %d (%.1f%%)\n",
            n_died, n_died/n_ml*100))
cat(sprintf("  Censored (excluded)     : %d\n", N_TOTAL-n_ml))
cat(sprintf("  Class balance (ML)      : %.1f:1\n",
            n_survived/n_died))
cat(sprintf("  Mean TMB                : %.3f mut/Mb\n",
            mean(merged$tmb_nonsynonymous, na.rm=TRUE)))
cat(sprintf("  Hypermutated (TMB>10)   : %d\n",
            sum(merged$tmb_nonsynonymous>10, na.rm=TRUE)))
cat("\n  Tables saved:\n")
cat("    05_cohort_summary.csv\n")
cat("    05_outcome_categories.csv\n")
cat("    05_tmb_by_stage.csv\n")
cat("    05_tmb_by_vital_status.csv\n")
cat("    05_gene_frequency_by_stage.csv\n")
cat("\nNext: source('scripts/06_generate_figures.R')\n")
cat("====================================================\n")
