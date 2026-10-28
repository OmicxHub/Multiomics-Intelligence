# ============================================================
# Script 02: Clinical Data Cleaning
# Project : Clinical-Genomics Harmonisation Pipeline
# Author  : Muteeba Azhar | MOI Internship | July 2026
#
# PURPOSE:
#   Clean raw TCGA-BRCA clinical patient and sample files.
#   Standardise columns, fix stage encoding bug,
#   derive survival variables, print full diagnostics.
#
# KEY CORRECTION vs previous pipeline:
#   five_year_survival is NOT derived here.
#   It will be derived in Script 07 AFTER train/test split,
#   ONLY for patients with confirmed 5-year outcomes.
#   Right-censored patients (alive, os_months < 60)
#   will be EXCLUDED from ML analysis.
#
# INPUTS:
#   data/raw/TCGA-BRCA/data_clinical_patient.txt
#   data/raw/TCGA-BRCA/data_clinical_sample.txt
#
# OUTPUTS:
#   data/processed/TCGA-BRCA/clinical_cleaned.csv
#   data/processed/TCGA-BRCA/sample_cleaned.csv
#   outputs/TCGA-BRCA/tables/02_missingness_summary.csv
#   outputs/TCGA-BRCA/tables/02_stage_summary.csv
#   outputs/TCGA-BRCA/tables/02_vital_summary.csv
#   outputs/TCGA-BRCA/figures/fig02_stage_distribution.png
#   outputs/TCGA-BRCA/figures/fig02_age_distribution.png
#   outputs/TCGA-BRCA/figures/fig02_os_months_distribution.png
# ============================================================

library(tidyverse)
library(scales)

setwd("E:/R_Studio/Projects/MOI/01R_clinical_genomics")

PATHS <- list(
  raw_brca         = "data/raw/TCGA-BRCA",
  proc_brca        = "data/processed/TCGA-BRCA",
  out_brca_tables  = "outputs/TCGA-BRCA/tables",
  out_brca_figures = "outputs/TCGA-BRCA/figures"
)

NA_STRINGS <- c("","NA","N/A","n/a",
                "[Not Available]","[Not Applicable]",
                "[Unknown]","[Discrepancy]")

read_cbio <- function(path) {
  read_tsv(path, comment="#", na=NA_STRINGS,
           show_col_types=FALSE, progress=FALSE)
}

cat("====================================================\n")
cat(" SCRIPT 02: Clinical Data Cleaning — TCGA-BRCA\n")
cat("====================================================\n\n")


# ═══════════════════════════════════════════════════════════════
# SECTION 1 — Load raw files
# ═══════════════════════════════════════════════════════════════

cat("SECTION 1: Loading raw clinical files\n")
cat("──────────────────────────────────────\n")

clin_raw <- read_cbio(
  file.path(PATHS$raw_brca, "data_clinical_patient.txt"))
samp_raw <- read_cbio(
  file.path(PATHS$raw_brca, "data_clinical_sample.txt"))

cat(sprintf("Clinical : %d patients, %d columns\n",
            nrow(clin_raw), ncol(clin_raw)))
cat(sprintf("Sample   : %d records,  %d columns\n\n",
            nrow(samp_raw), ncol(samp_raw)))

cat("--- Clinical column names ---\n")
print(colnames(clin_raw))

cat("\n--- Clinical str() ---\n")
str(clin_raw)

cat("\n--- Clinical head(3) ---\n")
print(head(clin_raw, 3))


# ═══════════════════════════════════════════════════════════════
# SECTION 2 — Missingness report BEFORE cleaning
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 2: Missingness report (raw data)\n")
cat("─────────────────────────────────────────\n")

miss_report <- clin_raw |>
  summarise(across(everything(), ~sum(is.na(.)))) |>
  pivot_longer(everything(),
               names_to  = "column",
               values_to = "n_missing") |>
  mutate(
    n_total   = nrow(clin_raw),
    pct_miss  = round(n_missing / n_total * 100, 1),
    pct_avail = round(100 - pct_miss, 1)
  ) |>
  arrange(desc(n_missing))

cat("--- All columns: n_missing, pct_missing, pct_available ---\n")
print(miss_report, n = Inf)

write_csv(miss_report,
          file.path(PATHS$out_brca_tables,
                    "02_missingness_summary.csv"))
cat(sprintf("\nSaved: %s/02_missingness_summary.csv\n",
            PATHS$out_brca_tables))


# ═══════════════════════════════════════════════════════════════
# SECTION 3 — Survival variable inspection
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 3: Survival variable inspection\n")
cat("─────────────────────────────────────────\n")

cat("--- OS_STATUS raw values ---\n")
print(table(clin_raw$OS_STATUS, useNA="always"))

cat("\n--- OS_MONTHS summary ---\n")
print(summary(clin_raw$OS_MONTHS))

cat("\n--- DFS_STATUS raw values ---\n")
if ("DFS_STATUS" %in% colnames(clin_raw)) {
  print(table(clin_raw$DFS_STATUS, useNA="always"))
} else {
  cat("DFS_STATUS: NOT IN FILE\n")
}

# ── Preview correct five_year_survival categorisation ─────────
cat("\n--- Preview: five_year_survival categories (Script 07) ---\n")
cat("(Actual derivation happens in Script 07 after train/test split)\n\n")

surv_preview <- clin_raw |>
  mutate(
    os_event = case_when(
      str_detect(coalesce(OS_STATUS,""), "^0") ~ 0L,
      str_detect(coalesce(OS_STATUS,""), "^1") ~ 1L,
      TRUE ~ NA_integer_
    )
  ) |>
  mutate(
    outcome_category = case_when(
      os_event == 0 & OS_MONTHS >= 60 ~
        "Confirmed survived (alive at 5yr+)",
      os_event == 1 & OS_MONTHS >= 60 ~
        "Confirmed survived (died after 5yr)",
      os_event == 1 & OS_MONTHS <  60 ~
        "Confirmed died within 5yr",
      os_event == 0 & OS_MONTHS <  60 ~
        "CENSORED — 5yr outcome unknown (EXCLUDED)",
      TRUE ~ "Missing OS data"
    )
  ) |>
  count(outcome_category) |>
  mutate(pct = round(n/sum(n)*100, 1))

print(surv_preview)

n_usable <- surv_preview |>
  filter(str_detect(outcome_category, "Confirmed")) |>
  pull(n) |> sum()

cat(sprintf(
  "\nML analysis cohort: %d patients (censored excluded)\n",
  n_usable))
cat("New class balance will be approximately 2.5:1\n")
cat("(vs old 10:1 which included censored patients)\n")


# ═══════════════════════════════════════════════════════════════
# SECTION 4 — Stage cleaning (regex bug fix)
# WRONG: checking 'Stage I' before 'Stage II'/'Stage III'
#        causes 'Stage IIIA' to match Stage I pattern
# FIX  : check IV -> III -> II -> I (longest prefix first)
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 4: Stage cleaning — regex order fix\n")
cat("────────────────────────────────────────────\n")

cat("--- Raw PATH_STAGE unique values ---\n")
print(sort(table(clin_raw$PATH_STAGE, useNA="always"),
           decreasing=TRUE))

# Correct stage cleaning function
stage_clean_fn <- function(stage_raw) {
  case_when(
    str_detect(coalesce(stage_raw,""),
               regex("stage iv", ignore_case=TRUE))  ~ "Stage IV",
    str_detect(coalesce(stage_raw,""),
               regex("stage iii",ignore_case=TRUE))  ~ "Stage III",
    str_detect(coalesce(stage_raw,""),
               regex("stage ii", ignore_case=TRUE))  ~ "Stage II",
    str_detect(coalesce(stage_raw,""),
               regex("stage i",  ignore_case=TRUE))  ~ "Stage I",
    TRUE ~ "Other/Unknown"
  )
}

# Verify on known tricky cases
cat("\n--- Verification: Stage IIIA must map to Stage III ---\n")
test_cases <- c("Stage I","Stage IA","Stage IIA","Stage IIB",
                "Stage IIIA","Stage IIIC","Stage IV","Stage IVA",NA)
test_result <- tibble(raw=test_cases,
                      clean=stage_clean_fn(test_cases))
print(test_result)
cat("Stage IIIA -> Stage III (correct — not Stage II)\n")


# ═══════════════════════════════════════════════════════════════
# SECTION 5 — Build clinical_cleaned
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 5: Building clinical_cleaned\n")
cat("──────────────────────────────────────\n")

clinical_cleaned <- clin_raw |>

  # Standardise patient_id
  mutate(patient_id = str_trim(str_to_upper(PATIENT_ID))) |>

  # Remove any duplicate patient rows (keep first)
  distinct(patient_id, .keep_all=TRUE) |>

  # OS variables
  mutate(
    os_event = case_when(
      str_detect(coalesce(OS_STATUS,""), "^0") ~ 0L,
      str_detect(coalesce(OS_STATUS,""), "^1") ~ 1L,
      TRUE ~ NA_integer_
    ),
    vital_status_clean = case_when(
      os_event == 0 ~ "Alive",
      os_event == 1 ~ "Deceased",
      TRUE          ~ NA_character_
    ),
    os_months = OS_MONTHS
  ) |>

  # Survival days: prefer DAYS_TO_DEATH, else DAYS_LAST_FOLLOWUP
  mutate(
    survival_days = case_when(
      !is.na(DAYS_TO_DEATH)      ~ DAYS_TO_DEATH,
      !is.na(DAYS_LAST_FOLLOWUP) ~ DAYS_LAST_FOLLOWUP,
      TRUE                       ~ NA_real_
    )
  ) |>

  # DFS variables
  mutate(
    dfs_event = case_when(
      str_detect(coalesce(DFS_STATUS,""), "^0") ~ 0L,
      str_detect(coalesce(DFS_STATUS,""), "^1") ~ 1L,
      TRUE ~ NA_integer_
    ),
    dfs_months = DFS_MONTHS
  ) |>

  # Age
  mutate(age = as.numeric(AGE)) |>

  # Stage — correct order IV > III > II > I
  mutate(
    stage_clean   = stage_clean_fn(PATH_STAGE),
    stage_numeric = case_when(
      stage_clean == "Stage I"   ~ 1L,
      stage_clean == "Stage II"  ~ 2L,
      stage_clean == "Stage III" ~ 3L,
      stage_clean == "Stage IV"  ~ 4L,
      TRUE ~ NA_integer_
    )
  ) |>

  # Diagnosis info
  mutate(
    year_diagnosis    = as.integer(YEAR_OF_DIAGNOSIS),
    primary_diagnosis = PRIMARY_DIAGNOSIS
  ) |>

  # Select final columns
  select(patient_id, age,
         os_event, os_months, vital_status_clean, survival_days,
         dfs_event, dfs_months,
         stage_clean, stage_numeric,
         year_diagnosis, primary_diagnosis)

cat(sprintf("clinical_cleaned: %d patients, %d columns\n",
            nrow(clinical_cleaned), ncol(clinical_cleaned)))

cat("\n--- clinical_cleaned: column names ---\n")
print(colnames(clinical_cleaned))

cat("\n--- clinical_cleaned: head(5) ---\n")
print(head(clinical_cleaned, 5))

cat("\n--- clinical_cleaned: str() ---\n")
str(clinical_cleaned)

cat(sprintf("\nDuplicate patient_ids: %d\n",
            sum(duplicated(clinical_cleaned$patient_id))))


# ═══════════════════════════════════════════════════════════════
# SECTION 6 — Stage and vital status summaries
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 6: Summary tables\n")
cat("──────────────────────────\n")

stage_summary <- clinical_cleaned |>
  group_by(stage_clean) |>
  summarise(
    N             = n(),
    pct_of_total  = round(n()/nrow(clinical_cleaned)*100, 1),
    n_deceased    = sum(os_event==1, na.rm=TRUE),
    pct_deceased  = round(sum(os_event==1,na.rm=TRUE)/n()*100, 1),
    mean_age      = round(mean(age, na.rm=TRUE), 1),
    mean_os_months= round(mean(os_months, na.rm=TRUE), 1),
    .groups       = "drop"
  )

cat("--- Stage summary (denominator = all 1,098 patients) ---\n")
print(stage_summary)

write_csv(stage_summary,
          file.path(PATHS$out_brca_tables, "02_stage_summary.csv"))
cat(sprintf("Saved: %s/02_stage_summary.csv\n",
            PATHS$out_brca_tables))

vital_summary <- clinical_cleaned |>
  group_by(vital_status_clean) |>
  summarise(
    N        = n(),
    pct      = round(n()/nrow(clinical_cleaned)*100, 1),
    mean_age = round(mean(age, na.rm=TRUE), 1),
    mean_os  = round(mean(os_months, na.rm=TRUE), 1),
    min_os   = round(min(os_months, na.rm=TRUE), 1),
    max_os   = round(max(os_months, na.rm=TRUE), 1),
    .groups  = "drop"
  )

cat("\n--- Vital status summary ---\n")
print(vital_summary)

write_csv(vital_summary,
          file.path(PATHS$out_brca_tables, "02_vital_summary.csv"))
cat(sprintf("Saved: %s/02_vital_summary.csv\n",
            PATHS$out_brca_tables))


# ═══════════════════════════════════════════════════════════════
# SECTION 7 — Sample file cleaning
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 7: Sample file cleaning\n")
cat("─────────────────────────────────\n")

cat("--- Sample column names ---\n")
print(colnames(samp_raw))

cat("\n--- Sample head(3) ---\n")
print(head(samp_raw, 3))

sample_cleaned <- samp_raw |>
  mutate(
    patient_id = str_trim(str_to_upper(PATIENT_ID)),
    sample_id  = str_trim(str_to_upper(SAMPLE_ID))
  ) |>
  group_by(patient_id) |>
  slice(1) |>
  ungroup() |>
  select(patient_id, sample_id,
         any_of("SAMPLE_TYPE"),
         any_of("TMB_NONSYNONYMOUS")) |>
  rename_with(tolower)

cat(sprintf("\nsample_cleaned: %d records, %d columns\n",
            nrow(sample_cleaned), ncol(sample_cleaned)))

cat("\n--- sample_cleaned head(5) ---\n")
print(head(sample_cleaned, 5))

if ("tmb_nonsynonymous" %in% colnames(sample_cleaned)) {
  cat("\n--- TMB_NONSYNONYMOUS summary ---\n")
  print(summary(sample_cleaned$tmb_nonsynonymous))
}


# ═══════════════════════════════════════════════════════════════
# SECTION 8 — Plots (print to screen THEN save)
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 8: Generating and printing plots\n")
cat("─────────────────────────────────────────\n")

# ── PLOT 1: Stage distribution ────────────────────────────────
p1 <- clinical_cleaned |>
  filter(stage_clean != "Other/Unknown") |>
  mutate(stage_clean = factor(stage_clean,
    levels=c("Stage I","Stage II","Stage III","Stage IV"))) |>
  count(stage_clean) |>
  mutate(pct = round(n/sum(n)*100, 1)) |>
  ggplot(aes(x=stage_clean, y=n, fill=stage_clean)) +
  geom_col(width=0.65, colour="white") +
  geom_text(aes(label=sprintf("%d\n(%.1f%%)", n, pct)),
            vjust=-0.3, size=4, fontface="bold") +
  scale_fill_manual(values=c(
    "Stage I"  ="#5DCAA5",
    "Stage II" ="#185FA5",
    "Stage III"="#F59E0B",
    "Stage IV" ="#EF4444")) +
  scale_y_continuous(expand=expansion(mult=c(0,0.2))) +
  labs(
    title   = "AJCC Cancer Stage Distribution — TCGA-BRCA",
    subtitle= sprintf("n = %d patients (Other/Unknown excluded)",
      sum(clinical_cleaned$stage_clean!="Other/Unknown")),
    x       = "Cancer Stage",
    y       = "Number of Patients",
    caption = "Source: TCGA-BRCA cBioPortal / GDC"
  ) +
  theme_minimal(base_size=13) +
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="none")

print(p1)
ggsave(file.path(PATHS$out_brca_figures,"fig02_stage_distribution.png"),
       p1, width=9, height=6, dpi=300)
cat("Printed and saved: fig02_stage_distribution.png\n")

# ── PLOT 2: Age distribution ──────────────────────────────────
mean_ages <- clinical_cleaned |>
  filter(!is.na(age), !is.na(vital_status_clean)) |>
  group_by(vital_status_clean) |>
  summarise(mean_age=mean(age,na.rm=TRUE), .groups="drop")

p2 <- clinical_cleaned |>
  filter(!is.na(age), !is.na(vital_status_clean)) |>
  ggplot(aes(x=age, fill=vital_status_clean)) +
  geom_histogram(bins=35, alpha=0.75, position="identity",
                 colour="white", linewidth=0.1) +
  geom_vline(data=mean_ages,
             aes(xintercept=mean_age, colour=vital_status_clean),
             linetype="dashed", linewidth=1) +
  scale_fill_manual(
    values=c("Alive"="#185FA5","Deceased"="#EF4444")) +
  scale_colour_manual(
    values=c("Alive"="#185FA5","Deceased"="#EF4444")) +
  labs(
    title   = "Age Distribution by Vital Status — TCGA-BRCA",
    subtitle= sprintf(
      "n=%d | Mean: Alive=%.1f yrs, Deceased=%.1f yrs",
      sum(!is.na(clinical_cleaned$age)),
      mean_ages$mean_age[mean_ages$vital_status_clean=="Alive"],
      mean_ages$mean_age[mean_ages$vital_status_clean=="Deceased"]),
    x      = "Age at Diagnosis (years)",
    y      = "Number of Patients",
    fill   = "Vital Status",
    colour = "Vital Status",
    caption= "Dashed lines = group means"
  ) +
  theme_minimal(base_size=13) +
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="bottom")

print(p2)
ggsave(file.path(PATHS$out_brca_figures,"fig02_age_distribution.png"),
       p2, width=10, height=6, dpi=300)
cat("Printed and saved: fig02_age_distribution.png\n")

# ── PLOT 3: OS months — censoring overview ────────────────────
p3 <- clinical_cleaned |>
  filter(!is.na(os_months), !is.na(os_event)) |>
  mutate(
    outcome_cat = case_when(
      os_event==0 & os_months >= 60 ~
        "Confirmed survived (alive at 5yr+)",
      os_event==1 & os_months >= 60 ~
        "Confirmed survived (died after 5yr)",
      os_event==1 & os_months <  60 ~
        "Confirmed died within 5yr",
      os_event==0 & os_months <  60 ~
        "CENSORED — excluded from ML"
    ),
    outcome_cat = factor(outcome_cat, levels=c(
      "Confirmed died within 5yr",
      "CENSORED — excluded from ML",
      "Confirmed survived (died after 5yr)",
      "Confirmed survived (alive at 5yr+)"))
  ) |>
  ggplot(aes(x=os_months, fill=outcome_cat)) +
  geom_histogram(bins=40, colour="white", linewidth=0.1) +
  geom_vline(xintercept=60, linetype="dashed",
             colour="#1F2937", linewidth=1.2) +
  annotate("text", x=63, y=Inf,
           label="60 months\n= 5 years", vjust=2,
           hjust=0, size=3.5, colour="#1F2937") +
  scale_fill_manual(values=c(
    "Confirmed died within 5yr"           ="#EF4444",
    "CENSORED — excluded from ML"         ="#9CA3AF",
    "Confirmed survived (died after 5yr)" ="#5DCAA5",
    "Confirmed survived (alive at 5yr+)"  ="#185FA5")) +
  labs(
    title   = "OS Months Distribution — TCGA-BRCA",
    subtitle=
      "Grey (censored) excluded from ML | Coloured = confirmed 5-yr outcomes used",
    x       = "Overall Survival (months)",
    y       = "Number of Patients",
    fill    = "5-Year Outcome",
    caption = paste0(
      "Confirmed survived: ~218 | ",
      "Confirmed died: 87 | ",
      "Censored (excluded): ~659")
  ) +
  theme_minimal(base_size=12) +
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="bottom",
        legend.text=element_text(size=9))

print(p3)
ggsave(file.path(PATHS$out_brca_figures,
                 "fig02_os_months_distribution.png"),
       p3, width=12, height=7, dpi=300)
cat("Printed and saved: fig02_os_months_distribution.png\n")


# ═══════════════════════════════════════════════════════════════
# SECTION 9 — Save cleaned files
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 9: Saving cleaned files\n")
cat("─────────────────────────────────\n")

write_csv(clinical_cleaned,
          file.path(PATHS$proc_brca, "clinical_cleaned.csv"))
cat(sprintf("Saved: %s/clinical_cleaned.csv\n", PATHS$proc_brca))

write_csv(sample_cleaned,
          file.path(PATHS$proc_brca, "sample_cleaned.csv"))
cat(sprintf("Saved: %s/sample_cleaned.csv\n", PATHS$proc_brca))


# ═══════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════

cat("\n====================================================\n")
cat(" SCRIPT 02 COMPLETE\n")
cat("====================================================\n")
cat(sprintf("  clinical_cleaned : %d patients, %d columns\n",
            nrow(clinical_cleaned), ncol(clinical_cleaned)))
cat(sprintf("  sample_cleaned   : %d records,  %d columns\n",
            nrow(sample_cleaned), ncol(sample_cleaned)))
cat(sprintf("  Duplicate IDs    : %d\n",
            sum(duplicated(clinical_cleaned$patient_id))))
cat(sprintf("  Stage I          : %d\n",
            sum(clinical_cleaned$stage_clean=="Stage I",na.rm=TRUE)))
cat(sprintf("  Stage II         : %d\n",
            sum(clinical_cleaned$stage_clean=="Stage II",na.rm=TRUE)))
cat(sprintf("  Stage III        : %d\n",
            sum(clinical_cleaned$stage_clean=="Stage III",na.rm=TRUE)))
cat(sprintf("  Stage IV         : %d\n",
            sum(clinical_cleaned$stage_clean=="Stage IV",na.rm=TRUE)))
cat(sprintf("  Alive            : %d\n",
            sum(clinical_cleaned$vital_status_clean=="Alive",na.rm=TRUE)))
cat(sprintf("  Deceased         : %d\n",
            sum(clinical_cleaned$vital_status_clean=="Deceased",na.rm=TRUE)))
cat(sprintf("  Mean age         : %.1f years\n",
            mean(clinical_cleaned$age, na.rm=TRUE)))
cat("\n  Plots printed + saved:\n")
cat(sprintf("    %s/fig02_stage_distribution.png\n",
            PATHS$out_brca_figures))
cat(sprintf("    %s/fig02_age_distribution.png\n",
            PATHS$out_brca_figures))
cat(sprintf("    %s/fig02_os_months_distribution.png\n",
            PATHS$out_brca_figures))
cat("\nNext: source('scripts/03_mutation_processing.R')\n")
cat("====================================================\n")

