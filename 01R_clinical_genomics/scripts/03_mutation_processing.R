# ============================================================
# Script 03: Mutation Data Processing
# Project : Clinical-Genomics Harmonisation Pipeline
# Author  : Muteeba Azhar | MOI Internship | July 2026
#
# PURPOSE:
#   Load raw MAF file, extract patient IDs from barcodes,
#   filter to coding mutations only, compute per-patient
#   mutation counts and tumour mutation burden (TMB),
#   identify top mutated genes, investigate the 964 vs 965
#   patient discrepancy.
#
# KEY FIX vs previous pipeline:
#   Document the 964 vs 965 discrepancy explicitly.
#   maftools reports 965 unique barcodes; our 12-char
#   extraction yields 964 unique patient IDs.
#   This is investigated and explained here.
#
# INPUTS:
#   data/raw/TCGA-BRCA/data_mutations.txt
#
# OUTPUTS:
#   data/processed/TCGA-BRCA/genomics_cleaned.csv
#   outputs/TCGA-BRCA/tables/03_variant_classification_summary.csv
#   outputs/TCGA-BRCA/tables/03_top_mutated_genes.csv
#   outputs/TCGA-BRCA/tables/03_mutation_counts_per_patient.csv
#   outputs/TCGA-BRCA/tables/03_barcode_discrepancy_report.csv
#   outputs/TCGA-BRCA/figures/fig03_variant_classification.png
#   outputs/TCGA-BRCA/figures/fig03_mutations_per_patient.png
#   outputs/TCGA-BRCA/figures/fig03_top_mutated_genes.png
#   outputs/TCGA-BRCA/figures/fig03_tmb_distribution.png
# ============================================================

library(tidyverse)
library(scales)
library(maftools)

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
cat(" SCRIPT 03: Mutation Data Processing — TCGA-BRCA\n")
cat("====================================================\n\n")


# ═══════════════════════════════════════════════════════════════
# SECTION 1 — Load raw MAF file
# ═══════════════════════════════════════════════════════════════

cat("SECTION 1: Loading raw MAF file\n")
cat("────────────────────────────────\n")

maf_raw <- read_cbio(
  file.path(PATHS$raw_brca, "data_mutations.txt"))

cat(sprintf("MAF raw: %d mutations, %d columns\n\n",
            nrow(maf_raw), ncol(maf_raw)))

cat("--- MAF column names ---\n")
print(colnames(maf_raw))

cat("\n--- MAF head(3) ---\n")
print(head(maf_raw, 3))

cat("\n--- MAF str() (first 15 columns) ---\n")
str(maf_raw[, 1:min(15, ncol(maf_raw))])

# ── Key columns check ─────────────────────────────────────────
cat("\n--- Required columns present? ---\n")
required_cols <- c("Hugo_Symbol","Variant_Classification",
                   "Tumor_Sample_Barcode","HGVSp_Short")
for (col in required_cols) {
  cat(sprintf("  %-30s : %s\n", col,
              ifelse(col %in% colnames(maf_raw), "YES", "NO")))
}


# ═══════════════════════════════════════════════════════════════
# SECTION 2 — Patient ID extraction + 964 vs 965 investigation
#
# TCGA barcodes are 28 characters long:
#   TCGA-AN-A046-01A-11D-A12B-09
#   ^^^^^^^^^^^^  <- first 12 chars = patient ID
#
# maftools may report 965 unique barcodes because:
#   1. One patient may have two tumour samples (different barcode
#      suffixes but same 12-char patient ID), OR
#   2. A barcode has a formatting anomaly
# We investigate here and document the finding.
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 2: Patient ID extraction + 964/965 investigation\n")
cat("──────────────────────────────────────────────────────────\n")

# maftools unique barcode count
n_unique_barcodes <- n_distinct(maf_raw$Tumor_Sample_Barcode)
cat(sprintf("Unique full barcodes (28-char) : %d\n",
            n_unique_barcodes))

# Extract 12-char patient IDs
maf_raw <- maf_raw |>
  mutate(patient_id = str_sub(Tumor_Sample_Barcode, 1, 12) |>
           str_trim() |> str_to_upper())

n_unique_patients <- n_distinct(maf_raw$patient_id)
cat(sprintf("Unique patient IDs (12-char)   : %d\n",
            n_unique_patients))
cat(sprintf("Discrepancy                    : %d\n",
            n_unique_barcodes - n_unique_patients))

# ── Find the patient(s) causing the discrepancy ───────────────
if (n_unique_barcodes != n_unique_patients) {
  cat("\n--- Investigating discrepancy ---\n")

  # Find patient IDs that appear with more than one barcode
  multi_barcode <- maf_raw |>
    group_by(patient_id) |>
    summarise(
      n_barcodes = n_distinct(Tumor_Sample_Barcode),
      barcodes   = paste(unique(Tumor_Sample_Barcode),
                         collapse=" | "),
      n_mutations= n(),
      .groups    = "drop"
    ) |>
    filter(n_barcodes > 1)

  cat(sprintf("Patients with >1 barcode: %d\n",
              nrow(multi_barcode)))
  print(multi_barcode)

  # Save discrepancy report
  write_csv(multi_barcode,
            file.path(PATHS$out_brca_tables,
                      "03_barcode_discrepancy_report.csv"))
  cat(sprintf("Saved: %s/03_barcode_discrepancy_report.csv\n",
              PATHS$out_brca_tables))

  cat("\nExplanation: maftools counts unique Tumor_Sample_Barcode\n")
  cat("values (full 28-char string). Our 12-char extraction maps\n")
  cat("multiple barcodes to the same patient ID when a patient has\n")
  cat("more than one tumour sample sequenced.\n")
  cat("Decision: use 12-char patient_id — this is the correct\n")
  cat("clinical unit. Both barcodes' mutations are KEPT.\n")
} else {
  cat("No discrepancy — barcodes and patient IDs match exactly.\n")
}


# ═══════════════════════════════════════════════════════════════
# SECTION 3 — Variant classification breakdown
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 3: Variant classification breakdown\n")
cat("────────────────────────────────────────────\n")

# Define coding mutations to KEEP
# These change the protein sequence and may have functional impact
CODING_TYPES <- c(
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "Splice_Site",
  "In_Frame_Del",
  "In_Frame_Ins",
  "Nonstop_Mutation",
  "Translation_Start_Site"
)

# All other types are non-coding and removed:
# Silent (synonymous), RNA, Intron, 3'UTR, 5'UTR,
# 3'Flank, 5'Flank, Splice_Region, IGR

var_class_summary <- maf_raw |>
  count(Variant_Classification, name="n_mutations") |>
  mutate(
    pct      = round(n_mutations / sum(n_mutations) * 100, 2),
    is_coding= Variant_Classification %in% CODING_TYPES,
    decision = if_else(is_coding, "KEPT (coding)", "REMOVED (non-coding)")
  ) |>
  arrange(desc(n_mutations))

cat("--- All variant classifications ---\n")
print(var_class_summary, n=Inf)

cat(sprintf("\nTotal mutations (raw)     : %d\n", nrow(maf_raw)))
cat(sprintf("Coding mutations (kept)   : %d (%.1f%%)\n",
    sum(var_class_summary$n_mutations[var_class_summary$is_coding]),
    sum(var_class_summary$n_mutations[var_class_summary$is_coding])/
      nrow(maf_raw)*100))
cat(sprintf("Non-coding (removed)      : %d (%.1f%%)\n",
    sum(var_class_summary$n_mutations[!var_class_summary$is_coding]),
    sum(var_class_summary$n_mutations[!var_class_summary$is_coding])/
      nrow(maf_raw)*100))

write_csv(var_class_summary,
          file.path(PATHS$out_brca_tables,
                    "03_variant_classification_summary.csv"))
cat(sprintf("Saved: %s/03_variant_classification_summary.csv\n",
            PATHS$out_brca_tables))


# ═══════════════════════════════════════════════════════════════
# SECTION 4 — Filter to coding mutations only
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 4: Filter to coding mutations\n")
cat("───────────────────────────────────────\n")

maf_coding <- maf_raw |>
  filter(Variant_Classification %in% CODING_TYPES)

cat(sprintf("After filter: %d coding mutations\n",
            nrow(maf_coding)))
cat(sprintf("Unique patients with coding mutations: %d\n",
            n_distinct(maf_coding$patient_id)))

cat("\n--- Coding mutations: str() ---\n")
str(maf_coding)


# ═══════════════════════════════════════════════════════════════
# SECTION 5 — Per-patient mutation counts and TMB
#
# TMB = Tumour Mutation Burden
# Standard denominator: 38 Mb (human exome capture target size)
# TMB > 10 mut/Mb = hypermutated (FDA threshold for
# pembrolizumab eligibility regardless of tumour type)
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 5: Per-patient mutation counts and TMB\n")
cat("──────────────────────────────────────────────\n")

EXOME_SIZE_MB <- 38  # standard exome size denominator

mutation_counts <- maf_coding |>
  group_by(patient_id) |>
  summarise(
    # Total coding mutations
    n_total_coding   = n(),

    # By variant type
    n_missense       = sum(Variant_Classification ==
                             "Missense_Mutation",    na.rm=TRUE),
    n_nonsense       = sum(Variant_Classification ==
                             "Nonsense_Mutation",    na.rm=TRUE),
    n_frameshift_del = sum(Variant_Classification ==
                             "Frame_Shift_Del",      na.rm=TRUE),
    n_frameshift_ins = sum(Variant_Classification ==
                             "Frame_Shift_Ins",      na.rm=TRUE),
    n_splice_site    = sum(Variant_Classification ==
                             "Splice_Site",          na.rm=TRUE),
    n_inframe_del    = sum(Variant_Classification ==
                             "In_Frame_Del",         na.rm=TRUE),
    n_inframe_ins    = sum(Variant_Classification ==
                             "In_Frame_Ins",         na.rm=TRUE),

    # Truncating = nonsense + frameshift (protein is cut short)
    n_truncating     = sum(Variant_Classification %in%
                             c("Nonsense_Mutation",
                               "Frame_Shift_Del",
                               "Frame_Shift_Ins"),   na.rm=TRUE),

    # Unique genes mutated per patient
    n_unique_genes   = n_distinct(Hugo_Symbol),

    # Tumour Mutation Burden
    tmb_nonsynonymous= round(n() / EXOME_SIZE_MB, 4),

    .groups = "drop"
  )

cat(sprintf("mutation_counts: %d patients, %d columns\n",
            nrow(mutation_counts), ncol(mutation_counts)))

cat("\n--- mutation_counts: head(5) ---\n")
print(head(mutation_counts, 5))

cat("\n--- mutation_counts: str() ---\n")
str(mutation_counts)

cat("\n--- TMB summary (mut/Mb) ---\n")
print(summary(mutation_counts$tmb_nonsynonymous))

cat("\n--- Total coding mutations per patient ---\n")
print(summary(mutation_counts$n_total_coding))

n_hyper <- sum(mutation_counts$tmb_nonsynonymous > 10)
cat(sprintf("\nHypermutated (TMB > 10): %d patients (%.1f%%)\n",
            n_hyper,
            n_hyper/nrow(mutation_counts)*100))
cat(sprintf("Median mutations/patient : %.0f\n",
            median(mutation_counts$n_total_coding)))
cat(sprintf("Mean   mutations/patient : %.1f\n",
            mean(mutation_counts$n_total_coding)))

write_csv(mutation_counts,
          file.path(PATHS$out_brca_tables,
                    "03_mutation_counts_per_patient.csv"))
cat(sprintf("Saved: %s/03_mutation_counts_per_patient.csv\n",
            PATHS$out_brca_tables))


# ═══════════════════════════════════════════════════════════════
# SECTION 6 — Top mutated genes
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 6: Top mutated genes\n")
cat("──────────────────────────────\n")

n_patients_total <- n_distinct(maf_coding$patient_id)

top_genes <- maf_coding |>
  group_by(Hugo_Symbol) |>
  summarise(
    n_mutations = n(),
    n_patients  = n_distinct(patient_id),
    pct_patients= round(n_distinct(patient_id)/
                        n_patients_total * 100, 1),
    .groups     = "drop"
  ) |>
  arrange(desc(n_patients))

cat(sprintf("Total unique genes mutated: %d\n",
            nrow(top_genes)))

cat("\n--- Top 30 genes by patients mutated ---\n")
print(head(top_genes, 30))

cat("\n--- Genes in >= 5%  of patients ---\n")
print(top_genes |> filter(pct_patients >= 5))

cat("\n--- Genes in >= 10% of patients ---\n")
print(top_genes |> filter(pct_patients >= 10))

write_csv(top_genes,
          file.path(PATHS$out_brca_tables,
                    "03_top_mutated_genes.csv"))
cat(sprintf("Saved: %s/03_top_mutated_genes.csv\n",
            PATHS$out_brca_tables))

# Select top 20 for feature flags in Script 07
top20_genes <- head(top_genes$Hugo_Symbol, 20)
cat(sprintf("\nTop 20 genes (will become binary features):\n"))
cat(paste(top20_genes, collapse=", "), "\n")


# ═══════════════════════════════════════════════════════════════
# SECTION 7 — Build genomics_cleaned
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 7: Building genomics_cleaned\n")
cat("──────────────────────────────────────\n")

# Select key columns for the cleaned genomics file
genomics_cleaned <- maf_coding |>
  select(
    patient_id,
    Hugo_Symbol,
    Variant_Classification,
    Variant_Type,
    any_of("HGVSp_Short"),
    any_of("HGVSc"),
    any_of("t_depth"),
    any_of("t_alt_count"),
    any_of("t_ref_count"),
    Tumor_Sample_Barcode
  )

# Add VAF column if depth columns are present
if (all(c("t_alt_count","t_depth") %in% colnames(genomics_cleaned))) {
  genomics_cleaned <- genomics_cleaned |>
    mutate(
      vaf = round(as.numeric(t_alt_count) /
                  as.numeric(t_depth), 4)
    )
  cat("VAF computed: t_alt_count / t_depth\n")
  cat("--- VAF summary ---\n")
  print(summary(genomics_cleaned$vaf))
} else {
  cat("VAF not computed: t_alt_count or t_depth not in MAF\n")
}

cat(sprintf("\ngenomics_cleaned: %d mutations, %d columns\n",
            nrow(genomics_cleaned), ncol(genomics_cleaned)))

cat("\n--- genomics_cleaned: head(5) ---\n")
print(head(genomics_cleaned, 5))

cat("\n--- genomics_cleaned: str() ---\n")
str(genomics_cleaned)

write_csv(genomics_cleaned,
          file.path(PATHS$proc_brca, "genomics_cleaned.csv"))
cat(sprintf("Saved: %s/genomics_cleaned.csv\n",
            PATHS$proc_brca))


# ═══════════════════════════════════════════════════════════════
# SECTION 8 — Plots (print to screen THEN save)
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 8: Generating and printing plots\n")
cat("─────────────────────────────────────────\n")

# ── PLOT 1: Variant classification (coding vs non-coding) ─────
p1 <- var_class_summary |>
  mutate(
    Variant_Classification = fct_reorder(
      Variant_Classification, n_mutations)
  ) |>
  ggplot(aes(x=Variant_Classification, y=n_mutations,
             fill=decision)) +
  geom_col(width=0.75, colour="white") +
  geom_text(aes(label=scales::comma(n_mutations)),
            hjust=-0.1, size=3) +
  coord_flip() +
  scale_fill_manual(
    values=c("KEPT (coding)"    = "#185FA5",
             "REMOVED (non-coding)" = "#9CA3AF")) +
  scale_y_continuous(
    labels=scales::comma,
    expand=expansion(mult=c(0, 0.25))) +
  labs(
    title   = "Variant Classification — TCGA-BRCA",
    subtitle= sprintf(
      "Total: %s | Kept (coding): %s (%.1f%%) | Removed: %s (%.1f%%)",
      scales::comma(nrow(maf_raw)),
      scales::comma(nrow(maf_coding)),
      nrow(maf_coding)/nrow(maf_raw)*100,
      scales::comma(nrow(maf_raw)-nrow(maf_coding)),
      (nrow(maf_raw)-nrow(maf_coding))/nrow(maf_raw)*100),
    x       = NULL,
    y       = "Number of Mutations",
    fill    = NULL,
    caption = "Blue = coding (protein-changing) | Grey = non-coding (excluded)"
  ) +
  theme_minimal(base_size=11) +
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="bottom")

print(p1)
ggsave(file.path(PATHS$out_brca_figures,
                 "fig03_variant_classification.png"),
       p1, width=11, height=8, dpi=300)
cat("Printed and saved: fig03_variant_classification.png\n")

# ── PLOT 2: Coding mutations per patient ──────────────────────
# Load clinical for stage colour
clinical_cleaned <- read_csv(
  file.path(PATHS$proc_brca, "clinical_cleaned.csv"),
  show_col_types=FALSE)

p2 <- mutation_counts |>
  left_join(clinical_cleaned |>
              select(patient_id, stage_clean),
            by="patient_id") |>
  mutate(stage_clean = factor(stage_clean,
    levels=c("Stage I","Stage II","Stage III",
             "Stage IV","Other/Unknown"))) |>
  ggplot(aes(x=log10(n_total_coding + 1),
             fill=stage_clean)) +
  geom_histogram(bins=50, colour="white",
                 linewidth=0.1, alpha=0.85) +
  scale_fill_manual(
    values=c("Stage I"  ="#5DCAA5",
             "Stage II" ="#185FA5",
             "Stage III"="#F59E0B",
             "Stage IV" ="#EF4444",
             "Other/Unknown"="#9CA3AF"),
    na.value="#9CA3AF") +
  scale_x_continuous(
    breaks=0:4,
    labels=c("1","10","100","1,000","10,000")) +
  labs(
    title   = "Coding Mutations per Patient — TCGA-BRCA",
    subtitle= sprintf(
      "n=%d patients | Median=%d | Mean=%.0f | log10 scale",
      nrow(mutation_counts),
      median(mutation_counts$n_total_coding),
      mean(mutation_counts$n_total_coding)),
    x       = "Number of Coding Mutations (log10 scale)",
    y       = "Number of Patients",
    fill    = "Cancer Stage",
    caption = "Silent mutations excluded | Source: cBioPortal/GDC"
  ) +
  theme_minimal(base_size=12) +
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="right")

print(p2)
ggsave(file.path(PATHS$out_brca_figures,
                 "fig03_mutations_per_patient.png"),
       p2, width=11, height=6, dpi=300)
cat("Printed and saved: fig03_mutations_per_patient.png\n")

# ── PLOT 3: Top 20 mutated genes ──────────────────────────────
p3 <- top_genes |>
  head(20) |>
  mutate(Hugo_Symbol = fct_reorder(Hugo_Symbol, pct_patients)) |>
  ggplot(aes(x=Hugo_Symbol, y=pct_patients)) +
  geom_col(fill="#185FA5", width=0.75, colour="white") +
  geom_text(aes(label=sprintf("%.1f%%", pct_patients)),
            hjust=-0.1, size=3.5, fontface="bold") +
  coord_flip() +
  scale_y_continuous(
    limits=c(0, max(top_genes$pct_patients[1:20])*1.25),
    labels=function(x) paste0(x,"%")) +
  labs(
    title   = "Top 20 Gene Mutation Frequency — TCGA-BRCA",
    subtitle= "% patients with at least one coding mutation in each gene",
    x       = NULL,
    y       = "% Patients Mutated",
    caption =
      "All 20 genes included as binary features (1=mutated, 0=not mutated)"
  ) +
  theme_minimal(base_size=12) +
  theme(plot.background=element_rect(fill="white",colour=NA))

print(p3)
ggsave(file.path(PATHS$out_brca_figures,
                 "fig03_top_mutated_genes.png"),
       p3, width=10, height=8, dpi=300)
cat("Printed and saved: fig03_top_mutated_genes.png\n")

# ── PLOT 4: TMB distribution ──────────────────────────────────
p4 <- mutation_counts |>
  ggplot(aes(x=tmb_nonsynonymous)) +
  geom_histogram(bins=60, fill="#185FA5", colour="white",
                 linewidth=0.1) +
  geom_vline(xintercept=10, linetype="dashed",
             colour="#EF4444", linewidth=1.1) +
  annotate("text", x=11, y=Inf,
           label=sprintf("TMB > 10\n(n=%d, %.1f%%)\nFDA threshold",
                         n_hyper,
                         n_hyper/nrow(mutation_counts)*100),
           vjust=2, hjust=0, size=3.5, colour="#EF4444") +
  scale_x_continuous(trans="log1p",
    breaks=c(0,1,2,5,10,25,50,100,150),
    labels=c("0","1","2","5","10","25","50","100","150")) +
  labs(
    title   = "Tumour Mutation Burden — TCGA-BRCA",
    subtitle= sprintf(
      "n=%d | Median=%.2f | Mean=%.2f | log1p scale",
      nrow(mutation_counts),
      median(mutation_counts$tmb_nonsynonymous),
      mean(mutation_counts$tmb_nonsynonymous)),
    x       = "TMB (mut/Mb, log1p scale)",
    y       = "Number of Patients",
    caption = sprintf(
      "Exome size = %d Mb | Red dashed = FDA immunotherapy threshold (10 mut/Mb)",
      EXOME_SIZE_MB)
  ) +
  theme_minimal(base_size=12) +
  theme(plot.background=element_rect(fill="white",colour=NA))

print(p4)
ggsave(file.path(PATHS$out_brca_figures,
                 "fig03_tmb_distribution.png"),
       p4, width=10, height=6, dpi=300)
cat("Printed and saved: fig03_tmb_distribution.png\n")


# ═══════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════

cat("\n====================================================\n")
cat(" SCRIPT 03 COMPLETE\n")
cat("====================================================\n")
cat(sprintf("  Raw mutations            : %d\n", nrow(maf_raw)))
cat(sprintf("  Coding mutations (kept)  : %d (%.1f%%)\n",
            nrow(maf_coding),
            nrow(maf_coding)/nrow(maf_raw)*100))
cat(sprintf("  Non-coding (removed)     : %d (%.1f%%)\n",
            nrow(maf_raw)-nrow(maf_coding),
            (nrow(maf_raw)-nrow(maf_coding))/nrow(maf_raw)*100))
cat(sprintf("  Unique full barcodes     : %d (maftools count)\n",
            n_unique_barcodes))
cat(sprintf("  Unique patient IDs       : %d (12-char extract)\n",
            n_unique_patients))
cat(sprintf("  Discrepancy explained    : multi-sample patients\n"))
cat(sprintf("  Unique genes mutated     : %d\n", nrow(top_genes)))
cat(sprintf("  Median mutations/patient : %d\n",
            as.integer(median(mutation_counts$n_total_coding))))
cat(sprintf("  Mean TMB (mut/Mb)        : %.2f\n",
            mean(mutation_counts$tmb_nonsynonymous)))
cat(sprintf("  Median TMB (mut/Mb)      : %.2f\n",
            median(mutation_counts$tmb_nonsynonymous)))
cat(sprintf("  Hypermutated (TMB>10)    : %d (%.1f%%)\n",
            n_hyper, n_hyper/nrow(mutation_counts)*100))
cat(sprintf("  Top 2 genes              : %s (%.1f%%), %s (%.1f%%)\n",
            top_genes$Hugo_Symbol[1], top_genes$pct_patients[1],
            top_genes$Hugo_Symbol[2], top_genes$pct_patients[2]))
cat("\n  Files saved:\n")
cat(sprintf("    %s/genomics_cleaned.csv\n",  PATHS$proc_brca))
cat(sprintf("    %s/03_variant_classification_summary.csv\n",
            PATHS$out_brca_tables))
cat(sprintf("    %s/03_top_mutated_genes.csv\n",
            PATHS$out_brca_tables))
cat(sprintf("    %s/03_mutation_counts_per_patient.csv\n",
            PATHS$out_brca_tables))
cat("\n  Plots printed + saved:\n")
cat(sprintf("    fig03_variant_classification.png\n"))
cat(sprintf("    fig03_mutations_per_patient.png\n"))
cat(sprintf("    fig03_top_mutated_genes.png\n"))
cat(sprintf("    fig03_tmb_distribution.png\n"))
cat("\nNext: source('scripts/04_match_patient_ids.R')\n")
cat("====================================================\n")

