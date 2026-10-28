# ============================================================
# Script 09a: EDA — CPTAC Breast Cancer 2020 (External Dataset)
# Project : Clinical-Genomics Harmonisation Pipeline
# Author  : Muteeba Azhar | MOI Internship | July 2026
#
# PURPOSE:
#   Before applying our trained model to the CPTAC dataset,
#   we must thoroughly understand:
#     1. What clinical variables are available?
#     2. Is OS_STATUS + OS_MONTHS present? (needed for outcome)
#     3. What is the mutation landscape?
#     4. Are our 20 training genes present in CPTAC?
#     5. How does CPTAC compare to TCGA-BRCA?
#     6. Is external validation feasible?
#
# THIS SCRIPT DOES NOT MODIFY ANY DATA.
# It is purely exploratory — all findings feed into Script 09b
# which applies the trained model.
#
# INPUT:
#   data/raw/TEST_brca_cptac_2020/data_clinical_patient.txt
#   data/raw/TEST_brca_cptac_2020/data_clinical_sample.txt
#   data/raw/TEST_brca_cptac_2020/data_mutations.txt
#
# ALSO LOADS (for comparison):
#   data/processed/TCGA-BRCA/clinical_genomics_merged.csv
#   outputs/TCGA-BRCA/tables/03_top_mutated_genes.csv
#
# OUTPUTS:
#   outputs/TCGA-BRCA/tables/09a_cptac_clinical_columns.csv
#   outputs/TCGA-BRCA/tables/09a_cptac_vs_tcga_comparison.csv
#   outputs/TCGA-BRCA/tables/09a_cptac_gene_overlap.csv
#   outputs/TCGA-BRCA/tables/09a_cptac_outcome_feasibility.csv
#   outputs/TCGA-BRCA/figures/fig09a_cptac_clinical_overview.png
#   outputs/TCGA-BRCA/figures/fig09a_cptac_mutation_landscape.png
#   outputs/TCGA-BRCA/figures/fig09a_cptac_vs_tcga_genes.png
#   outputs/TCGA-BRCA/figures/fig09a_cptac_os_distribution.png
# ============================================================

library(tidyverse)
library(scales)
library(maftools)

setwd("E:/R_Studio/Projects/MOI/01R_clinical_genomics")

PATHS <- list(
  raw_cptac        = "data/raw/TEST_brca_cptac_2020",
  proc_brca        = "data/processed/TCGA-BRCA",
  out_brca_tables  = "outputs/TCGA-BRCA/tables",
  out_brca_figures = "outputs/TCGA-BRCA/figures"
)

NA_STRINGS <- c("","NA","N/A","n/a",
                "[Not Available]","[Not Applicable]",
                "[Unknown]","[Discrepancy]")

read_cbio <- function(path, label="") {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path))
  }
  dat <- read_tsv(path, comment="#", na=NA_STRINGS,
                  show_col_types=FALSE, progress=FALSE)
  if (label != "")
    cat(sprintf("  Loaded %-30s → %d rows × %d cols\n",
                label, nrow(dat), ncol(dat)))
  dat
}

# Top 20 training genes
TOP20_GENES <- c(
  "TP53","PIK3CA","TTN","CDH1","GATA3","MUC16","KMT2C",
  "MAP3K1","HMCN1","FLG","RYR2","SYNE1","PTEN","USH2A",
  "NCOR1","ZFHX4","SPTA1","DMD","MUC17","NEB"
)

cat("====================================================\n")
cat(" SCRIPT 09a: EDA — CPTAC Breast Cancer 2020\n")
cat(" External validation candidate dataset\n")
cat("====================================================\n\n")


# ═══════════════════════════════════════════════════════════════
# SECTION 1 — Load all 3 CPTAC files
# ═══════════════════════════════════════════════════════════════

cat("SECTION 1: Loading CPTAC files\n")
cat("────────────────────────────────\n")

# Check files exist
files_needed <- c("data_clinical_patient.txt",
                  "data_clinical_sample.txt",
                  "data_mutations.txt")

cat("Checking files in:", PATHS$raw_cptac, "\n")
for (f in files_needed) {
  path <- file.path(PATHS$raw_cptac, f)
  exists_flag <- file.exists(path)
  size_kb <- if(exists_flag)
    round(file.size(path)/1024, 1) else NA
  cat(sprintf("  %-35s : %s  (%.1f KB)\n",
              f,
              if(exists_flag) "FOUND ✓" else "MISSING ✗",
              ifelse(is.na(size_kb), 0, size_kb)))
}
cat("\n")

cptac_clin <- read_cbio(
  file.path(PATHS$raw_cptac, "data_clinical_patient.txt"),
  "data_clinical_patient.txt")

cptac_samp <- read_cbio(
  file.path(PATHS$raw_cptac, "data_clinical_sample.txt"),
  "data_clinical_sample.txt")

cptac_maf  <- read_cbio(
  file.path(PATHS$raw_cptac, "data_mutations.txt"),
  "data_mutations.txt")


# ═══════════════════════════════════════════════════════════════
# SECTION 2 — Clinical patient file deep inspection
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 2: Clinical file — deep inspection\n")
cat("────────────────────────────────────────────\n")

cat("\n--- CPTAC clinical: all column names ---\n")
print(colnames(cptac_clin))

cat("\n--- CPTAC clinical: str() ---\n")
str(cptac_clin)

cat("\n--- CPTAC clinical: head(5) ---\n")
print(head(cptac_clin, 5))

# Save column list
clin_cols_df <- tibble(
  column_name = colnames(cptac_clin),
  data_type   = sapply(cptac_clin, class),
  n_missing   = sapply(cptac_clin,
                       function(x) sum(is.na(x))),
  pct_missing = round(sapply(cptac_clin,
                              function(x) mean(is.na(x)))*100, 1),
  n_unique    = sapply(cptac_clin, n_distinct),
  example_value = sapply(cptac_clin,
                          function(x) as.character(
                            head(na.omit(x), 1)))
)

cat("\n--- All clinical columns with completeness ---\n")
print(clin_cols_df, n=Inf)

write_csv(clin_cols_df,
          file.path(PATHS$out_brca_tables,
                    "09a_cptac_clinical_columns.csv"))
cat(sprintf("\nSaved: 09a_cptac_clinical_columns.csv\n"))


# ═══════════════════════════════════════════════════════════════
# SECTION 3 — KEY CHECK: Survival variables
# This determines whether external validation is possible
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 3: KEY CHECK — Survival variable availability\n")
cat("───────────────────────────────────────────────────────\n")
cat("For our pipeline we need: OS_STATUS + OS_MONTHS\n")
cat("These are required to compute five_year_survival.\n\n")

# Check for OS variables (may have different names)
all_cols_upper <- str_to_upper(colnames(cptac_clin))

check_vars <- c("OS_STATUS","OS_MONTHS","OS_MONTH",
                "OVERALL_SURVIVAL","SURVIVAL_STATUS",
                "SURVIVAL_MONTHS","AGE","PATH_STAGE",
                "STAGE","TUMOR_STAGE","AJCC","TMB",
                "DFS_STATUS","DFS_MONTHS","PFS_STATUS")

cat("--- Checking for key variables ---\n")
feasibility <- tibble(
  Variable_needed    = check_vars,
  Found_as           = NA_character_,
  Values_example     = NA_character_,
  N_available        = NA_integer_,
  Pct_available      = NA_real_
)

for (i in seq_along(check_vars)) {
  var <- check_vars[i]
  # Exact match first
  match_exact <- colnames(cptac_clin)[
    str_to_upper(colnames(cptac_clin)) == var]
  # Partial match
  match_partial <- colnames(cptac_clin)[
    str_detect(str_to_upper(colnames(cptac_clin)), var)]

  if (length(match_exact) > 0) {
    col <- match_exact[1]
    feasibility$Found_as[i] <- col
    vals <- cptac_clin[[col]]
    feasibility$N_available[i] <- sum(!is.na(vals))
    feasibility$Pct_available[i] <- round(mean(!is.na(vals))*100,1)
    feasibility$Values_example[i] <- paste(
      head(unique(na.omit(as.character(vals))), 3),
      collapse=" | ")
  } else if (length(match_partial) > 0) {
    col <- match_partial[1]
    feasibility$Found_as[i] <- paste0("(partial) ", col)
    vals <- cptac_clin[[col]]
    feasibility$N_available[i] <- sum(!is.na(vals))
    feasibility$Pct_available[i] <- round(mean(!is.na(vals))*100,1)
    feasibility$Values_example[i] <- paste(
      head(unique(na.omit(as.character(vals))), 3),
      collapse=" | ")
  } else {
    feasibility$Found_as[i] <- "NOT FOUND"
    feasibility$N_available[i] <- 0
    feasibility$Pct_available[i] <- 0
  }
}

cat("\n--- Variable availability check ---\n")
print(feasibility |>
        filter(N_available > 0 |
                 Variable_needed %in%
                 c("OS_STATUS","OS_MONTHS","AGE","PATH_STAGE")),
      n=Inf, width=Inf)

# VERDICT
has_os_status <- any(feasibility$N_available[
  feasibility$Variable_needed=="OS_STATUS"] > 0)
has_os_months <- any(feasibility$N_available[
  feasibility$Variable_needed %in% c("OS_MONTHS","OS_MONTH")] > 0)
has_age       <- any(feasibility$N_available[
  feasibility$Variable_needed=="AGE"] > 0)
has_stage     <- any(feasibility$N_available[
  feasibility$Variable_needed %in%
    c("PATH_STAGE","STAGE","TUMOR_STAGE","AJCC")] > 0)

cat("\n══ FEASIBILITY VERDICT ══\n")
cat(sprintf("  OS_STATUS available  : %s\n",
            if(has_os_status) "YES ✓" else "NO ✗"))
cat(sprintf("  OS_MONTHS available  : %s\n",
            if(has_os_months) "YES ✓" else "NO ✗"))
cat(sprintf("  AGE available        : %s\n",
            if(has_age) "YES ✓" else "NO ✗"))
cat(sprintf("  STAGE available      : %s\n",
            if(has_stage) "YES ✓" else "NO ✗"))

if (has_os_status & has_os_months) {
  cat("\n✅ FULL external validation POSSIBLE\n")
  cat("   five_year_survival can be computed.\n")
  cat("   Proceed with Script 09b.\n")
} else if (has_os_status | has_os_months) {
  cat("\n⚠️  PARTIAL survival data available\n")
  cat("   May need to adapt outcome definition.\n")
  cat("   Check what survival columns exist.\n")
} else {
  cat("\n❌ OS_STATUS / OS_MONTHS NOT FOUND\n")
  cat("   Alternative: apply model → predict probabilities\n")
  cat("   Cannot compute AUC (no true outcome to compare)\n")
  cat("   Can still: show predicted risk distribution,\n")
  cat("   gene flag comparison, feature distributions.\n")
}

write_csv(feasibility,
          file.path(PATHS$out_brca_tables,
                    "09a_cptac_outcome_feasibility.csv"))
cat(sprintf("\nSaved: 09a_cptac_outcome_feasibility.csv\n"))


# ═══════════════════════════════════════════════════════════════
# SECTION 4 — Sample file inspection
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 4: Sample file inspection\n")
cat("────────────────────────────────────\n")

cat("--- CPTAC sample: column names ---\n")
print(colnames(cptac_samp))

cat("\n--- CPTAC sample: head(5) ---\n")
print(head(cptac_samp, 5))

cat("\n--- Sample type distribution ---\n")
if ("SAMPLE_TYPE" %in% colnames(cptac_samp)) {
  print(table(cptac_samp$SAMPLE_TYPE, useNA="always"))
}

# TMB in sample file?
tmb_col <- colnames(cptac_samp)[
  str_detect(str_to_upper(colnames(cptac_samp)), "TMB")]
if (length(tmb_col) > 0) {
  cat(sprintf("\nTMB column found: %s\n", tmb_col[1]))
  print(summary(cptac_samp[[tmb_col[1]]]))
}


# ═══════════════════════════════════════════════════════════════
# SECTION 5 — Mutation file inspection
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 5: Mutation file inspection\n")
cat("──────────────────────────────────────\n")

cat(sprintf("Total mutations: %d\n", nrow(cptac_maf)))

cat("\n--- MAF column names ---\n")
print(colnames(cptac_maf))

cat("\n--- MAF: head(3) ---\n")
print(head(cptac_maf, 3))

cat("\n--- Unique patients in MAF ---\n")
n_patients_maf <- n_distinct(
  str_sub(cptac_maf$Tumor_Sample_Barcode, 1, 12))
cat(sprintf("Unique barcodes    : %d\n",
            n_distinct(cptac_maf$Tumor_Sample_Barcode)))
cat(sprintf("Unique patient IDs : %d\n", n_patients_maf))

cat("\n--- Variant classifications ---\n")
print(sort(table(cptac_maf$Variant_Classification),
           decreasing=TRUE))

# Coding mutations
CODING_TYPES <- c("Missense_Mutation","Nonsense_Mutation",
                  "Frame_Shift_Del","Frame_Shift_Ins",
                  "Splice_Site","In_Frame_Del","In_Frame_Ins",
                  "Nonstop_Mutation","Translation_Start_Site")

n_coding <- sum(cptac_maf$Variant_Classification %in%
                  CODING_TYPES)
cat(sprintf("\nCoding mutations   : %d (%.1f%%)\n",
            n_coding, n_coding/nrow(cptac_maf)*100))

# Top mutated genes in CPTAC
cptac_coding <- cptac_maf |>
  filter(Variant_Classification %in% CODING_TYPES) |>
  mutate(patient_id=str_sub(Tumor_Sample_Barcode,1,12))

n_cptac_pts <- n_distinct(cptac_coding$patient_id)

cptac_top_genes <- cptac_coding |>
  group_by(Hugo_Symbol) |>
  summarise(
    n_mutations = n(),
    n_patients  = n_distinct(patient_id),
    pct_patients= round(n_distinct(patient_id)/
                        n_cptac_pts*100, 1),
    .groups     = "drop"
  ) |>
  arrange(desc(n_patients))

cat(sprintf("\nTotal unique genes mutated: %d\n",
            nrow(cptac_top_genes)))

cat("\n--- Top 25 genes in CPTAC ---\n")
print(head(cptac_top_genes, 25))


# ═══════════════════════════════════════════════════════════════
# SECTION 6 — Gene overlap: our 20 training genes in CPTAC
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 6: Gene overlap — training genes in CPTAC\n")
cat("────────────────────────────────────────────────────\n")

gene_overlap <- tibble(Gene=TOP20_GENES) |>
  left_join(
    cptac_top_genes |> select(Hugo_Symbol,
                               cptac_n_patients=n_patients,
                               cptac_pct=pct_patients),
    by=c("Gene"="Hugo_Symbol")
  ) |>
  left_join(
    read_csv(file.path(PATHS$out_brca_tables,
                       "03_top_mutated_genes.csv"),
             show_col_types=FALSE) |>
      select(Hugo_Symbol,
             tcga_n_patients=n_patients,
             tcga_pct=pct_patients),
    by=c("Gene"="Hugo_Symbol")
  ) |>
  mutate(
    cptac_n_patients = coalesce(cptac_n_patients, 0L),
    cptac_pct        = coalesce(cptac_pct, 0),
    present_in_cptac = cptac_n_patients > 0
  )

cat("--- Training gene overlap with CPTAC ---\n")
cat("Gene           TCGA%    CPTAC%   In CPTAC?\n")
cat(strrep("-",50),"\n")
for (i in seq_len(nrow(gene_overlap))) {
  r <- gene_overlap[i,]
  cat(sprintf("%-15s  %5.1f%%   %5.1f%%   %s\n",
              r$Gene, r$tcga_pct, r$cptac_pct,
              if(r$present_in_cptac) "YES" else "NO (flag=0)"))
}
cat(sprintf("\nGenes present in CPTAC: %d / %d\n",
            sum(gene_overlap$present_in_cptac), 20))

write_csv(gene_overlap,
          file.path(PATHS$out_brca_tables,
                    "09a_cptac_gene_overlap.csv"))
cat("Saved: 09a_cptac_gene_overlap.csv\n")


# ═══════════════════════════════════════════════════════════════
# SECTION 7 — TCGA vs CPTAC comparison table
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 7: TCGA-BRCA vs CPTAC-BRCA comparison\n")
cat("──────────────────────────────────────────────────\n")

tcga_merged <- read_csv(
  file.path(PATHS$proc_brca, "clinical_genomics_merged.csv"),
  show_col_types=FALSE)

cptac_tmb <- cptac_coding |>
  group_by(patient_id) |>
  summarise(n_coding=n(),
            tmb=round(n()/38,4),
            .groups="drop")

comparison <- tibble(
  Metric = c(
    "N patients (clinical)",
    "N patients (with mutations)",
    "Total coding mutations",
    "Mean mutations/patient",
    "Median mutations/patient",
    "Mean TMB (mut/Mb)",
    "Median TMB (mut/Mb)",
    "Top gene (% patients)",
    "Second gene (% patients)",
    "AGE available",
    "OS_STATUS available",
    "OS_MONTHS available",
    "STAGE available"
  ),
  TCGA_BRCA = c(
    "1,098",
    "964",
    scales::comma(sum(tcga_merged$n_total_coding, na.rm=TRUE)),
    sprintf("%.1f", mean(tcga_merged$n_total_coding, na.rm=TRUE)),
    sprintf("%.0f", median(tcga_merged$n_total_coding, na.rm=TRUE)),
    sprintf("%.3f", mean(tcga_merged$tmb_nonsynonymous, na.rm=TRUE)),
    sprintf("%.3f", median(tcga_merged$tmb_nonsynonymous, na.rm=TRUE)),
    sprintf("PIK3CA / TP53 (34.3%%)"),
    sprintf("PIK3CA (34.0%%)"),
    "YES (98.5%)",
    "YES (99.9%)",
    "YES (99.9%)",
    "YES (98.9%)"
  ),
  CPTAC_BRCA_2020 = c(
    as.character(nrow(cptac_clin)),
    as.character(n_cptac_pts),
    scales::comma(n_coding),
    sprintf("%.1f", mean(cptac_tmb$n_coding)),
    sprintf("%.0f", median(cptac_tmb$n_coding)),
    sprintf("%.3f", mean(cptac_tmb$tmb)),
    sprintf("%.3f", median(cptac_tmb$tmb)),
    sprintf("%s (%.1f%%)",
            cptac_top_genes$Hugo_Symbol[1],
            cptac_top_genes$pct_patients[1]),
    sprintf("%s (%.1f%%)",
            cptac_top_genes$Hugo_Symbol[2],
            cptac_top_genes$pct_patients[2]),
    if(has_age) "YES" else "NO / CHECK",
    if(has_os_status) "YES" else "NO / CHECK",
    if(has_os_months) "YES" else "NO / CHECK",
    if(has_stage) "YES" else "NO / CHECK"
  )
)

cat("--- TCGA-BRCA vs CPTAC-BRCA 2020 ---\n")
print(comparison, n=Inf, width=Inf)

write_csv(comparison,
          file.path(PATHS$out_brca_tables,
                    "09a_cptac_vs_tcga_comparison.csv"))
cat("Saved: 09a_cptac_vs_tcga_comparison.csv\n")


# ═══════════════════════════════════════════════════════════════
# SECTION 8 — Plots (printed then saved)
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 8: EDA Plots\n")
cat("──────────────────────\n")

# ── PLOT 1: CPTAC clinical overview ──────────────────────────
# Age distribution (if available)
age_col <- colnames(cptac_clin)[
  str_to_upper(colnames(cptac_clin)) %in% c("AGE","AGE_AT_DIAGNOSIS")]

if (length(age_col) > 0) {
  age_data <- cptac_clin[[age_col[1]]] |>
    as.numeric() |> na.omit()

  p1 <- ggplot(tibble(age=age_data), aes(x=age)) +
    geom_histogram(bins=25, fill="#185FA5",
                   colour="white", linewidth=0.2) +
    geom_vline(xintercept=mean(age_data),
               linetype="dashed", colour="#EF4444",
               linewidth=1) +
    annotate("text", x=mean(age_data)+2, y=Inf,
             label=sprintf("Mean=%.1f", mean(age_data)),
             vjust=2, hjust=0, size=3.5, colour="#EF4444") +
    labs(title="CPTAC-BRCA 2020 — Age Distribution",
         subtitle=sprintf("n=%d | Mean=%.1f | Median=%.1f",
                          length(age_data),
                          mean(age_data),
                          median(age_data)),
         x="Age at Diagnosis",
         y="Number of Patients",
         caption="Source: CPTAC-BRCA 2020 via cBioPortal") +
    theme_minimal(base_size=12) +
    theme(plot.background=element_rect(fill="white",colour=NA))

  print(p1)
  ggsave(file.path(PATHS$out_brca_figures,
                   "fig09a_cptac_age.png"),
         p1, width=8, height=5, dpi=300)
  cat("Printed and saved: fig09a_cptac_age.png\n")
}

# ── PLOT 2: CPTAC top 20 mutated genes ───────────────────────
p2 <- cptac_top_genes |>
  head(20) |>
  mutate(Hugo_Symbol=fct_reorder(Hugo_Symbol,pct_patients),
         in_training=Hugo_Symbol %in% TOP20_GENES) |>
  ggplot(aes(x=Hugo_Symbol, y=pct_patients,
             fill=in_training)) +
  geom_col(width=0.75, colour="white") +
  geom_text(aes(label=sprintf("%.1f%%", pct_patients)),
            hjust=-0.1, size=3.5, fontface="bold") +
  coord_flip() +
  scale_fill_manual(
    values=c("TRUE"="#0F6E56","FALSE"="#9CA3AF"),
    labels=c("TRUE"="In training gene set (20 genes)",
             "FALSE"="Not in training set")) +
  scale_y_continuous(
    limits=c(0, max(cptac_top_genes$pct_patients[1:20])*1.3),
    labels=function(x) paste0(x,"%")) +
  labs(title="Top 20 Mutated Genes — CPTAC-BRCA 2020",
       subtitle=sprintf("n=%d patients | Green = in our 20 training genes",
                        n_cptac_pts),
       x=NULL, y="% Patients Mutated",
       fill=NULL,
       caption="CPTAC-BRCA 2020 | Coding mutations only") +
  theme_minimal(base_size=11) +
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="bottom")

print(p2)
ggsave(file.path(PATHS$out_brca_figures,
                 "fig09a_cptac_mutation_landscape.png"),
       p2, width=10, height=8, dpi=300)
cat("Printed and saved: fig09a_cptac_mutation_landscape.png\n")

# ── PLOT 3: Gene frequency comparison TCGA vs CPTAC ──────────
tcga_genes <- read_csv(
  file.path(PATHS$out_brca_tables,
            "03_top_mutated_genes.csv"),
  show_col_types=FALSE)

compare_genes_df <- gene_overlap |>
  pivot_longer(
    cols      = c(tcga_pct, cptac_pct),
    names_to  = "dataset",
    values_to = "pct_mutated"
  ) |>
  mutate(
    dataset = recode(dataset,
                     "tcga_pct"  = "TCGA-BRCA (n=964)",
                     "cptac_pct" = "CPTAC-BRCA 2020"),
    Gene = fct_reorder(Gene, pct_mutated, .fun=max)
  )

p3 <- ggplot(compare_genes_df,
             aes(x=Gene, y=pct_mutated, fill=dataset)) +
  geom_col(position="dodge", width=0.7,
           colour="white") +
  coord_flip() +
  scale_fill_manual(
    values=c("TCGA-BRCA (n=964)"="#185FA5",
             "CPTAC-BRCA 2020"="#0F6E56")) +
  scale_y_continuous(labels=function(x) paste0(x,"%")) +
  labs(title="Gene Mutation Frequency: TCGA-BRCA vs CPTAC-BRCA",
       subtitle="Our 20 training genes — % patients mutated in each dataset",
       x=NULL, y="% Patients Mutated",
       fill="Dataset",
       caption="Genes absent in CPTAC will have flag=0 for all CPTAC patients") +
  theme_minimal(base_size=11) +
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="bottom")

print(p3)
ggsave(file.path(PATHS$out_brca_figures,
                 "fig09a_cptac_vs_tcga_genes.png"),
       p3, width=11, height=9, dpi=300)
cat("Printed and saved: fig09a_cptac_vs_tcga_genes.png\n")

# ── PLOT 4: OS_MONTHS if available ───────────────────────────
os_months_col <- colnames(cptac_clin)[
  str_detect(str_to_upper(colnames(cptac_clin)),
             "OS_MONTH|OVERALL_SURVIVAL")]

os_status_col <- colnames(cptac_clin)[
  str_detect(str_to_upper(colnames(cptac_clin)),
             "OS_STATUS|SURVIVAL_STATUS|VITAL")]

if (length(os_months_col) > 0) {
  os_data <- tibble(
    os_months = as.numeric(cptac_clin[[os_months_col[1]]]),
    os_status = if(length(os_status_col)>0)
      cptac_clin[[os_status_col[1]]] else NA
  ) |> filter(!is.na(os_months))

  p4 <- ggplot(os_data, aes(x=os_months)) +
    geom_histogram(bins=35, fill="#185FA5",
                   colour="white", linewidth=0.1) +
    geom_vline(xintercept=60, linetype="dashed",
               colour="#EF4444", linewidth=1.2) +
    annotate("text", x=63, y=Inf,
             label="60 months\n= 5 years",
             vjust=2, hjust=0, size=3.5,
             colour="#EF4444") +
    labs(title="CPTAC-BRCA 2020 — OS Months Distribution",
         subtitle=sprintf(
           "n=%d patients | Patients with >60 months: %d",
           nrow(os_data),
           sum(os_data$os_months >= 60, na.rm=TRUE)),
         x="Overall Survival (months)",
         y="Number of Patients",
         caption="Red dashed = 5-year threshold (60 months)") +
    theme_minimal(base_size=12) +
    theme(plot.background=element_rect(fill="white",colour=NA))

  print(p4)
  ggsave(file.path(PATHS$out_brca_figures,
                   "fig09a_cptac_os_distribution.png"),
         p4, width=9, height=6, dpi=300)
  cat("Printed and saved: fig09a_cptac_os_distribution.png\n")
} else {
  cat("OS_MONTHS not found — skipping OS distribution plot\n")
}


# ═══════════════════════════════════════════════════════════════
# FINAL SUMMARY + DECISION TREE
# ═══════════════════════════════════════════════════════════════

cat("\n====================================================\n")
cat(" SCRIPT 09a COMPLETE — EDA SUMMARY\n")
cat("====================================================\n")
cat(sprintf("  CPTAC clinical patients : %d\n", nrow(cptac_clin)))
cat(sprintf("  CPTAC with mutations    : %d\n", n_cptac_pts))
cat(sprintf("  Total coding mutations  : %d\n", n_coding))
cat(sprintf("  Unique genes mutated    : %d\n", nrow(cptac_top_genes)))
cat(sprintf("  Training genes in CPTAC : %d / 20\n",
            sum(gene_overlap$present_in_cptac)))
cat(sprintf("  OS_STATUS available     : %s\n",
            if(has_os_status) "YES" else "NO"))
cat(sprintf("  OS_MONTHS available     : %s\n",
            if(has_os_months) "YES" else "NO"))
cat(sprintf("  AGE available           : %s\n",
            if(has_age) "YES" else "NO"))
cat(sprintf("  STAGE available         : %s\n",
            if(has_stage) "YES" else "NO"))

cat("\n══ WHAT NEXT (Script 09b) ══\n\n")

if (has_os_status & has_os_months) {
  cat("SCENARIO 1 — FULL EXTERNAL VALIDATION:\n")
  cat("  ✅ Compute five_year_survival on CPTAC\n")
  cat("  ✅ Apply trained model → predict probabilities\n")
  cat("  ✅ Compute AUC, Recall, Specificity, CM\n")
  cat("  ✅ Compare CPTAC AUC vs TCGA test set AUC\n")
} else {
  cat("SCENARIO 2 — PREDICTION ONLY (no outcome available):\n")
  cat("  ✅ Apply trained model → predict probabilities\n")
  cat("  ✅ Show predicted risk distribution in CPTAC\n")
  cat("  ✅ Compare gene flag patterns TCGA vs CPTAC\n")
  cat("  ✅ Report limitations in final report\n")
  cat("\n  This is still scientifically valid:\n")
  cat("  'We applied the model to an independent cohort\n")
  cat("   and report predicted risk scores. Formal AUC\n")
  cat("   evaluation was not possible due to limited\n")
  cat("   survival data in the CPTAC dataset.'\n")
}

cat("\n  Plots saved:\n")
cat("    fig09a_cptac_age.png\n")
cat("    fig09a_cptac_mutation_landscape.png\n")
cat("    fig09a_cptac_vs_tcga_genes.png\n")
cat("    fig09a_cptac_os_distribution.png\n")
cat("\nNext: source('scripts/09b_test_model_on_cptac.R')\n")
cat("====================================================\n")
