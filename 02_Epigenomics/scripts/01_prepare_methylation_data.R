#!/usr/bin/env Rscript
# ============================================================
# 01_prepare_methylation_data.R
# Query GDC for TCGA-BRCA methylation metadata, match against
# Month 1 cohort, save manifest + summary plots.
#
# Run from Ubuntu terminal with:
#   Rscript scripts/01_prepare_methylation_data.R
# ============================================================

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(ggplot2)
})

# ---- 0. Working directory (Ubuntu/WSL path to the same E: drive) ----
setwd("/mnt/e/R_Studio/Projects/MOI/Month_2")
cat("Working directory set to:", getwd(), "\n")

# ---- 1. Confirm Month 1 data is reachable before doing anything else ----
month1_path <- "../01_clinical_genomics/data/processed/clinical_genomics_merged.csv"
if (!file.exists(month1_path)) {
  stop("Cannot find Month 1 merged dataset at: ", month1_path,
       "\nCheck that 01_clinical_genomics/ sits alongside Month_2/ and the path is correct.")
}
cat("Month 1 dataset found.\n")

# ---- 2. Query GDC (metadata only - no download yet) ----
cat("Querying GDC for TCGA-BRCA methylation metadata...\n")
query_meth <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "DNA Methylation",
  platform = "Illumina Human Methylation 450",
  data.type = "Methylation Beta Value"
)
results <- getResults(query_meth)
results$patient_id <- results$cases.submitter_id
cat("Total methylation files found:", nrow(results), "\n")

# ---- 3. Load Month 1 cohort for matching ----
month1 <- read.csv(month1_path)
month1_ids <- unique(month1$patient_id)
cat("Month 1 patients:", length(month1_ids), "\n")

# ---- 4. Matching summary table ----
meth_patient_ids <- unique(results$patient_id)
overlap_all <- intersect(month1_ids, meth_patient_ids)

tumor_only <- results[results$sample_type == "Primary Tumor", ]
overlap_tumor <- intersect(month1_ids, unique(tumor_only$patient_id))

id_matching_summary <- data.frame(
  metric = c("Month1 patients", "Methylation patients (all)",
             "Overlap (all sample types)", "Methylation files (Primary Tumor)",
             "Overlap (Primary Tumor only)"),
  value = c(length(month1_ids), length(meth_patient_ids),
            length(overlap_all), nrow(tumor_only), length(overlap_tumor))
)
write.csv(id_matching_summary, "outputs/tables/methylation_id_matching_summary.csv", row.names = FALSE)
cat("\n===== ID MATCHING SUMMARY =====\n")
print(id_matching_summary)
cat("================================\n\n")

# ---- 5. Plot: sample type breakdown ----
p1 <- ggplot(results, aes(x = sample_type, fill = sample_type)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  labs(title = "TCGA-BRCA Methylation Files by Sample Type",
       x = "Sample Type", y = "Number of Files") +
  theme_minimal() +
  theme(legend.position = "none")
ggsave("outputs/figures/methylation_sample_type_breakdown.png", p1, width = 6, height = 5)
cat("Saved: outputs/figures/methylation_sample_type_breakdown.png\n")

# ---- 6. Plot: overlap with Month 1 cohort ----
overlap_df <- data.frame(
  category = c("Month 1 only\n(no methylation)", "Overlap\n(both)", "Methylation only\n(no Month 1 match)"),
  count = c(length(month1_ids) - length(overlap_all),
            length(overlap_all),
            length(meth_patient_ids) - length(overlap_all))
)
p2 <- ggplot(overlap_df, aes(x = category, y = count, fill = category)) +
  geom_col() +
  geom_text(aes(label = count), vjust = -0.5) +
  labs(title = "Patient Overlap: Month 1 Cohort vs Methylation Data",
       x = "", y = "Number of Patients") +
  theme_minimal() +
  theme(legend.position = "none")
ggsave("outputs/figures/methylation_month1_overlap.png", p2, width = 6, height = 4)
cat("Saved: outputs/figures/methylation_month1_overlap.png\n")

# ---- 7. Save matched sample manifest (metadata only - not beta values yet) ----
matched_manifest <- results[results$patient_id %in% month1_ids &
                             results$sample_type == "Primary Tumor", ]
write.csv(matched_manifest, "data/processed/methylation_matched_manifest.csv", row.names = FALSE)

cat("\nScript 01 complete.\n")
cat("Matched manifest saved:", nrow(matched_manifest), "tumor samples.\n")
cat("-> data/processed/methylation_matched_manifest.csv\n")
