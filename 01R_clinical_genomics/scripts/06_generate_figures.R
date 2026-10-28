# ============================================================
# Script 06: Mutation Landscape Figures
# Project : Clinical-Genomics Harmonisation Pipeline
# Author  : Muteeba Azhar | MOI Internship | July 2026
#
# PURPOSE:
#   Generate 8 publication-quality figures describing the
#   TCGA-BRCA mutation landscape. All plots are printed to
#   the screen before being saved to file.
#
# CAUTIOUS LANGUAGE:
#   Patterns are described as "consistent with" not "confirms".
#   No formal statistical tests are performed in this script.
#   TP53/PIK3CA pattern described as "consistent with mutual
#   exclusivity" — formal test (Fisher's exact) not done.
#   SBS1 interpretation: "consistent with" COSMIC SBS1 —
#   formal signature decomposition not performed.
#
# INPUTS:
#   data/processed/TCGA-BRCA/clinical_genomics_merged.csv
#   data/processed/TCGA-BRCA/genomics_cleaned.csv
#   data/raw/TCGA-BRCA/data_mutations.txt  (maftools)
#   outputs/TCGA-BRCA/tables/03_top_mutated_genes.csv
#   outputs/TCGA-BRCA/tables/03_mutation_counts_per_patient.csv
#
# OUTPUTS (all printed to screen THEN saved):
#   outputs/TCGA-BRCA/figures/fig06_top_mutated_genes.png
#   outputs/TCGA-BRCA/figures/fig06_variant_classification.png
#   outputs/TCGA-BRCA/figures/fig06_mutations_per_patient.png
#   outputs/TCGA-BRCA/figures/fig06_tmb_by_stage.png
#   outputs/TCGA-BRCA/figures/fig06_maf_summary.png     (maftools)
#   outputs/TCGA-BRCA/figures/fig06_oncoplot_top20.png  (maftools)
#   outputs/TCGA-BRCA/figures/fig06_titv.png            (maftools)
#   outputs/TCGA-BRCA/figures/fig06_vaf_distribution.png
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

cat("====================================================\n")
cat(" SCRIPT 06: Mutation Landscape Figures\n")
cat(" All plots: printed to screen THEN saved.\n")
cat(" Language: cautious — 'consistent with' not 'proves'.\n")
cat("====================================================\n\n")


# ═══════════════════════════════════════════════════════════════
# SECTION 1 — Load data
# ═══════════════════════════════════════════════════════════════

cat("SECTION 1: Loading data\n")
cat("────────────────────────\n")

merged <- read_csv(
  file.path(PATHS$proc_brca, "clinical_genomics_merged.csv"),
  show_col_types=FALSE)

genomics_cleaned <- read_csv(
  file.path(PATHS$proc_brca, "genomics_cleaned.csv"),
  show_col_types=FALSE)

top_genes <- read_csv(
  file.path(PATHS$out_brca_tables, "03_top_mutated_genes.csv"),
  show_col_types=FALSE)

mut_counts <- read_csv(
  file.path(PATHS$out_brca_tables,
            "03_mutation_counts_per_patient.csv"),
  show_col_types=FALSE)

var_class <- read_csv(
  file.path(PATHS$out_brca_tables,
            "03_variant_classification_summary.csv"),
  show_col_types=FALSE)

N_MATCHED <- nrow(merged)

cat(sprintf("merged          : %d patients\n", N_MATCHED))
cat(sprintf("genomics_cleaned: %d coding mutations\n",
            nrow(genomics_cleaned)))
cat(sprintf("top_genes       : %d unique genes\n",
            nrow(top_genes)))

# Load raw MAF for maftools (needs original file)
cat("\nLoading raw MAF for maftools...\n")
maf_raw_path <- file.path(PATHS$raw_brca, "data_mutations.txt")

# maftools reads cBioPortal files directly
brca_maf <- tryCatch({
  read.maf(maf = maf_raw_path,
           clinicalData = merged,
           isTCGA = TRUE,
           verbose = FALSE)
}, error = function(e) {
  cat(sprintf("maftools load note: %s\n", e$message))
  cat("maftools figures will be skipped if load fails.\n")
  NULL
})

if (!is.null(brca_maf)) {
  cat(sprintf("maftools MAF loaded: %d samples\n",
              length(unique(getSampleSummary(brca_maf)$Tumor_Sample_Barcode))))
}


# ═══════════════════════════════════════════════════════════════
# SECTION 2 — FIGURE 1: Top 20 mutated genes
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 2: Figure 1 — Top 20 mutated genes\n")
cat("────────────────────────────────────────────\n")

fig1 <- top_genes |>
  head(20) |>
  mutate(Hugo_Symbol = fct_reorder(Hugo_Symbol, pct_patients)) |>
  ggplot(aes(x=Hugo_Symbol, y=pct_patients)) +
  geom_col(fill="#185FA5", width=0.75, colour="white") +
  geom_text(aes(label=sprintf("%.1f%%", pct_patients)),
            hjust=-0.1, size=3.5, fontface="bold",
            colour="#1F2937") +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, max(top_genes$pct_patients[1:20]) * 1.3),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title   = "Top 20 Gene Mutation Frequency — TCGA-BRCA",
    subtitle= sprintf(
      "n=%d matched patients | %% with >=1 coding mutation",
      N_MATCHED),
    x       = NULL,
    y       = "% Patients Mutated",
    caption = paste0(
      "TP53 and PIK3CA show a pattern consistent with mutual exclusivity",
      " (formal Fisher's exact test not performed).\n",
      "All 20 genes included as binary features (1=mutated, 0=not mutated).",
      " Source: cBioPortal / GDC")
  ) +
  theme_minimal(base_size=12) +
  theme(plot.background=element_rect(fill="white",colour=NA),
        plot.caption=element_text(size=8, colour="#6B7280"))

print(fig1)
ggsave(file.path(PATHS$out_brca_figures,
                 "fig06_top_mutated_genes.png"),
       fig1, width=10, height=8, dpi=300)
cat("Printed and saved: fig06_top_mutated_genes.png\n")


# ═══════════════════════════════════════════════════════════════
# SECTION 3 — FIGURE 2: Variant classification breakdown
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 3: Figure 2 — Variant classification\n")
cat("──────────────────────────────────────────────\n")

fig2 <- var_class |>
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
    values=c("KEPT (coding)"        = "#185FA5",
             "REMOVED (non-coding)" = "#9CA3AF")) +
  scale_y_continuous(
    labels=scales::comma,
    expand=expansion(mult=c(0, 0.25))
  ) +
  labs(
    title   = "Variant Classification Breakdown — TCGA-BRCA",
    subtitle= sprintf(
      "Total: %s | Kept (coding): %s (%.1f%%) | Removed: %s",
      scales::comma(sum(var_class$n_mutations)),
      scales::comma(sum(var_class$n_mutations[var_class$is_coding])),
      sum(var_class$n_mutations[var_class$is_coding])/
        sum(var_class$n_mutations)*100,
      scales::comma(sum(var_class$n_mutations[!var_class$is_coding]))),
    x       = NULL,
    y       = "Number of Mutations",
    fill    = NULL,
    caption = "Blue = coding (protein-changing, kept) | Grey = non-coding (excluded). Source: cBioPortal / GDC"
  ) +
  theme_minimal(base_size=11) +
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="bottom")

print(fig2)
ggsave(file.path(PATHS$out_brca_figures,
                 "fig06_variant_classification.png"),
       fig2, width=11, height=8, dpi=300)
cat("Printed and saved: fig06_variant_classification.png\n")


# ═══════════════════════════════════════════════════════════════
# SECTION 4 — FIGURE 3: Coding mutations per patient
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 4: Figure 3 — Mutations per patient\n")
cat("─────────────────────────────────────────────\n")

fig3 <- mut_counts |>
  left_join(merged |> select(patient_id, stage_clean),
            by="patient_id") |>
  mutate(stage_clean = factor(stage_clean,
    levels=c("Stage I","Stage II","Stage III",
             "Stage IV","Other/Unknown"))) |>
  ggplot(aes(x=log10(n_total_coding + 1),
             fill=stage_clean)) +
  geom_histogram(bins=50, colour="white",
                 linewidth=0.1, alpha=0.85) +
  scale_fill_manual(
    values=c("Stage I"      ="#5DCAA5",
             "Stage II"     ="#185FA5",
             "Stage III"    ="#F59E0B",
             "Stage IV"     ="#EF4444",
             "Other/Unknown"="#9CA3AF"),
    na.value="#9CA3AF") +
  scale_x_continuous(
    breaks=0:4,
    labels=c("1","10","100","1,000","10,000")) +
  labs(
    title   = "Coding Mutations per Patient — TCGA-BRCA",
    subtitle= sprintf(
      "n=%d patients | Median=%d | Mean=%.0f | log10 scale",
      nrow(mut_counts),
      median(mut_counts$n_total_coding),
      mean(mut_counts$n_total_coding)),
    x       = "Number of Coding Mutations (log10 scale)",
    y       = "Number of Patients",
    fill    = "Cancer Stage",
    caption = "Silent mutations excluded | Source: cBioPortal / GDC"
  ) +
  theme_minimal(base_size=12) +
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="right")

print(fig3)
ggsave(file.path(PATHS$out_brca_figures,
                 "fig06_mutations_per_patient.png"),
       fig3, width=11, height=6, dpi=300)
cat("Printed and saved: fig06_mutations_per_patient.png\n")


# ═══════════════════════════════════════════════════════════════
# SECTION 5 — FIGURE 4: TMB by cancer stage
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 5: Figure 4 — TMB by stage\n")
cat("─────────────────────────────────────\n")

fig4 <- merged |>
  filter(stage_clean != "Other/Unknown") |>
  mutate(stage_clean = factor(stage_clean,
    levels=c("Stage I","Stage II","Stage III","Stage IV"))) |>
  ggplot(aes(x=stage_clean, y=tmb_nonsynonymous,
             colour=stage_clean)) +
  geom_boxplot(outlier.shape=NA, width=0.5,
               fill="white", linewidth=0.8) +
  geom_jitter(alpha=0.25, width=0.2, size=0.9) +
  geom_hline(yintercept=10, linetype="dashed",
             colour="#EF4444", linewidth=0.9) +
  annotate("text", x=0.55, y=11.5,
           label="TMB=10 (FDA threshold)",
           size=3.2, colour="#EF4444", hjust=0) +
  scale_y_continuous(trans="log1p",
    breaks=c(0,1,2,5,10,25,50,100,150)) +
  scale_colour_manual(
    values=c("Stage I"  ="#5DCAA5",
             "Stage II" ="#185FA5",
             "Stage III"="#F59E0B",
             "Stage IV" ="#EF4444")) +
  labs(
    title   = "Tumour Mutation Burden by Stage — TCGA-BRCA",
    subtitle= "Descriptive only — no formal statistical test performed",
    x       = "Cancer Stage",
    y       = "TMB (mut/Mb, log1p scale)",
    caption = paste0(
      "No clear linear stage-TMB relationship observed in this cohort,\n",
      "which is consistent with published BRCA literature. ",
      "Exome size = 38 Mb | Source: cBioPortal / GDC")
  ) +
  theme_minimal(base_size=12) +
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="none")

print(fig4)
ggsave(file.path(PATHS$out_brca_figures,
                 "fig06_tmb_by_stage.png"),
       fig4, width=9, height=6, dpi=300)
cat("Printed and saved: fig06_tmb_by_stage.png\n")


# ═══════════════════════════════════════════════════════════════
# SECTION 6 — maftools FIGURES 5-7
# MAF Summary dashboard, Oncoplot, Ti/Tv
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 6: maftools figures (5, 6, 7)\n")
cat("───────────────────────────────────────\n")

if (!is.null(brca_maf)) {

  # ── FIGURE 5: MAF Summary dashboard ──────────────────────────
  cat("Generating Figure 5: MAF Summary dashboard...\n")
  png(file.path(PATHS$out_brca_figures, "fig06_maf_summary.png"),
      width=2800, height=2200, res=300)
  plotmafSummary(
    maf         = brca_maf,
    rmOutlier   = TRUE,
    addStat     = "median",
    dashboard   = TRUE,
    titvRaw     = FALSE
  )
  dev.off()
  # Also display in RStudio plots pane
  plotmafSummary(maf=brca_maf, rmOutlier=TRUE,
                 addStat="median", dashboard=TRUE,
                 titvRaw=FALSE)
  cat("Printed and saved: fig06_maf_summary.png\n")

  # ── FIGURE 6: Oncoplot top 20 genes ──────────────────────────
  cat("Generating Figure 6: Oncoplot top 20 genes...\n")
  cat("Note: TP53/PIK3CA pattern is consistent with mutual\n")
  cat("exclusivity — formal co-occurrence test not performed.\n")
  png(file.path(PATHS$out_brca_figures,
                "fig06_oncoplot_top20.png"),
      width=3200, height=2400, res=300)
  oncoplot(
    maf    = brca_maf,
    top    = 20,
    title  = "TCGA-BRCA — Top 20 Mutated Genes",
    showTumorSampleBarcodes=FALSE
  )
  dev.off()
  oncoplot(maf=brca_maf, top=20,
           title="TCGA-BRCA — Top 20 Mutated Genes",
           showTumorSampleBarcodes=FALSE)
  cat("Printed and saved: fig06_oncoplot_top20.png\n")

  # ── FIGURE 7: Ti/Tv ratio ─────────────────────────────────────
  cat("Generating Figure 7: Ti/Tv ratio...\n")
  cat("Note: C>T dominance is consistent with COSMIC SBS1\n")
  cat("(age-related deamination). Formal signature decomposition\n")
  cat("(e.g. deconstructSigs) not performed here.\n")
  titv_data <- titv(maf=brca_maf, plot=FALSE,
                    useSyn=FALSE)
  png(file.path(PATHS$out_brca_figures, "fig06_titv.png"),
      width=2800, height=1800, res=300)
  plotTiTv(res=titv_data)
  dev.off()
  plotTiTv(res=titv_data)
  cat("Printed and saved: fig06_titv.png\n")

} else {
  cat("maftools figures skipped — MAF load failed.\n")
  cat("Check that data_mutations.txt is in data/raw/TCGA-BRCA/\n")
}


# ═══════════════════════════════════════════════════════════════
# SECTION 7 — FIGURE 8: VAF distribution
# CAUTIOUS LANGUAGE: "consistent with" not "confirms"
# ═══════════════════════════════════════════════════════════════

cat("\nSECTION 7: Figure 8 — VAF distribution\n")
cat("─────────────────────────────────────────\n")

if ("vaf" %in% colnames(genomics_cleaned) &&
    sum(!is.na(genomics_cleaned$vaf)) > 0) {

  vaf_data <- genomics_cleaned |>
    filter(!is.na(vaf), vaf >= 0, vaf <= 1)

  cat(sprintf("VAF values available: %d mutations\n",
              nrow(vaf_data)))
  cat("--- VAF summary ---\n")
  print(summary(vaf_data$vaf))

  fig8 <- ggplot(vaf_data, aes(x=vaf)) +
    geom_histogram(bins=50, fill="#185FA5",
                   colour="white", linewidth=0.1) +
    scale_x_continuous(labels=scales::percent_format()) +
    labs(
      title   = "Variant Allele Frequency — TCGA-BRCA",
      subtitle= sprintf(
        "n=%s coding mutations | Mean=%.3f | Median=%.3f",
        scales::comma(nrow(vaf_data)),
        mean(vaf_data$vaf, na.rm=TRUE),
        median(vaf_data$vaf, na.rm=TRUE)),
      x       = "Variant Allele Frequency (VAF)",
      y       = "Number of Mutations",
      caption = paste0(
        "VAF = t_alt_count / t_depth.\n",
        "The peak at 15-25% is consistent with heterozygous ",
        "somatic mutations in bulk tumour sequencing,\n",
        "though formal clonality analysis was not performed. ",
        "Source: cBioPortal / GDC")
    ) +
    theme_minimal(base_size=12) +
    theme(plot.background=element_rect(fill="white",colour=NA),
          plot.caption=element_text(size=8, colour="#6B7280"))

  print(fig8)
  ggsave(file.path(PATHS$out_brca_figures,
                   "fig06_vaf_distribution.png"),
         fig8, width=10, height=6, dpi=300)
  cat("Printed and saved: fig06_vaf_distribution.png\n")

} else {
  cat("VAF not available in genomics_cleaned — skipping.\n")
}


# ═══════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════

cat("\n====================================================\n")
cat(" SCRIPT 06 COMPLETE\n")
cat("====================================================\n")
cat("  All figures printed to screen before saving.\n\n")
cat("  ggplot2 figures:\n")
cat("    fig06_top_mutated_genes.png\n")
cat("    fig06_variant_classification.png\n")
cat("    fig06_mutations_per_patient.png\n")
cat("    fig06_tmb_by_stage.png\n")
cat("    fig06_vaf_distribution.png\n")
cat("\n  maftools figures:\n")
cat("    fig06_maf_summary.png\n")
cat("    fig06_oncoplot_top20.png\n")



# ── MAFTOOLS FIGURES ─────────────────────────────────────────
cat("\n[5] Loading MAF for maftools figures...\n")

maf <- read.maf(
  maf     = file.path(RAW_DIR,"data_mutations.txt"),
  isTCGA  = TRUE,
  useAll  = TRUE,
  verbose = FALSE
)
cat(sprintf("    MAF loaded: %d mutations\n",nrow(maf@data)))


# ── FIGURE 5: MAF summary ────────────────────────────────────
cat("\n[6] Figure 5: MAF summary dashboard...\n")

# maftools base R plot: shows in Plots pane automatically
# dev.copy saves the current Plots pane plot to PNG
plotmafSummary(maf=maf, rmOutlier=TRUE,
               addStat="median", dashboard=TRUE, titvRaw=FALSE)

dev.copy(png, file.path(OUT_F,"fig06_maf_summary.png"),
         width=3200, height=2600, res=300)
dev.off()
cat(sprintf("    Saved: %s/fig06_maf_summary.png\n",OUT_F))


# ── FIGURE 6: Oncoplot ───────────────────────────────────────
cat("\n[7] Figure 6: Oncoplot top 20...\n")

oncoplot(maf=maf, top=20, draw_titv=TRUE,
         titleText=sprintf("TCGA-BRCA — Top 20 Mutated Genes (n=%d)",
                           nrow(getSampleSummary(maf))))

dev.copy(png, file.path(OUT_F,"fig06_oncoplot_top20.png"),
         width=4200, height=3200, res=300)
dev.off()
cat(sprintf("    Saved: %s/fig06_oncoplot_top20.png\n",OUT_F))


# ── FIGURE 7: Ti/Tv ──────────────────────────────────────────
cat("\n[8] Figure 7: Ti/Tv ratio...\n")

titv_data <- titv(maf=maf, plot=FALSE, useSyn=TRUE)
plotTiTv(res=titv_data)

dev.copy(png, file.path(OUT_F,"fig06_titv.png"),
         width=3200, height=2200, res=300)
dev.off()
cat(sprintf("    Saved: %s/fig06_titv.png\n",OUT_F))


# ── FIGURE 8: VAF ────────────────────────────────────────────
cat("\n[9] Figure 8: VAF distribution...\n")

flds <- getFields(maf)
cat(sprintf("    t_AF present      : %s\n","t_AF" %in% flds))
cat(sprintf("    t_alt_count present: %s\n","t_alt_count" %in% flds))

if("t_AF" %in% flds){
  plotVaf(maf=maf, vafCol="t_AF")
  dev.copy(png,file.path(OUT_F,"fig06_vaf.png"),width=3200,height=2200,res=300)
  dev.off()
  cat("    Used: t_AF\n")
} else if(all(c("t_alt_count","t_depth") %in% flds)){
  vaf_df <- maf@data |> as.data.frame() |>
    filter(!is.na(t_alt_count),!is.na(t_depth),t_depth>0) |>
    mutate(vaf=as.numeric(t_alt_count)/as.numeric(t_depth))
  cat(sprintf("    n mutations: %d | Mean VAF: %.3f | Median: %.3f\n",
              nrow(vaf_df),mean(vaf_df$vaf,na.rm=TRUE),
              median(vaf_df$vaf,na.rm=TRUE)))
  p_vaf <- ggplot(vaf_df,aes(x=vaf))+
    geom_histogram(bins=60,fill="#185FA5",colour="white",linewidth=0.2)+
    scale_x_continuous(limits=c(0,1),labels=scales::percent_format())+
    labs(title=sprintf("VAF Distribution — %s",DATASET),
         subtitle=sprintf("n=%s mutations | t_alt_count/t_depth",
                          format(nrow(vaf_df),big.mark=",")),
         x="Variant Allele Frequency",y="Number of Mutations")+
    theme_minimal(base_size=13)+
    theme(plot.background=element_rect(fill="white",colour=NA))
  print(p_vaf)
  ggsave(file.path(OUT_F,"fig06_vaf.png"),p_vaf,width=9,height=7,dpi=300)
  cat("    Used: t_alt_count / t_depth\n")
} else {
  cat("    VAF columns not available\n")
}
cat(sprintf("    Saved: %s/fig06_vaf.png\n",OUT_F))

