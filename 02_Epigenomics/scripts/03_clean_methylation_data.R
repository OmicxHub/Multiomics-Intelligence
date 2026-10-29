#!/usr/bin/env Rscript
# ============================================================
# 03_clean_methylation_data.R
#
# Quality control (QC), probe filtering, probe annotation summaries and
# exploratory analysis of the TCGA-BRCA 450K beta-value matrix.
# (Month 2, Week 2: "DNA methylation preprocessing, probe filtering, QC
#  and annotation".)
#
# WHAT THIS SCRIPT DOES, IN ORDER
#   1. Loads the assembled object from script 02b
#   2. Sample QC       - missingness, methylation level, outlier detection,
#                        batch (plate) effects
#   3. Probe filtering - removes non-CpG, unmapped, sex-chromosome and
#                        quality-masked probes, then probes with too many
#                        missing values (a step-by-step "waterfall")
#   4. Probe annotation summaries - where the kept CpGs sit in the genome
#                        (CpG islands, promoters, chromosomes) and their
#                        methylation level in each context
#   5. Replicate concordance - the 6 patients with 2-3 tumour samples give a
#                        free test of technical reproducibility
#   6. PCA + covariate association - which biological or technical factors
#                        drive the main axes of variation
#   7. Unsupervised clustering + heatmaps, compared with PAM50 subtypes and
#                        the published TCGA methylation clusters
#   8. Global methylation by PAM50 subtype
#   9. Saves the cleaned object for the next scripts
#
# INPUT : data/processed/methylation_beta_matrix_raw.rds   (script 02b)
# OUTPUT: data/processed/methylation_beta_matrix_cleaned.rds
#           SummarizedExperiment - QC-passed samples (replicates still
#           included; collapsed to one per patient in script 04), filtered
#           probes, per-probe stats in rowData, per-sample QC in colData
#         data/processed/methylation_probe_stats.csv
#         outputs/figures/03_qc/  - 24 numbered figures (PNG, 300 dpi)
#         outputs/tables/03_qc/   - summary tables (CSV)
#         outputs/reports/03_qc_report.txt - everything printed on screen
#
# MEMORY: the raw matrix is ~2.5 GB. All statistics are computed with
# matrixStats directly on that matrix (no copies). Only the very last step
# makes one filtered copy (~2 GB) to save, so peak memory is ~5 GB.
# Close web browsers before running.
#
# Run from the Ubuntu terminal with:
#   Rscript scripts/03_clean_methylation_data.R
# Test on the pilot object (outputs get a _pilot40 suffix):
#   Rscript scripts/03_clean_methylation_data.R \
#       --input=data/processed/methylation_beta_matrix_raw_pilot40.rds --tag=pilot40
# ============================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)  # container: beta matrix + probe/sample annotation
  library(matrixStats)           # fast row/column statistics without copying the matrix
  library(data.table)            # fast table building and CSV writing
  library(ggplot2)               # all standard plots
  library(pheatmap)              # clustered heatmaps
  library(RColorBrewer)          # colour palettes
})

# ============================================================
# 0. SETTINGS, PATHS AND HELPER FUNCTIONS
# ============================================================
setwd("/mnt/e/R_Studio/Projects/MOI/Month_2")

# ---- Command-line flags (--input=..., --tag=...) ----
args    <- commandArgs(trailingOnly = TRUE)
getflag <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}
INPUT <- getflag("input", "data/processed/methylation_beta_matrix_raw.rds")
TAG   <- getflag("tag", "")              # e.g. "pilot40" keeps test outputs separate
sfx   <- if (nzchar(TAG)) paste0("_", TAG) else ""

# ---- Output locations ----
FIG_DIR     <- paste0("outputs/figures/03_qc", sfx)
TAB_DIR     <- paste0("outputs/tables/03_qc", sfx)
REPORT      <- paste0("outputs/reports/03_qc_report", sfx, ".txt")
OUT_RDS     <- paste0("data/processed/methylation_beta_matrix_cleaned", sfx, ".rds")
PROBE_STATS <- paste0("data/processed/methylation_probe_stats", sfx, ".csv")
for (d in c(FIG_DIR, TAB_DIR, dirname(REPORT), dirname(OUT_RDS)))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- QC thresholds (change here, everything below follows) ----
SAMPLE_NA_MAX  <- 10     # % missing among analysis probes above which a sample FAILS
                         #   (same idea as the ChAMP default: >10% failed probes)
PROBE_NA_MAX   <- 10     # % of samples missing above which a probe is REMOVED
                         #   (keeps probes measured in >= 90% of samples; the few
                         #    remaining gaps are imputed later for modelling)
OUTLIER_Z      <- 3      # robust z-score beyond which a sample is FLAGGED (kept,
                         #   but listed for inspection)
TOP_VAR_N      <- 5000   # most variable probes used for PCA / clustering
HEATMAP_N      <- 1000   # most variable probes drawn in the probe heatmap
DENSITY_PROBES <- 20000  # random probes used for density curves (speed)
N_UNRELATED    <- 300    # random unrelated sample pairs for the replicate comparison
SEED           <- 42     # same seed as Month 1, for reproducibility
set.seed(SEED)

# ---- Send everything printed to the screen AND to the report file ----
report_con <- file(REPORT, open = "wt")
sink(report_con, split = TRUE)

# ---- Timing / section headers ----
run_start <- Sys.time()
fmt_time <- function(secs) {
  if (secs < 60) sprintf("%.0f s", secs) else sprintf("%.1f min", secs / 60)
}
step <- function(title) {
  el <- as.numeric(difftime(Sys.time(), run_start, units = "secs"))
  cat(sprintf("\n[%s | %s] ================ %s ================\n",
              format(Sys.time(), "%H:%M:%S"), fmt_time(el), title))
}

# ---- Keep a list of every file written (saved as an index at the end) ----
outputs_written <- list()
record <- function(path, what) {
  outputs_written[[length(outputs_written) + 1]] <<-
    data.frame(type = what, file = path, stringsAsFactors = FALSE)
}

# ---- Consistent plot style: large readable text, legend at the bottom ----
theme_qc <- theme_minimal(base_size = 12) +
  theme(plot.title       = element_text(face = "bold", size = 14),
        plot.subtitle    = element_text(colour = "grey30", size = 10.5),
        plot.caption     = element_text(colour = "grey45", size = 9, hjust = 0),
        legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        strip.text       = element_text(face = "bold"),
        plot.margin      = margin(10, 15, 10, 10))

# Saves a ggplot at 300 dpi on a white background and records it.
save_plot <- function(p, name, width, height) {
  path <- file.path(FIG_DIR, name)
  ggsave(path, p, width = width, height = height, dpi = 300, bg = "white")
  record(path, "figure")
  cat("  saved figure:", path, "\n")
}

# Writes a CSV table and records it.
save_table <- function(df, name, dir = TAB_DIR) {
  path <- file.path(dir, name)
  fwrite(df, path)
  record(path, "table")
  cat("  saved table :", path, "\n")
}

# Prints a data.frame without row names (cleaner console output).
show <- function(df) print(as.data.frame(df), row.names = FALSE)

# Robust z-score: like a z-score but uses median/MAD, so outliers do not
# distort the reference they are measured against.
robust_z <- function(x) {
  s <- mad(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - median(x, na.rm = TRUE)) / s
}

# Adjusted Rand index: agreement between two groupings (1 = identical,
# 0 = no better than chance).
ari <- function(a, b) {
  ok <- !is.na(a) & !is.na(b); tab <- table(a[ok], b[ok]); n <- sum(tab)
  if (n < 2) return(NA_real_)
  s_ij <- sum(choose(tab, 2)); s_a <- sum(choose(rowSums(tab), 2))
  s_b  <- sum(choose(colSums(tab), 2)); expct <- s_a * s_b / choose(n, 2)
  (s_ij - expct) / ((s_a + s_b) / 2 - expct)
}

fmt_n   <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_p   <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-4, formatC(p, format = "e", digits = 1),
                                                     formatC(p, format = "f", digits = 4)))

# ---- Colour palettes ----
PAM50_LEVELS <- c("LumA", "LumB", "Her2", "Basal", "Normal", "Unknown")
PAL_PAM50 <- c(LumA = "#2C7BB6", LumB = "#74ADD1", Her2 = "#F46D43",
               Basal = "#D73027", Normal = "#1A9850", Unknown = "grey70")
PAL_STATUS <- c(Pass = "#185FA5", Flagged = "#F2A900", Fail = "#EF4444")  # Month 1 colours
PAL_CONTEXT <- c(Island = "#1B9E77", Shore = "#D95F02", Shelf = "#7570B3", OpenSea = "#E7298A")

cat("Script 03 - methylation QC, filtering and exploration\n")
cat("Input  :", INPUT, "\n")
cat("Tag    :", ifelse(nzchar(TAG), TAG, "(none - real run)"), "\n")
cat("Started:", format(run_start), "\n")

# ============================================================
# 1. LOAD THE ASSEMBLED OBJECT
# ============================================================
step("1. Load data")
if (!file.exists(INPUT)) stop("Input not found: ", INPUT, "\nRun scripts/02b_assemble_beta_matrix.R first.")

se   <- readRDS(INPUT)                              # SummarizedExperiment from script 02b
beta <- assay(se, "beta")                           # probes x samples matrix (no copy made)
rd   <- as.data.frame(rowData(se))                  # probe annotation
cd   <- as.data.frame(colData(se))                  # sample annotation
rm(se); invisible(gc())                             # drop the container; the matrix stays

# PAM50: the TCGA table uses the text "NA" for a few patients; treat it,
# and patients missing from the table, as "Unknown".
cd$pam50[cd$pam50 %in% c("NA", "")] <- NA
cd$pam50_label <- factor(ifelse(is.na(cd$pam50), "Unknown", cd$pam50), levels = PAM50_LEVELS)

# Published TCGA DNA-methylation cluster (from TCGAquery_subtype in 02b).
mc <- cd$meth_cluster_tcga
cd$meth_cluster_label <- ifelse(is.na(mc) | mc %in% c("NA", ""), "Unknown", paste0("TCGA-", mc))

n_p <- nrow(beta); n_s <- ncol(beta)
overall_na <- sum(colCounts(beta, value = NA)) / (as.numeric(n_p) * n_s) * 100
rng <- c(min(colMins(beta, na.rm = TRUE)), max(colMaxs(beta, na.rm = TRUE)))  # no matrix copy

overview <- data.frame(
  metric = c("Probes", "Samples", "Patients", "Replicate patients (>1 sample)",
             "Overall missing values (%)", "Minimum beta", "Maximum beta",
             paste0("PAM50: ", PAM50_LEVELS)),
  value  = c(fmt_n(n_p), n_s, length(unique(cd$patient_id)),
             length(unique(cd$patient_id[cd$is_replicate_patient])),
             sprintf("%.2f", overall_na), sprintf("%.4f", rng[1]), sprintf("%.4f", rng[2]),
             as.integer(table(cd$pam50_label)[PAM50_LEVELS]))) # stored as text: mixed number formats
cat("\n--- Dataset overview ---\n"); show(overview)
save_table(overview, "03_00_dataset_overview.csv")

# ============================================================
# 2. DEFINE THE "ANALYSIS PROBE" SET (design-level filters)
#    These filters depend only on the probe itself, not on the data, so
#    they are applied first and used for sample QC as well.
# ============================================================
step("2. Design-level probe filters")
is_cg     <- rd$probe_type == "cg"                         # real CpG sites (drop ch/rs/ct)
is_mapped <- !is.na(rd$chr)                                # has an hg38 position
is_auto   <- rd$chr %in% paste0("chr", 1:22)               # autosomes only: X/Y methylation
                                                           #   reflects sex / X-inactivation
not_mask  <- !rd$mask_general                              # not SNP-affected / cross-reactive
                                                           #   (Zhou et al. 2017 general mask)
design_ok  <- is_cg & is_mapped & is_auto & not_mask
design_idx <- which(design_ok)
cat("Probes passing design-level filters:", fmt_n(length(design_idx)), "of", fmt_n(n_p), "\n")

# ============================================================
# 3. SAMPLE QC
# ============================================================
step("3. Sample QC")

# ---- 3a. Per-sample metrics ----
# Some analysis probes are missing in (almost) EVERY sample: the GDC SeSAMe
# pipeline masks them systematically. That is a property of the PROBE, not
# of the sample, so those probes are left out of the sample metrics -
# otherwise every sample starts with a ~7-8% "missing" baseline and good
# samples look like failures. (They are removed later by the probe filter.)
probe_na_all <- rowCounts(beta, rows = design_idx, value = NA) / n_s * 100
qc_idx <- design_idx[probe_na_all < 95]                   # probes that normally work
cat("Analysis probes masked in >= 95% of samples (systematic, excluded from sample QC):",
    fmt_n(length(design_idx) - length(qc_idx)), "\n")
cat("Probes used for sample QC:", fmt_n(length(qc_idx)), "\n")

na_pct_total <- colCounts(beta, rows = design_idx, value = NA) / length(design_idx) * 100  # for reference
na_pct       <- colCounts(beta, rows = qc_idx,     value = NA) / length(qc_idx) * 100     # sample-specific failures
mean_beta    <- colMeans2(beta, rows = qc_idx, na.rm = TRUE)
median_beta  <- colMedians(beta, rows = qc_idx, na.rm = TRUE)

# Share of CpGs that are unmethylated (<0.2) or methylated (>0.8). Done in
# blocks of 40 samples so the temporary TRUE/FALSE matrices stay small.
pct_hypo <- pct_hyper <- numeric(n_s)
for (cc in split(seq_len(n_s), ceiling(seq_len(n_s) / 40))) {
  blk   <- beta[qc_idx, cc, drop = FALSE]
  n_obs <- colSums(!is.na(blk))
  pct_hypo[cc]  <- colSums(blk < 0.2, na.rm = TRUE) / n_obs * 100
  pct_hyper[cc] <- colSums(blk > 0.8, na.rm = TRUE) / n_obs * 100
}
rm(blk); invisible(gc())

# Correlation of each sample with the "median sample" profile, on a random
# subset of probes. A low value means the sample looks unlike the others.
dens_idx <- sort(sample(qc_idx, min(DENSITY_PROBES, length(qc_idx))))
dens_mat <- beta[dens_idx, , drop = FALSE]                 # small: 20,000 x samples
med_prof <- rowMedians(dens_mat, na.rm = TRUE)
cor_to_median <- apply(dens_mat, 2, function(v) cor(v, med_prof, use = "complete.obs"))

# ---- 3b. Pass / Flag / Fail rules ----
z_na   <- robust_z(na_pct)
z_mean <- robust_z(mean_beta)
z_cor  <- robust_z(cor_to_median)
fail   <- na_pct > SAMPLE_NA_MAX                           # hard failure: too many failed probes
flag   <- !fail & (abs(z_na) > OUTLIER_Z | abs(z_mean) > OUTLIER_Z | z_cor < -OUTLIER_Z)
status <- ifelse(fail, "Fail", ifelse(flag, "Flagged", "Pass"))
reason <- mapply(function(f, zn, zm, zc) {
  r <- c(if (f) paste0("missing > ", SAMPLE_NA_MAX, "%"),
         if (!f && abs(zn) > OUTLIER_Z) "unusual missingness",
         if (!f && abs(zm) > OUTLIER_Z) "unusual mean methylation",
         if (!f && zc < -OUTLIER_Z)     "low correlation to other samples")
  if (length(r)) paste(r, collapse = "; ") else ""
}, fail, z_na, z_mean, z_cor)

sample_qc <- data.frame(
  barcode = cd$barcode, patient_id = cd$patient_id, pam50 = cd$pam50_label,
  plate = cd$plate, tss = cd$tss, vial = cd$vial,
  replicate_patient = cd$is_replicate_patient,
  pct_missing = round(na_pct, 3),                  # sample-specific failures (used for QC)
  pct_missing_incl_systematic = round(na_pct_total, 3),
  mean_beta = round(mean_beta, 4),
  median_beta = round(median_beta, 4), pct_unmethylated = round(pct_hypo, 2),
  pct_methylated = round(pct_hyper, 2), cor_to_median = round(cor_to_median, 4),
  z_missing = round(z_na, 2), z_mean_beta = round(z_mean, 2), z_cor = round(z_cor, 2),
  qc_status = status, qc_reason = reason, stringsAsFactors = FALSE)
save_table(sample_qc, "03_01_sample_qc_metrics.csv")

keep_sample <- status != "Fail"                            # Flagged samples are KEPT
ks_idx      <- which(keep_sample)

cat(sprintf("\nSample-specific missing values : median %.2f%%, range %.2f%% - %.2f%%\n",
            median(na_pct), min(na_pct), max(na_pct)))
cat(sprintf("(including systematic masks    : median %.2f%%, range %.2f%% - %.2f%%)\n",
            median(na_pct_total), min(na_pct_total), max(na_pct_total)))
cat(sprintf("Mean beta per sample      : median %.3f, range %.3f - %.3f\n",
            median(mean_beta), min(mean_beta), max(mean_beta)))
cat(sprintf("Correlation to median     : median %.3f, minimum %.3f\n",
            median(cor_to_median), min(cor_to_median)))
cat("\nQC status:\n"); print(table(factor(status, levels = c("Pass", "Flagged", "Fail"))))
if (any(status != "Pass")) {
  cat("\nSamples that failed or were flagged:\n")
  show(sample_qc[status != "Pass", c("barcode", "pam50", "plate", "pct_missing",
                                     "mean_beta", "cor_to_median", "qc_status", "qc_reason")])
  save_table(sample_qc[status != "Pass", ], "03_02_samples_failed_or_flagged.csv")
}

# ---- 3c. Batch view: summary per processing plate ----
plate_summary <- as.data.frame(data.table(sample_qc)[, .(
  n_samples = .N,
  median_pct_missing = round(median(pct_missing), 2),
  median_mean_beta   = round(median(mean_beta), 4),
  n_flag_or_fail     = sum(qc_status != "Pass")), by = plate][order(-n_samples)])
cat("\nProcessing plates (batches):", nrow(plate_summary), "\n")
show(head(plate_summary, 10)); if (nrow(plate_summary) > 10) cat("  ... (full table saved)\n")
save_table(plate_summary, "03_03_plate_summary.csv")

# ---- 3d. Figures ----
# Fig 1: missingness of every sample, sorted - shows the threshold and any failures.
qdf <- sample_qc[order(sample_qc$pct_missing), ]
qdf$rank <- seq_len(nrow(qdf))
p <- ggplot(qdf, aes(rank, pct_missing, colour = qc_status)) +
  geom_hline(yintercept = SAMPLE_NA_MAX, linetype = "dashed", colour = PAL_STATUS[["Fail"]]) +
  annotate("text", x = 1, y = SAMPLE_NA_MAX, label = paste0("fail threshold ", SAMPLE_NA_MAX, "%"),
           hjust = 0, vjust = -0.6, size = 3.5, colour = PAL_STATUS[["Fail"]]) +
  geom_point(size = 1.6, alpha = 0.85) +
  scale_colour_manual(values = PAL_STATUS, name = "QC status") +
  scale_y_continuous(limits = c(0, max(SAMPLE_NA_MAX * 1.25, max(qdf$pct_missing) * 1.1))) +
  labs(title = "Missing values per sample (sample-specific failures)",
       subtitle = sprintf("%s autosomal, unmasked CpG probes that work in most samples; samples sorted low to high",
                          fmt_n(length(qc_idx))),
       x = "Samples (ranked)", y = "Missing values (%)") + theme_qc
save_plot(p, "03_01_sample_missingness.png", 9, 5.5)

# Fig 2: missingness by processing plate - a plate with high values is a batch problem.
plate_order <- plate_summary$plate[order(plate_summary$median_pct_missing)]
p <- ggplot(sample_qc, aes(factor(plate, levels = plate_order), pct_missing)) +
  geom_boxplot(fill = "#DCE9F5", colour = "#185FA5", outlier.size = 0.8) +
  labs(title = "Missing values by processing plate (batch)",
       subtitle = "Plates ordered by median; a plate sitting clearly higher would indicate a batch problem",
       x = "Processing plate", y = "Missing values (%)") +
  theme_qc + theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8))
save_plot(p, "03_02_missingness_by_plate.png", 12, 5.5)

# Fig 3: beta-value distribution of every sample, one panel per PAM50 subtype.
# Healthy 450K data is bimodal: peaks near 0 (unmethylated) and near 1.
dens_df <- rbindlist(lapply(seq_len(n_s), function(j) {
  v <- dens_mat[, j]; v <- v[!is.na(v)]
  d <- density(v, from = 0, to = 1, n = 128)
  data.table(barcode = colnames(dens_mat)[j], x = d$x, y = d$y)
}))
dens_df[, pam50 := cd$pam50_label[match(barcode, cd$barcode)]]
dens_med <- dens_df[, .(y = median(y)), by = .(pam50, x)]
p <- ggplot(dens_df, aes(x, y, group = barcode, colour = pam50)) +
  geom_line(alpha = 0.18, linewidth = 0.3) +
  geom_line(data = dens_med, aes(x, y, group = pam50), colour = "black", linewidth = 0.8,
            inherit.aes = FALSE) +
  facet_wrap(~ pam50, ncol = 3) +
  scale_colour_manual(values = PAL_PAM50, guide = "none") +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(title = "Beta-value distribution of each sample, by PAM50 subtype",
       subtitle = sprintf("One line per sample (%s random analysis probes); black = median curve",
                          fmt_n(length(dens_idx))),
       x = "Beta value (0 = unmethylated, 1 = fully methylated)", y = "Density") + theme_qc
save_plot(p, "03_03_beta_density_by_pam50.png", 11, 7)
rm(dens_df, dens_med); invisible(gc())

# Fig 4: missingness vs overall methylation - outliers stand apart from the cloud.
p <- ggplot(sample_qc, aes(pct_missing, mean_beta, colour = qc_status)) +
  geom_point(size = 1.8, alpha = 0.8) +
  scale_colour_manual(values = PAL_STATUS, name = "QC status") +
  labs(title = "Sample QC overview",
       subtitle = "Each point is a sample; flagged samples are listed in 03_02_samples_failed_or_flagged.csv",
       x = "Missing values (%)", y = "Mean beta (analysis probes)") + theme_qc
save_plot(p, "03_04_sample_qc_scatter.png", 8, 6)

# ============================================================
# 4. PROBE FILTERING (step-by-step waterfall)
# ============================================================
step("4. Probe filtering")

# % of QC-passed samples in which each analysis probe is missing.
probe_na_design <- rowCounts(beta, rows = design_idx, cols = ks_idx, value = NA) / length(ks_idx) * 100
na_ok      <- probe_na_design <= PROBE_NA_MAX
keep_idx   <- design_idx[na_ok]                           # FINAL retained probes (row numbers)

cat(sprintf("Analysis probes never missing           : %s\n", fmt_n(sum(probe_na_design == 0))))
cat(sprintf("Analysis probes missing in ALL samples  : %s\n", fmt_n(sum(probe_na_design == 100))))
cat(sprintf("Analysis probes missing in > %d%% samples: %s (removed)\n",
            PROBE_NA_MAX, fmt_n(sum(!na_ok))))

waterfall <- data.frame(
  step = c("0. All probes on the array",
           "1. CpG probes only (drop ch / rs / control)",
           "2. Mapped to hg38",
           "3. Autosomes only (drop chrX / chrY / chrM)",
           "4. Not quality-masked (SNP / cross-reactive)",
           paste0("5. Missing in <= ", PROBE_NA_MAX, "% of samples")),
  n_remaining = c(n_p, sum(is_cg), sum(is_cg & is_mapped), sum(is_cg & is_mapped & is_auto),
                  length(design_idx), length(keep_idx)))
waterfall$n_removed     <- c(0, -diff(waterfall$n_remaining))
waterfall$pct_remaining <- round(100 * waterfall$n_remaining / n_p, 1)
cat("\n--- Probe filtering waterfall ---\n"); show(waterfall)
save_table(waterfall, "03_04_probe_filtering_waterfall.csv")

# Fig 5: waterfall - how many probes each filter removes.
wf <- waterfall
wf$step  <- factor(wf$step, levels = rev(wf$step))
wf$label <- ifelse(wf$n_removed > 0,
                   sprintf("%s  (-%s)", fmt_n(wf$n_remaining), fmt_n(wf$n_removed)),
                   fmt_n(wf$n_remaining))
p <- ggplot(wf, aes(n_remaining, step)) +
  geom_col(fill = "#185FA5", width = 0.65) +
  geom_text(aes(label = label), hjust = -0.05, size = 3.6) +
  scale_x_continuous(labels = fmt_n, expand = expansion(mult = c(0, 0.32))) +
  labs(title = "Probe filtering: probes remaining after each step",
       subtitle = "Number in brackets = probes removed by that step",
       x = "Probes remaining", y = NULL) + theme_qc
save_plot(p, "03_05_probe_filtering_waterfall.png", 11, 5)

# Fig 6: distribution of probe missingness (log scale - most probes have none).
p <- ggplot(data.frame(pct = probe_na_design), aes(pct)) +
  geom_histogram(bins = 50, fill = "#0F6E56", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = PROBE_NA_MAX, linetype = "dashed", colour = PAL_STATUS[["Fail"]]) +
  annotate("text", x = PROBE_NA_MAX, y = Inf, label = paste0(" removal threshold ", PROBE_NA_MAX, "%"),
           hjust = 0, vjust = 1.5, size = 3.5, colour = PAL_STATUS[["Fail"]]) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 10, 100, 1e3, 1e4, 1e5), labels = fmt_n) +
  labs(title = "How often each analysis probe is missing",
       subtitle = "Across QC-passed samples; y-axis is logarithmic",
       x = "Samples in which the probe is missing (%)", y = "Number of probes") + theme_qc
save_plot(p, "03_06_probe_missingness_histogram.png", 9, 5.5)

# ============================================================
# 5. PROBE STATISTICS AND ANNOTATION SUMMARIES
# ============================================================
step("5. Probe annotation summaries")

# Per-probe mean, SD and missingness over retained probes and passing samples.
p_mean <- rowMeans2(beta, rows = keep_idx, cols = ks_idx, na.rm = TRUE)
p_sd   <- rowSds(beta,   rows = keep_idx, cols = ks_idx, na.rm = TRUE)
p_na   <- probe_na_design[na_ok]

# Per-probe median. NOTE: matrixStats 1.5.0 rowMedians() silently IGNORES
# its `rows =` argument on this data (verified: it returned all 486,427
# probes instead of the subset). So medians are computed on explicit blocks
# of 20,000 probes instead - each block is a small temporary copy.
block_row_medians <- function(m, rows, cols, block = 20000) {
  out <- numeric(length(rows))
  for (b in split(seq_along(rows), ceiling(seq_along(rows) / block)))
    out[b] <- rowMedians(m[rows[b], cols, drop = FALSE], na.rm = TRUE)
  out
}
p_med <- block_row_medians(beta, keep_idx, ks_idx)
stopifnot(length(p_med) == length(keep_idx))   # guard against the bug above

kept_rd <- rd[keep_idx, c("probe_id", "chr", "pos_hg38", "gene", "dist_to_tss",
                          "cgi_context", "region", "design_type")]
probe_stats <- cbind(kept_rd, mean_beta = round(p_mean, 4), median_beta = round(p_med, 4),
                     sd_beta = round(p_sd, 4), pct_missing = round(p_na, 2))
fwrite(probe_stats, PROBE_STATS); record(PROBE_STATS, "processed data")
cat("  saved probe statistics:", PROBE_STATS, "(", fmt_n(nrow(probe_stats)), "probes )\n")

# ---- CpG-island context and genomic region: all CpG probes vs retained ----
ctx_levels <- c("Island", "Shore", "Shelf", "OpenSea")
reg_levels <- c("Promoter (TSS +/-1500)", "Gene body", "Intergenic")
compare_sets <- function(col, levels) {
  before <- table(factor(rd[[col]][is_cg & is_mapped], levels = levels))
  after  <- table(factor(rd[[col]][keep_idx], levels = levels))
  data.frame(category = levels,
             all_cpg_probes = as.integer(before), retained = as.integer(after),
             pct_of_all = round(100 * as.integer(before) / sum(before), 1),
             pct_of_retained = round(100 * as.integer(after) / sum(after), 1),
             pct_category_kept = round(100 * as.integer(after) / pmax(as.integer(before), 1), 1))
}
ctx_tab <- compare_sets("cgi_context", ctx_levels)
reg_tab <- compare_sets("region", reg_levels)
cat("\n--- CpG-island context (all mapped CpG probes vs retained) ---\n"); show(ctx_tab)
cat("\n--- Genomic region ---\n"); show(reg_tab)
save_table(rbind(cbind(annotation = "CpG-island context", ctx_tab),
                 cbind(annotation = "Genomic region", reg_tab)), "03_05_probe_annotation_summary.csv")

chr_tab <- as.data.frame(table(chr = factor(kept_rd$chr, levels = paste0("chr", 1:22))))
names(chr_tab)[2] <- "retained_probes"
save_table(chr_tab, "03_06_probes_per_chromosome.csv")

# ---- Methylation level by genomic context (the biological sanity check) ----
by_ctx <- data.table(context = kept_rd$cgi_context, region = kept_rd$region,
                     design = kept_rd$design_type, mean_beta = p_mean)
# Summarises probe methylation within each level of one annotation column.
ctx_summary <- function(dt, col, label) {
  s <- dt[, .(n_probes = .N,
              median_of_probe_means   = round(median(mean_beta), 3),
              pct_mostly_unmethylated = round(100 * mean(mean_beta < 0.2), 1),
              pct_mostly_methylated   = round(100 * mean(mean_beta > 0.8), 1)), by = col]
  setnames(s, col, "level")
  s[, grouping := label]
  setcolorder(s, c("grouping", "level"))
  s[order(level)]
}
beta_ctx <- rbind(ctx_summary(by_ctx, "context", "CpG-island context"),
                  ctx_summary(by_ctx, "region",  "Genomic region"),
                  ctx_summary(by_ctx, "design",  "Infinium design"))
cat("\n--- Methylation level by context (retained probes) ---\n"); show(beta_ctx)
save_table(beta_ctx, "03_07_methylation_by_genomic_context.csv")

# Fig 7 / 8: composition before vs after filtering (percent labels on bars).
comp_plot <- function(tab, title, levels, pal = NULL) {
  long <- rbind(data.frame(set = "All mapped CpG probes", category = tab$category,
                           n = tab$all_cpg_probes, pct = tab$pct_of_all),
                data.frame(set = "Retained after QC", category = tab$category,
                           n = tab$retained, pct = tab$pct_of_retained))
  long$category <- factor(long$category, levels = levels)
  ggplot(long, aes(category, pct, fill = set)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_text(aes(label = paste0(pct, "%")), position = position_dodge(width = 0.8),
              vjust = -0.4, size = 3.3) +
    scale_fill_manual(values = c("All mapped CpG probes" = "grey70", "Retained after QC" = "#185FA5"),
                      name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = title, x = NULL, y = "Share of probes (%)") + theme_qc
}
p <- comp_plot(ctx_tab, "CpG-island context of probes, before vs after filtering", ctx_levels) +
  labs(subtitle = "Shore = up to 2 kb from an island, Shelf = 2-4 kb, OpenSea = further away")
save_plot(p, "03_07_probe_cgi_context.png", 9, 5.5)
p <- comp_plot(reg_tab, "Genomic region of probes, before vs after filtering", reg_levels) +
  labs(subtitle = "Promoter = within 1,500 bp of a transcription start site")
save_plot(p, "03_08_probe_genomic_region.png", 9, 5.5)

# Fig 9: retained probes per chromosome.
p <- ggplot(chr_tab, aes(chr, retained_probes)) +
  geom_col(fill = "#0F6E56", width = 0.75) +
  scale_y_continuous(labels = fmt_n, expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Retained probes per chromosome",
       subtitle = "Autosomes only; larger and gene-rich chromosomes carry more probes",
       x = NULL, y = "Retained probes") +
  theme_qc + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p, "03_09_probes_per_chromosome.png", 10, 5)

# Fig 10 / 11: probe mean methylation by context (violin = distribution, box = quartiles).
violin_plot <- function(dt, xcol, levels, title, subtitle, pal) {
  d <- data.frame(x = factor(dt[[xcol]], levels = levels), y = dt$mean_beta)
  med <- aggregate(y ~ x, d, median)
  ggplot(d, aes(x, y, fill = x)) +
    geom_violin(scale = "width", colour = NA, alpha = 0.75) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
    geom_text(data = med, aes(x, 1.06, label = sprintf("median %.2f", y)),
              inherit.aes = FALSE, size = 3.4) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Mean beta across samples") + theme_qc
}
p <- violin_plot(by_ctx, "context", ctx_levels, "Methylation level by CpG-island context",
                 "Expected biology: CpG islands mostly unmethylated, open sea mostly methylated",
                 PAL_CONTEXT)
save_plot(p, "03_10_methylation_by_cgi_context.png", 9, 5.5)
p <- violin_plot(by_ctx, "region", reg_levels, "Methylation level by genomic region",
                 "Promoters are typically unmethylated (active genes); gene bodies methylated",
                 c("Promoter (TSS +/-1500)" = "#1B9E77", "Gene body" = "#7570B3", "Intergenic" = "#D95F02"))
save_plot(p, "03_11_methylation_by_region.png", 9, 5.5)

# Fig 12: Infinium I vs II design - known technical difference in beta distributions.
p <- ggplot(data.frame(design = paste("Infinium", kept_rd$design_type), median_beta = p_med),
            aes(median_beta, colour = design, fill = design)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  scale_colour_manual(values = c("Infinium I" = "#D95F02", "Infinium II" = "#185FA5"), name = NULL) +
  scale_fill_manual(values = c("Infinium I" = "#D95F02", "Infinium II" = "#185FA5"), name = NULL) +
  labs(title = "Beta-value distribution by probe design (Infinium I vs II)",
       subtitle = "Mostly reflects WHERE each design sits (Type I probes are concentrated in CpG islands),\nplus the known narrower dynamic range of Type II probes",
       x = "Median beta across samples", y = "Density") + theme_qc
save_plot(p, "03_12_beta_by_probe_design.png", 9, 5.5)

# Fig 13: probe variability - the most variable probes carry the biological signal.
sd_cut <- sort(p_sd, decreasing = TRUE)[min(TOP_VAR_N, length(p_sd))]
p <- ggplot(data.frame(sd = p_sd), aes(sd)) +
  geom_histogram(bins = 80, fill = "#7570B3", colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = sd_cut, linetype = "dashed", colour = PAL_STATUS[["Fail"]]) +
  annotate("text", x = sd_cut, y = Inf, label = paste0(" top ", fmt_n(TOP_VAR_N), " most variable"),
           hjust = 0, vjust = 1.5, size = 3.5, colour = PAL_STATUS[["Fail"]]) +
  scale_y_continuous(labels = fmt_n) +
  labs(title = "Variability of retained probes across samples",
       subtitle = "Most CpGs barely vary between tumours; a minority vary strongly",
       x = "Standard deviation of beta", y = "Number of probes") + theme_qc
save_plot(p, "03_13_probe_variability.png", 9, 5.5)

# ============================================================
# 6. REPLICATE CONCORDANCE (technical reproducibility)
# ============================================================
step("6. Replicate concordance")
pair_stats <- function(i, j) {
  a <- beta[keep_idx, i]; b <- beta[keep_idx, j]; ok <- !is.na(a) & !is.na(b)
  c(r = cor(a[ok], b[ok]), mean_abs_diff = mean(abs(a[ok] - b[ok])),
    pct_diff_gt_0.2 = 100 * mean(abs(a[ok] - b[ok]) > 0.2), n_probes = sum(ok))
}
rep_pat <- unique(cd$patient_id[cd$is_replicate_patient & keep_sample])
rep_pairs <- do.call(rbind, lapply(rep_pat, function(pt) {
  idx <- which(cd$patient_id == pt & keep_sample)
  if (length(idx) < 2) return(NULL)
  cmb <- combn(idx, 2)
  do.call(rbind, lapply(seq_len(ncol(cmb)), function(k) {
    s <- pair_stats(cmb[1, k], cmb[2, k])
    data.frame(group = "Same patient (replicates)", patient_id = pt,
               sample_1 = cd$barcode[cmb[1, k]], sample_2 = cd$barcode[cmb[2, k]],
               r = s[["r"]], mean_abs_diff = s[["mean_abs_diff"]],
               pct_diff_gt_0.2 = s[["pct_diff_gt_0.2"]], n_probes = s[["n_probes"]])
  }))
}))

# Random pairs of DIFFERENT patients as the comparison group.
unrel <- list(); tries <- 0
while (length(unrel) < N_UNRELATED && tries < N_UNRELATED * 20 && length(ks_idx) > 2) {
  tries <- tries + 1
  pr <- sample(ks_idx, 2)
  if (cd$patient_id[pr[1]] == cd$patient_id[pr[2]]) next
  s <- pair_stats(pr[1], pr[2])
  unrel[[length(unrel) + 1]] <- data.frame(
    group = "Different patients", patient_id = NA, sample_1 = cd$barcode[pr[1]],
    sample_2 = cd$barcode[pr[2]], r = s[["r"]], mean_abs_diff = s[["mean_abs_diff"]],
    pct_diff_gt_0.2 = s[["pct_diff_gt_0.2"]], n_probes = s[["n_probes"]])
}
unrel <- do.call(rbind, unrel)
conc <- rbind(rep_pairs, unrel)

if (!is.null(rep_pairs) && nrow(rep_pairs) > 0) {
  cat("Replicate pairs:", nrow(rep_pairs), "from", length(rep_pat), "patients\n")
  show(transform(rep_pairs[, c("patient_id", "r", "mean_abs_diff", "pct_diff_gt_0.2")],
                 r = round(r, 4), mean_abs_diff = round(mean_abs_diff, 4),
                 pct_diff_gt_0.2 = round(pct_diff_gt_0.2, 2)))
  cat(sprintf("\nMedian correlation - replicates: %.3f | different patients: %.3f\n",
              median(rep_pairs$r), median(unrel$r)))
  save_table(conc, "03_08_replicate_concordance.csv")

  # Fig 14: replicate pairs should correlate far more strongly than random pairs.
  conc$group <- factor(conc$group, levels = c("Same patient (replicates)", "Different patients"))
  p <- ggplot(conc, aes(group, r, colour = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.45, colour = "grey35") +
    geom_jitter(width = 0.15, height = 0, size = 1.6, alpha = 0.6) +
    scale_colour_manual(values = c("Same patient (replicates)" = "#EF4444",
                                   "Different patients" = "#185FA5"), guide = "none") +
    labs(title = "Technical reproducibility: replicate samples vs unrelated tumours",
         subtitle = sprintf("Pearson r across %s retained probes | median replicates %.3f vs unrelated %.3f",
                            fmt_n(length(keep_idx)), median(rep_pairs$r), median(unrel$r)),
         x = NULL, y = "Correlation between the two samples (r)") + theme_qc
  save_plot(p, "03_14_replicate_concordance.png", 8, 6)

  # Fig 15: one representative replicate pair (the median one) as a scatter.
  mid <- rep_pairs[order(rep_pairs$r), ][ceiling(nrow(rep_pairs) / 2), ]
  i1 <- match(mid$sample_1, cd$barcode); i2 <- match(mid$sample_2, cd$barcode)
  sc <- data.frame(a = beta[keep_idx, i1], b = beta[keep_idx, i2]); sc <- sc[complete.cases(sc), ]
  p <- ggplot(sc, aes(a, b)) +
    geom_bin_2d(bins = 120) +
    geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
    scale_fill_viridis_c(trans = "log10", name = "Probes") +
    coord_equal() +
    labs(title = paste("Replicate pair:", mid$patient_id),
         subtitle = sprintf("%s  vs  %s\nr = %.3f; red line = perfect agreement",
                            mid$sample_1, mid$sample_2, mid$r),
         x = "Beta - sample 1", y = "Beta - sample 2") +
    theme_qc + theme(legend.position = "right")
  save_plot(p, "03_15_replicate_pair_scatter.png", 8, 7.5)
  rm(sc)
} else {
  cat("No replicate pairs among QC-passed samples - section skipped.\n")
}

# ============================================================
# 7. PCA ON THE MOST VARIABLE PROBES
# ============================================================
step("7. PCA")
top_order <- order(p_sd, decreasing = TRUE)
top_idx   <- keep_idx[top_order[seq_len(min(TOP_VAR_N, length(top_order)))]]
X <- beta[top_idx, ks_idx, drop = FALSE]                  # small: 5,000 x samples
# Fill the few missing values with that PROBE's mean (row mean), so PCA can run.
na_pos <- which(is.na(X), arr.ind = TRUE)
if (nrow(na_pos)) X[na_pos] <- rowMeans(X, na.rm = TRUE)[na_pos[, 1]]
cat("Top variable probes used:", fmt_n(nrow(X)), "| missing values imputed:", fmt_n(nrow(na_pos)), "\n")

top_tab <- cbind(rd[top_idx, c("probe_id", "chr", "pos_hg38", "gene", "cgi_context", "region")],
                 sd_beta = round(p_sd[top_order[seq_len(nrow(X))]], 4),
                 mean_beta = round(p_mean[top_order[seq_len(nrow(X))]], 4))
save_table(top_tab, "03_09_top_variable_probes.csv")

pca <- prcomp(t(X), center = TRUE, scale. = FALSE)       # samples as rows; betas share one scale
var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)
n_pc    <- min(10, ncol(pca$x))
pc_df   <- data.frame(barcode = colnames(X), pca$x[, seq_len(n_pc), drop = FALSE])
cdk     <- cd[ks_idx, ]
pc_df$pam50 <- cdk$pam50_label; pc_df$meth_cluster <- cdk$meth_cluster_label
pc_df$plate <- cdk$plate

var_tab <- data.frame(PC = paste0("PC", seq_along(var_exp)), pct_variance = round(var_exp, 2),
                      cumulative_pct = round(cumsum(var_exp), 2))
cat("\nVariance explained (first 10 PCs):\n"); show(head(var_tab, 10))
save_table(head(var_tab, 20), "03_10_pca_variance_explained.csv")
save_table(pc_df, "03_11_pca_scores.csv")

# ---- Which factors explain each PC? (adjusted R-squared from a linear model) ----
# Adjusted R2 penalises factors with many levels (e.g. 34 plates), making
# biology and batch comparable. p-value = overall F-test.
covars <- list(
  "PAM50 subtype"            = ifelse(cdk$pam50_label == "Unknown", NA, as.character(cdk$pam50_label)),
  "TCGA methylation cluster" = ifelse(cdk$meth_cluster_label == "Unknown", NA, cdk$meth_cluster_label),
  "Pathology"                = ifelse(cdk$pathology %in% c("NA", ""), NA, cdk$pathology),
  "Processing plate (batch)" = cdk$plate,
  "Tissue source site"       = cdk$tss,
  "Vial (A/B)"               = cdk$vial,
  "Sample missingness (%)"   = na_pct[ks_idx])
assoc <- rbindlist(lapply(names(covars), function(cv) {
  v <- covars[[cv]]
  rbindlist(lapply(seq_len(n_pc), function(k) {
    y <- pca$x[, k]; ok <- !is.na(v)
    if (is.character(v)) {
      f <- factor(v[ok]); if (nlevels(f) < 2) return(NULL); fit <- lm(y[ok] ~ f)
    } else {
      if (length(unique(v[ok])) < 3) return(NULL); fit <- lm(y[ok] ~ v[ok])
    }
    s  <- summary(fit)
    pv <- if (is.null(s$fstatistic)) NA else
      pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)
    data.table(covariate = cv, PC = paste0("PC", k),
               adj_r2 = max(0, s$adj.r.squared), p_value = unname(pv))
  }))
}))
assoc_wide <- dcast(assoc, covariate ~ PC, value.var = "adj_r2")
cat("\nAdjusted R2 of each factor on each PC:\n")
show(assoc_wide[, c("covariate", intersect(paste0("PC", 1:n_pc), names(assoc_wide))), with = FALSE][
  , lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])
save_table(assoc, "03_12_pca_covariate_association.csv")

# Fig 16: scree plot - how much variation each PC captures.
sc_df <- head(var_tab, 20); sc_df$PC <- factor(sc_df$PC, levels = sc_df$PC)
p <- ggplot(sc_df, aes(PC, pct_variance)) +
  geom_col(fill = "#185FA5", width = 0.7) +
  geom_text(data = head(sc_df, 10), aes(label = sprintf("%.1f", pct_variance)), vjust = -0.4, size = 3.1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Variance explained by each principal component",
       subtitle = sprintf("PCA of the %s most variable probes; first 10 PCs together explain %.1f%%",
                          fmt_n(nrow(X)), sum(head(var_exp, 10))),
       x = NULL, y = "Variance explained (%)") +
  theme_qc + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p, "03_16_pca_scree.png", 10, 5)

# Fig 17 / 19: PC1 vs PC2 coloured by PAM50 subtype and by TCGA methylation cluster.
axis_lab <- function(k) sprintf("PC%d (%.1f%% of variance)", k, var_exp[k])
pca_scatter <- function(col, pal, title, legend_title) {
  ggplot(pc_df, aes(PC1, PC2, colour = .data[[col]])) +
    geom_point(size = 1.9, alpha = 0.8) +
    scale_colour_manual(values = pal, name = legend_title, drop = TRUE) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3.5, alpha = 1))) +
    labs(title = title, x = axis_lab(1), y = axis_lab(2)) + theme_qc
}
p <- pca_scatter("pam50", PAL_PAM50, "PCA of tumour methylation, coloured by PAM50 subtype", "PAM50") +
  labs(subtitle = "Separation by colour = methylation differs systematically between subtypes")
save_plot(p, "03_17_pca_by_pam50.png", 9, 7)

mc_levels <- sort(unique(pc_df$meth_cluster))
mc_levels <- c(setdiff(mc_levels, "Unknown"), intersect("Unknown", mc_levels))
pal_mc    <- setNames(c(brewer.pal(max(3, length(mc_levels)), "Set1")[seq_len(sum(mc_levels != "Unknown"))],
                        if ("Unknown" %in% mc_levels) "grey70"), mc_levels)
pc_df$meth_cluster <- factor(pc_df$meth_cluster, levels = mc_levels)
p <- pca_scatter("meth_cluster", pal_mc, "PCA coloured by published TCGA methylation cluster",
                 "TCGA cluster") +
  labs(subtitle = "Published TCGA-BRCA methylation clusters (Berger et al. 2018, Cancer Cell, TCGA Pan-Gyn study)")
save_plot(p, "03_19_pca_by_tcga_methylation_cluster.png", 9, 7)

# Fig 18: PC1 and PC2 per processing plate - checks whether batch drives the main axes.
pl_long <- rbind(data.frame(PC = axis_lab(1), plate = pc_df$plate, score = pc_df$PC1),
                 data.frame(PC = axis_lab(2), plate = pc_df$plate, score = pc_df$PC2))
pl_ord  <- names(sort(tapply(pc_df$PC1, pc_df$plate, median)))
pl_long$plate <- factor(pl_long$plate, levels = pl_ord)
p <- ggplot(pl_long, aes(plate, score)) +
  geom_boxplot(fill = "#DCE9F5", colour = "#185FA5", outlier.size = 0.6) +
  facet_wrap(~ PC, ncol = 1, scales = "free_y") +
  labs(title = "Principal component scores by processing plate (batch check)",
       subtitle = "Plates ordered by median PC1; similar boxes across plates = little batch effect",
       x = "Processing plate", y = "PC score") +
  theme_qc + theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8))
save_plot(p, "03_18_pca_by_plate.png", 12, 8)

# Fig 20: heatmap of adjusted R2 - biology vs batch at a glance.
assoc$stars <- ifelse(is.na(assoc$p_value), "", ifelse(assoc$p_value < 0.001, "***",
               ifelse(assoc$p_value < 0.01, "**", ifelse(assoc$p_value < 0.05, "*", ""))))
assoc$PC <- factor(assoc$PC, levels = paste0("PC", 1:n_pc))
assoc$covariate <- factor(assoc$covariate, levels = rev(names(covars)))
p <- ggplot(assoc, aes(PC, covariate, fill = adj_r2)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = paste0(sprintf("%.2f", adj_r2), stars),
                colour = adj_r2 > 0.6 * max(0.05, max(adj_r2))), size = 3.1) +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +  # readable on dark tiles
  scale_fill_gradient(low = "white", high = "#B2182B", limits = c(0, max(0.05, max(assoc$adj_r2))),
                      name = "Adjusted R2") +
  labs(title = "What drives each principal component?",
       subtitle = "Adjusted R2 of each factor; * p<0.05, ** p<0.01, *** p<0.001",
       x = NULL, y = NULL) +
  theme_qc + theme(panel.grid = element_blank(), legend.position = "right")
save_plot(p, "03_20_pca_covariate_association.png", 11, 5)

# ============================================================
# 8. UNSUPERVISED CLUSTERING AND HEATMAPS
# ============================================================
step("8. Clustering and heatmaps")

# Samples are clustered on 1 - correlation of their top-variable-probe
# profiles (Ward's method). The number of clusters is set to the number of
# published TCGA methylation clusters present, so the two can be compared.
cmat <- cor(X)
hc   <- hclust(as.dist(1 - cmat), method = "ward.D2")
n_tcga <- length(setdiff(unique(cdk$meth_cluster_label), "Unknown"))
K    <- if (n_tcga >= 2 && n_tcga <= 8) n_tcga else 5
K    <- min(K, ncol(X) - 1)
clus <- paste0("MC", cutree(hc, k = K))   # MC = methylation cluster (this analysis)
cat("Hierarchical clustering into K =", K, "clusters\n"); print(table(cluster = clus))

tab_pam <- table(cluster = clus, PAM50 = droplevels(cdk$pam50_label))
tab_mc  <- table(cluster = clus, TCGA_cluster = cdk$meth_cluster_label)
ari_mc  <- ari(clus, ifelse(cdk$meth_cluster_label == "Unknown", NA, cdk$meth_cluster_label))
ari_pam <- ari(clus, ifelse(cdk$pam50_label == "Unknown", NA, as.character(cdk$pam50_label)))
cat("\nClusters vs PAM50:\n"); print(tab_pam)
cat("\nClusters vs published TCGA methylation clusters:\n"); print(tab_mc)
cat(sprintf("\nAdjusted Rand index vs TCGA methylation clusters: %.3f\n", ari_mc))
cat(sprintf("Adjusted Rand index vs PAM50 subtypes          : %.3f\n", ari_pam))
cat("(1 = identical grouping, 0 = no better than chance)\n")
save_table(as.data.frame.matrix(tab_pam) |> (\(d) cbind(cluster = rownames(d), d))(),
           "03_13_clusters_vs_pam50.csv")
save_table(as.data.frame.matrix(tab_mc) |> (\(d) cbind(cluster = rownames(d), d))(),
           "03_14_clusters_vs_tcga_methylation_clusters.csv")

# Annotation bars shared by both heatmaps (sample names are hidden - too many to read).
ann_col <- data.frame(PAM50 = droplevels(cdk$pam50_label),
                      TCGA_cluster = factor(cdk$meth_cluster_label, levels = mc_levels),
                      Cluster = factor(clus), row.names = colnames(X))
# Cluster colours are deliberately different from the PAM50, TCGA-cluster and
# CpG-context palettes, so no colour means two things in the same heatmap.
PAL_CLUSTER_BASE <- c("#404040", "#8C510A", "#E6AB02", "#08306B", "#01665E",
                      "#B8E186", "#F1B6DA", "#DFC27D")
pal_cl  <- setNames(PAL_CLUSTER_BASE[seq_len(K)], paste0("MC", seq_len(K)))
ann_colors <- list(PAM50 = PAL_PAM50[levels(ann_col$PAM50)],
                   TCGA_cluster = pal_mc[levels(ann_col$TCGA_cluster)],
                   Cluster = pal_cl)

# Fig 21: sample-sample correlation heatmap.
brk <- seq(quantile(cmat, 0.01), 1, length.out = 101)
pheatmap(cmat, cluster_rows = hc, cluster_cols = hc, show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = ann_col, annotation_colors = ann_colors, treeheight_row = 0,
         color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100), breaks = brk,
         border_color = NA, fontsize = 10,
         main = sprintf("Sample-sample correlation (%s most variable probes)", fmt_n(nrow(X))),
         filename = file.path(FIG_DIR, "03_21_sample_correlation_heatmap.png"), width = 11, height = 10)
record(file.path(FIG_DIR, "03_21_sample_correlation_heatmap.png"), "figure")
cat("  saved figure:", file.path(FIG_DIR, "03_21_sample_correlation_heatmap.png"), "\n")

# Fig 22: the most variable probes x samples (beta values).
hm_idx <- seq_len(min(HEATMAP_N, nrow(X)))
ann_row <- data.frame(CpG_context = factor(rd$cgi_context[top_idx[hm_idx]], levels = ctx_levels),
                      row.names = rownames(X)[hm_idx])
pheatmap(X[hm_idx, ], cluster_cols = hc, clustering_method = "ward.D2",
         show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = ann_col, annotation_row = ann_row,
         annotation_colors = c(ann_colors, list(CpG_context = PAL_CONTEXT)),
         color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
         breaks = seq(0, 1, length.out = 101), border_color = NA, fontsize = 10,
         treeheight_row = 25, annotation_names_row = FALSE,   # legend already names it main = sprintf("Beta values of the %s most variable probes", fmt_n(length(hm_idx))),
         filename = file.path(FIG_DIR, "03_22_top_variable_probes_heatmap.png"), width = 11, height = 10)
record(file.path(FIG_DIR, "03_22_top_variable_probes_heatmap.png"), "figure")
cat("  saved figure:", file.path(FIG_DIR, "03_22_top_variable_probes_heatmap.png"), "\n")

# Fig 23: PAM50 make-up of each methylation cluster.
comp <- as.data.frame(tab_pam); names(comp) <- c("cluster", "PAM50", "n")
comp <- comp[comp$n > 0, ]
comp$share <- ave(comp$n, comp$cluster, FUN = function(z) z / sum(z))
p <- ggplot(comp, aes(cluster, n, fill = PAM50)) +
  geom_col(position = "fill", width = 0.7) +
  geom_text(aes(label = ifelse(share >= 0.08, n, "")), position = position_fill(vjust = 0.5),
            size = 3.3, colour = "white", fontface = "bold") +
  scale_fill_manual(values = PAL_PAM50, name = "PAM50") +
  scale_y_continuous(labels = function(z) paste0(z * 100, "%")) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(title = "PAM50 subtype composition of each methylation cluster",
       subtitle = sprintf("Numbers = samples; agreement with PAM50 ARI = %.2f, with TCGA clusters ARI = %.2f",
                          ari_pam, ari_mc),
       x = "Methylation cluster (this analysis)", y = "Share of samples") + theme_qc
save_plot(p, "03_23_cluster_pam50_composition.png", 9, 6)

# ============================================================
# 9. GLOBAL METHYLATION BY PAM50 SUBTYPE
# ============================================================
step("9. Global methylation by subtype")
isl_idx <- keep_idx[rd$cgi_context[keep_idx] == "Island"]
sea_idx <- keep_idx[rd$cgi_context[keep_idx] == "OpenSea"]
glob <- data.frame(barcode = cdk$barcode, pam50 = cdk$pam50_label,
                   mean_all     = colMeans2(beta, rows = keep_idx, cols = ks_idx, na.rm = TRUE),
                   mean_island  = colMeans2(beta, rows = isl_idx,  cols = ks_idx, na.rm = TRUE),
                   mean_opensea = colMeans2(beta, rows = sea_idx,  cols = ks_idx, na.rm = TRUE))
known <- glob$pam50 != "Unknown"
kw <- sapply(c("mean_all", "mean_island", "mean_opensea"), function(m)
  if (length(unique(glob$pam50[known])) > 1) kruskal.test(glob[[m]][known], droplevels(glob$pam50[known]))$p.value else NA)
glob_tab <- as.data.frame(data.table(glob)[, .(n = .N,
  median_all = round(median(mean_all), 4), median_island = round(median(mean_island), 4),
  median_opensea = round(median(mean_opensea), 4)), by = pam50][order(pam50)])
cat("Median per-sample methylation by PAM50:\n"); show(glob_tab)
cat(sprintf("Kruskal-Wallis p (differs between subtypes?) - all: %s | islands: %s | open sea: %s\n",
            fmt_p(kw[1]), fmt_p(kw[2]), fmt_p(kw[3])))
save_table(glob_tab, "03_15_global_methylation_by_pam50.csv")

# Fig 24: global methylation per subtype, three genomic compartments side by side.
g_long <- rbind(
  data.frame(pam50 = glob$pam50, value = glob$mean_all,
             metric = sprintf("All retained CpGs\n(Kruskal p = %s)", fmt_p(kw[1]))),
  data.frame(pam50 = glob$pam50, value = glob$mean_island,
             metric = sprintf("CpG islands\n(Kruskal p = %s)", fmt_p(kw[2]))),
  data.frame(pam50 = glob$pam50, value = glob$mean_opensea,
             metric = sprintf("Open sea\n(Kruskal p = %s)", fmt_p(kw[3]))))
g_long$metric <- factor(g_long$metric, levels = unique(g_long$metric))
p <- ggplot(g_long, aes(pam50, value, fill = pam50)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85) +
  geom_jitter(width = 0.15, height = 0, size = 0.7, alpha = 0.45, colour = "grey20") +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = PAL_PAM50, guide = "none") +
  labs(title = "Average methylation of each tumour, by PAM50 subtype",
       subtitle = "Island hypermethylation = CIMP-like pattern; open-sea loss = global hypomethylation",
       x = NULL, y = "Mean beta per sample") +
  theme_qc + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p, "03_24_global_methylation_by_pam50.png", 13, 5.8)

# ============================================================
# 10. SAVE THE CLEANED OBJECT
# ============================================================
step("10. Save cleaned object")

# Sample annotation for the kept samples, plus QC metrics, cluster and PCs.
cd_out <- cbind(cd[ks_idx, ],
                sample_qc[ks_idx, c("pct_missing", "pct_missing_incl_systematic",
                                    "mean_beta", "median_beta", "pct_unmethylated",
                                    "pct_methylated", "cor_to_median", "qc_status", "qc_reason")],
                meth_cluster_03 = clus,
                glob[, c("mean_all", "mean_island", "mean_opensea")],
                pc_df[, paste0("PC", seq_len(min(5, n_pc)))])
rownames(cd_out) <- cd$barcode[ks_idx]

# Probe annotation for the kept probes, plus per-probe statistics.
rd_out <- cbind(rd[keep_idx, ], mean_beta = p_mean, median_beta = p_med, sd_beta = p_sd,
                pct_missing = p_na)

# Free everything not needed before the one large copy.
rm(X, cmat, dens_mat, pca, probe_na_design); invisible(gc())
cat("Creating filtered matrix:", fmt_n(length(keep_idx)), "probes x", length(ks_idx), "samples...\n")
beta_clean <- beta[keep_idx, ks_idx, drop = FALSE]        # the only full-size copy
rm(beta); invisible(gc())

se_clean <- SummarizedExperiment(assays = list(beta = beta_clean), rowData = rd_out, colData = cd_out)
S4Vectors::metadata(se_clean) <- list(
  description   = "QC-filtered TCGA-BRCA HM450 beta values (script 03)",
  input         = INPUT,
  thresholds    = list(SAMPLE_NA_MAX = SAMPLE_NA_MAX, PROBE_NA_MAX = PROBE_NA_MAX,
                       OUTLIER_Z = OUTLIER_Z, TOP_VAR_N = TOP_VAR_N),
  probe_filters = waterfall,
  note          = "Replicate samples still included - collapse to one sample per patient in script 04",
  created       = format(Sys.time()), R_version = R.version.string)

fp <- sum(beta_clean[, 1], na.rm = TRUE)                   # fingerprint to check the saved file
dims_saved <- dim(se_clean)
saveRDS(se_clean, OUT_RDS, compress = FALSE); record(OUT_RDS, "processed data")
cat(sprintf("Saved %s (%.1f MB)\n", OUT_RDS, file.size(OUT_RDS) / 1024^2))
rm(se_clean, beta_clean); invisible(gc())

chk <- readRDS(OUT_RDS)
ok_reload <- identical(dim(chk), dims_saved) && isTRUE(all.equal(sum(assay(chk)[, 1], na.rm = TRUE), fp))
cat("Reload check - dimensions", paste(dim(chk), collapse = " x "), "| values match:", ok_reload, "\n")
if (!ok_reload) stop("Saved cleaned object failed the reload check.")
rm(chk); invisible(gc())

# ============================================================
# 11. KEY INSIGHTS AND OUTPUT INDEX
# ============================================================
step("11. Key insights")
top_pc1 <- assoc[PC == "PC1"][order(-adj_r2)][1]
isl_med <- beta_ctx[grouping == "CpG-island context" & level == "Island", median_of_probe_means]
sea_med <- beta_ctx[grouping == "CpG-island context" & level == "OpenSea", median_of_probe_means]
known_g <- glob_tab$pam50 != "Unknown"
hi_sub  <- as.character(glob_tab$pam50[known_g][which.max(glob_tab$median_island[known_g])])

cat("\n====================== KEY INSIGHTS ======================\n")
cat(sprintf("* Samples: %d of %d passed QC (%d flagged but kept, %d failed and removed).\n",
            length(ks_idx), n_s, sum(status == "Flagged"), sum(status == "Fail")))
cat(sprintf("* Probes: %s of %s retained (%.1f%%) after design, mask and missingness filters.\n",
            fmt_n(length(keep_idx)), fmt_n(n_p), 100 * length(keep_idx) / n_p))
cat(sprintf("* Biology check: CpG islands median beta %.2f vs open sea %.2f - the expected pattern.\n",
            isl_med, sea_med))
if (!is.null(rep_pairs) && nrow(rep_pairs) > 0)
  cat(sprintf("* Reproducibility: replicate pairs r = %.3f (median) vs unrelated tumours r = %.3f.\n",
              median(rep_pairs$r), median(unrel$r)))
cat(sprintf("* PCA: PC1 explains %.1f%% of variance; its strongest factor is '%s' (adj. R2 = %.2f).\n",
            var_exp[1], as.character(top_pc1$covariate), top_pc1$adj_r2))
cat(sprintf("* Clustering: agreement with published TCGA methylation clusters ARI = %.2f; with PAM50 ARI = %.2f.\n",
            ari_mc, ari_pam))
cat(sprintf("* Subtypes: highest median CpG-island methylation in %s; the difference between subtypes is %s (Kruskal p = %s).\n",
            hi_sub, ifelse(!is.na(kw[2]) && kw[2] < 0.05, "statistically significant", "NOT statistically significant"),
            fmt_p(kw[2])))
cat("* Next: script 04 collapses the replicate samples to one per patient before merging with Month 1.\n")
cat("==========================================================\n")

idx <- do.call(rbind, outputs_written)
idx$size_kb <- round(file.size(idx$file) / 1024, 1)
save_table(idx, "03_99_outputs_index.csv")
cat("\nFiles written:", nrow(idx), "(", sum(idx$type == "figure"), "figures,",
    sum(idx$type == "table"), "tables )\n")
cat("Total run time:", fmt_time(as.numeric(difftime(Sys.time(), run_start, units = "secs"))), "\n")
cat("Report saved to:", REPORT, "\n")
cat("\nScript 03 complete.\n")

sink()
close(report_con)
