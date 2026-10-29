#!/usr/bin/env Rscript
# ============================================================
# 05_integrate_epigenomics.R
#
# Integrates the methylation layer with the Month 1 clinical-genomic data
# to build the merged clinical-genomic-epigenomic dataset.
# (Internship brief: "13_integrate_epigenomics"; Week 3.)
#
# WHY METHYLATION IS SUMMARISED FIRST
#   380,811 probe values per patient cannot sit in a patient table. They are
#   condensed into a small set of meaningful, patient-level features:
#     - global methylation in each genomic context (islands, shores, ...)
#     - the main methylation patterns (principal components, PC1-PC10)
#     - the methylation cluster from script 03 (MC1-MC5)
#     - core-promoter methylation of 11 breast-cancer genes, plus a
#       BRCA1 promoter hypermethylation flag
#
# HOW A GENE'S PROMOTER IS DEFINED (tested carefully - see Part J notes)
#   1. The gene's MAIN protein-coding start site (TSS) = the TSS shared by
#      the most protein-coding transcripts (GENCODE v36, reconstructed from
#      the Zhou lab probe manifest).
#   2. Core promoter = 500 bp upstream to 300 bp downstream of that TSS.
#   3. Feature = mean beta of all quality-passed probes in that window
#      (genes with fewer than 2 probes are skipped).
#   A wider window or "island" labels were rejected: for BRCA1 they picked up
#   a neighbouring, normally-methylated CpG island (NBR2) and made BRCA1 look
#   methylated in 97% of tumours instead of the true ~2%.
#
# WHAT THIS SCRIPT DOES, IN ORDER
#   1. Loads the patient-level methylation object (script 04) and Month 1
#   2. Builds the methylation features listed above
#   3. Merges them with the Month 1 table and checks the merge
#   4. Cross-omics analyses:
#        a. BRCA1 promoter methylation vs BRCA1 mutation
#        b. Methylation clusters vs gene mutations
#        c. Methylation clusters vs clinical features
#        d. Survival by methylation cluster (Kaplan-Meier, Cox)
#        e. Global methylation vs age and mutation burden
#        f. Gene promoter methylation by PAM50 subtype
#   5. Saves the merged datasets, feature dictionary, tables and figures
#
# INPUT : data/processed/methylation_beta_matrix_patient.rds    (script 04)
#         data/processed/matched_patient_ids.csv                (script 04)
#         ../01_clinical_genomics/data/processed/clinical_genomics_merged.csv
#         data/reference/HM450.hg38.manifest.gencode.v36.tsv.gz (script 02b)
# OUTPUT: data/processed/clinical_genomic_epigenomic_merged.csv  (671 patients)
#         data/processed/clinical_genomic_epigenomic_all.csv     (all 964)
#         data/processed/methylation_patient_features.csv        (features only)
#         outputs/figures/05_integration/  - 7 figures (PNG, 300 dpi)
#         outputs/tables/05_integration/   - summary tables (CSV)
#         outputs/reports/05_integration_report.txt - everything printed
#
# MEMORY: ~2.5 GB (loads the 2 GB patient-level matrix; no full copies).
#
# Run from the Ubuntu terminal with:
#   Rscript scripts/05_integrate_epigenomics.R
# Test run whose outputs get a suffix (nothing real is overwritten):
#   Rscript scripts/05_integrate_epigenomics.R --tag=test
# ============================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)  # the methylation container
  library(matrixStats)           # fast statistics on the big matrix, no copies
  library(data.table)            # fast tables and CSV writing
  library(ggplot2)               # plots
  library(pheatmap)              # clustered heatmap
  library(survival)              # Kaplan-Meier curves and Cox models
})

# ============================================================
# 0. SETTINGS, PATHS AND HELPER FUNCTIONS
# ============================================================
setwd("/mnt/e/R_Studio/Projects/MOI/Month_2")

args    <- commandArgs(trailingOnly = TRUE)
getflag <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}
TAG <- getflag("tag", "")
sfx <- if (nzchar(TAG)) paste0("_", TAG) else ""

# ---- Inputs ----
INPUT         <- getflag("input", "data/processed/methylation_beta_matrix_patient.rds")
MASTER_IDS    <- "data/processed/matched_patient_ids.csv"
MONTH1_MERGED <- "../01_clinical_genomics/data/processed/clinical_genomics_merged.csv"
GENCODE_MANIFEST <- "data/reference/HM450.hg38.manifest.gencode.v36.tsv.gz"

# ---- Outputs ----
FIG_DIR   <- paste0("outputs/figures/05_integration", sfx)
TAB_DIR   <- paste0("outputs/tables/05_integration", sfx)
REPORT    <- paste0("outputs/reports/05_integration_report", sfx, ".txt")
OUT_MERGED   <- paste0("data/processed/clinical_genomic_epigenomic_merged", sfx, ".csv")
OUT_ALL      <- paste0("data/processed/clinical_genomic_epigenomic_all", sfx, ".csv")
OUT_FEATURES <- paste0("data/processed/methylation_patient_features", sfx, ".csv")
for (d in c(FIG_DIR, TAB_DIR, dirname(REPORT))) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Analysis settings ----
PROMOTER_GENES <- c("BRCA1", "ESR1", "PGR", "GSTP1", "RASSF1", "CDH1",
                    "CDKN2A", "MLH1", "APC", "RARB", "PITX2", "FOXA1")
UPSTREAM_BP   <- 500     # core promoter: 500 bp upstream of the main TSS ...
DOWNSTREAM_BP <- 300     # ... to 300 bp downstream
MIN_PROBES    <- 2       # a gene needs at least 2 probes in its window
HYPER_BETA    <- 0.2     # BRCA1 promoter "hypermethylated" above this (clear gap in the data)
TOP_VAR_N     <- 5000    # most variable probes for the patient-level PCA
N_PC          <- 10      # principal components kept as features
MIN_MUTATED   <- 10      # genes need >= 10 mutated patients to be tested against clusters
SEED          <- 42
set.seed(SEED)

# ---- Everything printed goes to the screen AND the report file ----
report_con <- file(REPORT, open = "wt")
sink(report_con, split = TRUE)

run_start <- Sys.time()
fmt_time  <- function(secs) if (secs < 60) sprintf("%.0f s", secs) else sprintf("%.1f min", secs / 60)
step <- function(title) {
  el <- as.numeric(difftime(Sys.time(), run_start, units = "secs"))
  cat(sprintf("\n[%s | %s] ================ %s ================\n",
              format(Sys.time(), "%H:%M:%S"), fmt_time(el), title))
}
outputs_written <- list()
record <- function(path, what) {
  outputs_written[[length(outputs_written) + 1]] <<- data.frame(type = what, file = path, stringsAsFactors = FALSE)
}
theme_qc <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(colour = "grey30", size = 10.5),
        legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"), plot.margin = margin(10, 15, 10, 10))
save_plot <- function(p, name, width, height) {
  path <- file.path(FIG_DIR, name)
  ggsave(path, p, width = width, height = height, dpi = 300, bg = "white")
  record(path, "figure"); cat("  saved figure:", path, "\n")
}
save_table <- function(df, name, dir = TAB_DIR) {
  path <- file.path(dir, name); fwrite(df, path); record(path, "table"); cat("  saved table :", path, "\n")
}
show  <- function(df) print(as.data.frame(df), row.names = FALSE)
fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-4, formatC(p, format = "e", digits = 1),
                                                   formatC(p, format = "f", digits = 4)))
stars <- function(q) ifelse(is.na(q), "", ifelse(q < 0.001, "***", ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", ""))))

PAM50_LEVELS <- c("LumA", "LumB", "Her2", "Basal", "Normal", "Unknown")
PAL_PAM50 <- c(LumA = "#2C7BB6", LumB = "#74ADD1", Her2 = "#F46D43",
               Basal = "#D73027", Normal = "#1A9850", Unknown = "grey70")
PAL_CLUSTER <- c(MC1 = "#404040", MC2 = "#8C510A", MC3 = "#E6AB02", MC4 = "#08306B",
                 MC5 = "#01665E", MC6 = "#B8E186", MC7 = "#F1B6DA", MC8 = "#DFC27D")   # same as script 03

cat("Script 05 - clinical-genomic-epigenomic integration\n")
cat("Tag    :", ifelse(nzchar(TAG), TAG, "(none - real run)"), "\n")
cat("Started:", format(run_start), "\n")

# ============================================================
# 1. LOAD THE DATA
# ============================================================
step("1. Load data")
if (!file.exists(INPUT)) stop("Input not found: ", INPUT, "\nRun scripts/04_match_harmonise_ids.R first.")
se <- readRDS(INPUT)
B  <- assay(se, withDimnames = FALSE)     # the beta matrix (no copy; see script 04 note)
rd <- as.data.frame(rowData(se))          # probe annotation
cd <- as.data.frame(colData(se))          # one row per patient
patients <- colnames(se)
cat("Methylation:", fmt_n(nrow(se)), "probes x", ncol(se), "patients\n")

m1 <- fread(MONTH1_MERGED)
cat("Month 1   :", nrow(m1), "patients x", ncol(m1), "columns\n")
stopifnot(!any(duplicated(m1$patient_id)), all(patients %in% m1$patient_id))
cat("All methylation patients are present in Month 1 (checked).\n")

# ============================================================
# 2. BUILD PATIENT-LEVEL METHYLATION FEATURES
# ============================================================
step("2. Methylation features")

# ---- 2a. Global methylation in each genomic context ----
# colMeans2(..., rows = ) averages each patient's column over a probe subset
# without copying the matrix (verified against plain R in script 03).
ctx_rows <- function(col, value) which(rd[[col]] == value)
feat <- data.frame(patient_id = patients,
  meth_global_all        = colMeans2(B, na.rm = TRUE),
  meth_global_island     = colMeans2(B, rows = ctx_rows("cgi_context", "Island"),  na.rm = TRUE),
  meth_global_shore      = colMeans2(B, rows = ctx_rows("cgi_context", "Shore"),   na.rm = TRUE),
  meth_global_shelf      = colMeans2(B, rows = ctx_rows("cgi_context", "Shelf"),   na.rm = TRUE),
  meth_global_opensea    = colMeans2(B, rows = ctx_rows("cgi_context", "OpenSea"), na.rm = TRUE),
  meth_global_promoter   = colMeans2(B, rows = ctx_rows("region", "Promoter (TSS +/-1500)"), na.rm = TRUE),
  meth_global_genebody   = colMeans2(B, rows = ctx_rows("region", "Gene body"),  na.rm = TRUE),
  meth_global_intergenic = colMeans2(B, rows = ctx_rows("region", "Intergenic"), na.rm = TRUE),
  stringsAsFactors = FALSE)
cat("Global methylation features: 8 (all CpGs, 4 island contexts, 3 genomic regions)\n")

# ---- 2b. Methylation principal components (patient level) ----
psd     <- rowSds(B, na.rm = TRUE)                               # variability of every probe
top_idx <- order(psd, decreasing = TRUE)[seq_len(TOP_VAR_N)]
X       <- B[top_idx, , drop = FALSE]                            # small copy: 5,000 x patients
na_pos  <- which(is.na(X), arr.ind = TRUE)
if (nrow(na_pos)) X[na_pos] <- rowMeans(X, na.rm = TRUE)[na_pos[, 1]]   # fill gaps with the PROBE mean
pca     <- prcomp(t(X), center = TRUE, scale. = FALSE)
var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pcs     <- as.data.frame(pca$x[, seq_len(N_PC)]); names(pcs) <- paste0("meth_PC", seq_len(N_PC))
feat    <- cbind(feat, pcs)
pc_var  <- data.frame(PC = paste0("meth_PC", seq_len(N_PC)), pct_variance = round(var_exp[seq_len(N_PC)], 2),
                      cumulative_pct = round(cumsum(var_exp)[seq_len(N_PC)], 2))
cat(sprintf("Methylation PCs: %d kept; PC1 %.1f%%, PC2 %.1f%%; first %d explain %.1f%%\n",
            N_PC, var_exp[1], var_exp[2], N_PC, sum(var_exp[seq_len(N_PC)])))
save_table(pc_var, "05_01_methylation_pc_variance.csv")
rm(X, na_pos); invisible(gc())

# ---- 2c. Cluster, subtype and QC covariates carried over from scripts 02b-04 ----
feat$meth_cluster      <- cd$meth_cluster_03                      # MC1-MC5 (script 03)
feat$tcga_meth_cluster <- ifelse(cd$meth_cluster_tcga %in% c("NA", ""), NA, cd$meth_cluster_tcga)
feat$pam50_subtype     <- ifelse(cd$pam50 %in% c("NA", ""), NA, cd$pam50)
feat$meth_barcode      <- cd$methylation_barcode
feat$meth_plate        <- cd$plate                                # processing batch
feat$meth_qc_status    <- cd$qc_status
feat$meth_pct_missing  <- round(cd$pct_missing, 3)

# ---- 2d. Core-promoter methylation of key breast-cancer genes ----
# Step 1: for each gene, rebuild every protein-coding transcript's TSS from
# the manifest (TSS = probe position -/+ distance, depending on strand; this
# convention was verified: it gives one constant TSS per transcript in
# 99.98% of transcripts).
man <- fread(GENCODE_MANIFEST, select = c("CpG_chrm", "CpG_beg", "geneNames", "transcriptTypes",
                                          "transcriptIDs", "distToTSS"))
man <- man[grepl(paste0("(^|;)(", paste(PROMOTER_GENES, collapse = "|"), ")(;|$)"), geneNames)]
tx  <- man[, .(gene = unlist(strsplit(geneNames, ";", fixed = TRUE)),
               ttype = unlist(strsplit(transcriptTypes, ";", fixed = TRUE)),
               tid = unlist(strsplit(transcriptIDs, ";", fixed = TRUE)),
               dist = as.integer(unlist(strsplit(distToTSS, ";", fixed = TRUE)))),
           by = .(chr = CpG_chrm, pos = CpG_beg + 1L)]
tx  <- tx[gene %in% PROMOTER_GENES & ttype == "protein_coding"]
data("gencode.v36.annotation.genes", package = "TCGAbiolinksGUI.data", envir = environment())
strand_tab <- as.data.table(as.data.frame(gencode.v36.annotation.genes))[, .(gene = gene_name, strand = as.character(strand))]
strand_tab <- strand_tab[!duplicated(gene)]
tx  <- merge(tx, strand_tab, by = "gene")
tx[, tss := ifelse(strand == "+", pos - dist, pos + dist)]
tx  <- unique(tx[, .(gene, chr, strand, tid, tss)])

# Step 2: the main TSS = the group of TSSs (within 200 bp) shared by the most
# protein-coding transcripts.
major <- tx[, {
  t <- sort(tss); grp <- cumsum(c(1, diff(t) > 200)); cnt <- tabulate(grp); best <- which.max(cnt)
  .(main_tss = as.integer(median(t[grp == best])), transcripts_at_main_tss = cnt[best],
    protein_coding_transcripts = length(t))
}, by = .(gene, chr, strand)]

# Step 3: probes in the core window; feature = their mean beta per patient.
prom_info <- list()
for (k in seq_len(nrow(major))) {
  x   <- major[k]
  win <- if (x$strand == "+") c(x$main_tss - UPSTREAM_BP, x$main_tss + DOWNSTREAM_BP) else
                              c(x$main_tss - DOWNSTREAM_BP, x$main_tss + UPSTREAM_BP)
  pr  <- which(rd$chr == x$chr & rd$pos_hg38 >= win[1] & rd$pos_hg38 <= win[2])
  used <- length(pr) >= MIN_PROBES
  if (used) feat[[paste0("meth_prom_", x$gene)]] <- colMeans2(B, rows = pr, na.rm = TRUE)
  prom_info[[x$gene]] <- data.frame(gene = x$gene, chr = x$chr, strand = x$strand, main_tss = x$main_tss,
    transcripts_at_main_tss = x$transcripts_at_main_tss, protein_coding_transcripts = x$protein_coding_transcripts,
    window = sprintf("%s:%d-%d", x$chr, win[1], win[2]), n_probes = length(pr),
    used = ifelse(used, "yes", paste0("no (<", MIN_PROBES, " probes)")),
    probes = paste(rd$probe_id[pr], collapse = ";"),
    median_beta = if (used) round(median(feat[[paste0("meth_prom_", x$gene)]]), 3) else NA)
}
prom_info <- rbindlist(prom_info)[order(gene)]
cat("\nCore-promoter windows (main protein-coding TSS, -", UPSTREAM_BP, "/+", DOWNSTREAM_BP, " bp):\n", sep = "")
show(prom_info[, .(gene, window, transcripts_at_main_tss, n_probes, used, median_beta)])
save_table(prom_info, "05_02_gene_promoter_definitions.csv")
prom_genes <- prom_info$gene[prom_info$used == "yes"]
prom_cols  <- paste0("meth_prom_", prom_genes)

feat$meth_BRCA1_hypermethylated <- as.integer(feat$meth_prom_BRCA1 > HYPER_BETA)
cat(sprintf("\nBRCA1 promoter hypermethylated (beta > %.1f): %d of %d patients (%.1f%%)\n",
            HYPER_BETA, sum(feat$meth_BRCA1_hypermethylated), nrow(feat), 100 * mean(feat$meth_BRCA1_hypermethylated)))

num_feat <- setdiff(names(feat)[sapply(feat, is.numeric)], "meth_pct_missing")
feat[num_feat] <- lapply(feat[num_feat], function(v) round(v, 5))
fwrite(feat, OUT_FEATURES); record(OUT_FEATURES, "processed data")
cat("  saved features:", OUT_FEATURES, "(", nrow(feat), "patients x", ncol(feat) - 1, "features )\n")

# ---- Feature dictionary: what every new column means ----
dict <- rbind(
  data.frame(feature = "meth_global_all", description = "Mean beta over all QC-passed CpGs", type = "numeric 0-1"),
  data.frame(feature = paste0("meth_global_", c("island", "shore", "shelf", "opensea")),
             description = paste("Mean beta over", c("CpG-island", "shore", "shelf", "open-sea"), "CpGs"), type = "numeric 0-1"),
  data.frame(feature = paste0("meth_global_", c("promoter", "genebody", "intergenic")),
             description = paste("Mean beta over", c("promoter (TSS +/-1500)", "gene-body", "intergenic"), "CpGs"), type = "numeric 0-1"),
  data.frame(feature = paste0("meth_PC", seq_len(N_PC)),
             description = paste0("Methylation principal component ", seq_len(N_PC), " (", round(var_exp[seq_len(N_PC)], 1),
                                   "% of variance; ", fmt_n(TOP_VAR_N), " most variable probes)"), type = "numeric"),
  data.frame(feature = "meth_cluster", description = "Unsupervised methylation cluster from script 03 (MC1-MC5)", type = "category"),
  data.frame(feature = "tcga_meth_cluster", description = "Published TCGA methylation cluster (Berger et al. 2018)", type = "category"),
  data.frame(feature = "pam50_subtype", description = "PAM50 subtype (TCGA)", type = "category"),
  data.frame(feature = c("meth_barcode", "meth_plate", "meth_qc_status", "meth_pct_missing"),
             description = c("Barcode of the methylation sample used", "Processing plate (batch) of that sample",
                             "QC status from script 03 (Pass / Flagged)", "Sample-specific missing values (%)"),
             type = c("text", "category", "category", "numeric")),
  data.frame(feature = prom_cols,
             description = paste0("Core-promoter methylation of ", prom_genes, " (mean beta, ",
                                  prom_info$n_probes[match(prom_genes, prom_info$gene)], " probes)"), type = "numeric 0-1"),
  data.frame(feature = "meth_BRCA1_hypermethylated",
             description = sprintf("1 if BRCA1 core-promoter beta > %.1f (promoter silenced)", HYPER_BETA), type = "0/1"))
save_table(dict, "05_03_feature_dictionary.csv")

# ============================================================
# 3. MERGE WITH THE MONTH 1 CLINICAL-GENOMIC TABLE
# ============================================================
step("3. Merge and harmonise")
clash <- intersect(setdiff(names(feat), "patient_id"), names(m1))
if (length(clash)) stop("Column name clash between methylation features and Month 1: ", paste(clash, collapse = ", "))

# Negative follow-up time is impossible - set to missing rather than keep an error.
n_neg <- sum(m1$OS_MONTHS < 0, na.rm = TRUE)
if (n_neg) { m1$OS_MONTHS[m1$OS_MONTHS < 0] <- NA; cat("OS_MONTHS < 0 set to missing for", n_neg, "patient(s)\n") }

master <- fread(MASTER_IDS)[, .(patient_id, has_methylation, methylation_status)]
all_tab <- merge(m1, master, by = "patient_id", all.x = TRUE, sort = FALSE)
all_tab <- merge(all_tab, as.data.table(feat), by = "patient_id", all.x = TRUE, sort = FALSE)
merged  <- all_tab[has_methylation == TRUE]

merge_summary <- data.frame(
  item  = c("Month 1 patients", "Month 1 columns", "Methylation features added", "Patients with methylation (merged table)",
            "Patients without methylation (kept in the 'all' table)", "Duplicated patient IDs after merging",
            "Merged table columns"),
  value = c(nrow(m1), ncol(m1), ncol(feat) - 1, nrow(merged), nrow(all_tab) - nrow(merged),
            sum(duplicated(all_tab$patient_id)), ncol(merged)))
show(merge_summary)
save_table(merge_summary, "05_04_merge_summary.csv")
stopifnot(nrow(all_tab) == nrow(m1), nrow(merged) == length(patients), !any(duplicated(all_tab$patient_id)))

fwrite(merged, OUT_MERGED);  record(OUT_MERGED, "processed data")
fwrite(all_tab, OUT_ALL);    record(OUT_ALL, "processed data")
n_cols_saved <- ncol(merged)            # helper columns added below are NOT in the saved file
cat("  saved merged dataset:", OUT_MERGED, "(", nrow(merged), "x", ncol(merged), ")\n")
cat("  saved full dataset  :", OUT_ALL, "(", nrow(all_tab), "x", ncol(all_tab), ")\n")

# Completeness of the key variables in the merged table, by data layer.
key_vars <- data.frame(
  variable = c("AGE", "stage_simple", "os_event", "OS_MONTHS", "DFS_MONTHS", "tmb_nonsynonymous",
               "TP53_mutated", "PIK3CA_mutated", "BRCA1_mutated", "pam50_subtype", "tcga_meth_cluster",
               "meth_cluster", "meth_global_island", "meth_PC1", "meth_prom_BRCA1", "meth_prom_GSTP1"),
  layer = c(rep("Clinical", 5), rep("Genomic", 4), "Subtype", rep("Epigenomic", 6)))
key_vars <- key_vars[key_vars$variable %in% names(merged), ]
key_vars$pct_complete <- sapply(key_vars$variable, function(v) {
  x <- merged[[v]]; round(100 * mean(!is.na(x) & !(as.character(x) %in% c("", "NA"))), 1) })
save_table(key_vars, "05_05_variable_completeness.csv")

key_vars$variable <- factor(key_vars$variable, levels = rev(key_vars$variable))
p <- ggplot(key_vars, aes(pct_complete, variable, fill = layer)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(pct_complete, "%")), hjust = -0.1, size = 3.3) +
  scale_fill_manual(values = c(Clinical = "#0F6E56", Genomic = "#8C510A", Subtype = "#74ADD1", Epigenomic = "#185FA5"),
                    name = "Data layer") +
  scale_x_continuous(limits = c(0, 115), breaks = seq(0, 100, 25)) +
  labs(title = "The merged clinical-genomic-epigenomic dataset",
       subtitle = sprintf("%d patients x %d columns; share of patients with each key variable available",
                          nrow(merged), ncol(merged)),
       x = "Patients with a value (%)", y = NULL) + theme_qc
save_plot(p, "05_01_merged_dataset_completeness.png", 9, 6.5)

# Analysis-ready factors used below.
merged[, meth_cluster := factor(meth_cluster, levels = sort(unique(meth_cluster)))]
merged[, pam50_label := factor(ifelse(is.na(pam50_subtype), "Unknown", pam50_subtype), levels = PAM50_LEVELS)]
merged[, histology := ifelse(PRIMARY_DIAGNOSIS == "Infiltrating Ductal Carcinoma", "Ductal",
                      ifelse(PRIMARY_DIAGNOSIS == "Infiltrating Lobular Carcinoma", "Lobular", "Other"))]

# ============================================================
# 4a. BRCA1: PROMOTER METHYLATION vs MUTATION
#     Two ways to switch BRCA1 off: a damaging mutation or promoter
#     methylation. They are expected to be (almost) mutually exclusive.
# ============================================================
step("4a. BRCA1 methylation vs mutation")
tb_brca1 <- table(Mutated = factor(merged$BRCA1_mutated == 1, levels = c(FALSE, TRUE), labels = c("No", "Yes")),
                  Hypermethylated = factor(merged$meth_BRCA1_hypermethylated == 1, levels = c(FALSE, TRUE),
                                           labels = c("No", "Yes")))
fisher_brca1 <- fisher.test(tb_brca1)
cat("BRCA1 mutation x BRCA1 promoter hypermethylation:\n"); print(tb_brca1)
cat(sprintf("Fisher's exact test p = %s (odds ratio %.2f)\n", fmt_p(fisher_brca1$p.value), fisher_brca1$estimate))
by_sub <- merged[, .(patients = .N, brca1_hypermethylated = sum(meth_BRCA1_hypermethylated),
                     pct = round(100 * mean(meth_BRCA1_hypermethylated), 1)), by = pam50_label][order(pam50_label)]
cat("\nBRCA1 promoter hypermethylation by PAM50:\n"); show(by_sub)
brca1_out <- rbind(
  data.frame(group = "BRCA1 mutated", hypermethylated = tb_brca1["Yes", "Yes"], total = sum(tb_brca1["Yes", ])),
  data.frame(group = "BRCA1 not mutated", hypermethylated = tb_brca1["No", "Yes"], total = sum(tb_brca1["No", ])),
  data.frame(group = paste("PAM50", by_sub$pam50_label), hypermethylated = by_sub$brca1_hypermethylated, total = by_sub$patients))
brca1_out$pct <- round(100 * brca1_out$hypermethylated / brca1_out$total, 1)
brca1_out$fisher_p_mutation_vs_methylation <- c(signif(fisher_brca1$p.value, 3), rep(NA, nrow(brca1_out) - 1))
save_table(brca1_out, "05_06_brca1_methylation_vs_mutation.csv")

# Fig 3: BRCA1 promoter methylation per patient, by subtype, mutated patients highlighted.
bd <- data.frame(pam50 = merged$pam50_label, beta = merged$meth_prom_BRCA1,
                 status = ifelse(merged$BRCA1_mutated == 1, "BRCA1 mutated", "BRCA1 not mutated"))
bd <- bd[order(bd$status == "BRCA1 mutated"), ]            # draw mutated patients on top
# The points are drawn before the text label so the x-axis is set up as
# categories (subtypes); a numeric label first would make it a number axis.
p <- ggplot(bd, aes(pam50, beta)) +
  geom_hline(yintercept = HYPER_BETA, linetype = "dashed", colour = "#EF4444") +
  geom_jitter(aes(colour = status, size = status), width = 0.25, height = 0, alpha = 0.75) +
  annotate("text", x = 0.55, y = HYPER_BETA, label = paste("hypermethylated above", HYPER_BETA),
           hjust = 0, vjust = -0.6, size = 3.4, colour = "#EF4444") +
  scale_x_discrete() +
  scale_colour_manual(values = c("BRCA1 not mutated" = "grey55", "BRCA1 mutated" = "#EF4444"), name = NULL) +
  scale_size_manual(values = c("BRCA1 not mutated" = 1.1, "BRCA1 mutated" = 2.6), name = NULL) +
  labs(title = "BRCA1 promoter methylation versus BRCA1 mutation",
       subtitle = sprintf(paste0("Hypermethylated: %d tumours, %d of them Basal. Only %d of %d BRCA1-mutated tumours is also\n",
                                 "methylated - the two rarely co-occur, but cases are too few for a formal test (Fisher p = %s)"),
                          sum(merged$meth_BRCA1_hypermethylated), by_sub$brca1_hypermethylated[by_sub$pam50_label == "Basal"],
                          tb_brca1["Yes", "Yes"], sum(tb_brca1["Yes", ]), fmt_p(fisher_brca1$p.value)),
       x = "PAM50 subtype", y = "BRCA1 core-promoter methylation (beta)") + theme_qc
save_plot(p, "05_03_brca1_methylation_vs_mutation.png", 9.5, 6)

# ============================================================
# 4b. METHYLATION CLUSTERS vs GENE MUTATIONS
# ============================================================
step("4b. Clusters vs mutations")
mut_cols <- grep("_mutated$", names(merged), value = TRUE)
mut_cols <- mut_cols[sapply(mut_cols, function(v) sum(merged[[v]] == 1, na.rm = TRUE) >= MIN_MUTATED)]
cl_mut <- rbindlist(lapply(mut_cols, function(v) {
  x  <- merged[[v]] == 1
  tb <- table(merged$meth_cluster, x)
  p  <- tryCatch(suppressWarnings(chisq.test(tb, simulate.p.value = TRUE, B = 5000)$p.value), error = function(e) NA)
  pct <- tapply(x, merged$meth_cluster, function(z) round(100 * mean(z, na.rm = TRUE), 1))
  data.table(gene = sub("_mutated$", "", v), n_mutated = sum(x, na.rm = TRUE),
             overall_pct = round(100 * mean(x, na.rm = TRUE), 1), t(pct), p_value = p)
}))
cl_mut[, q_value := p.adjust(p_value, method = "BH")]
setorder(cl_mut, p_value)
cat("% of patients with each mutation, by methylation cluster (BH-adjusted q):\n")
show(cl_mut[, lapply(.SD, function(z) if (is.numeric(z)) signif(z, 3) else z)])
save_table(cl_mut, "05_07_clusters_vs_mutations.csv")

# Fig 4: mutation frequency heatmap (rows = genes, columns = clusters).
cl_levels <- levels(merged$meth_cluster)
mlong <- melt(cl_mut, id.vars = c("gene", "q_value"), measure.vars = cl_levels, variable.name = "cluster", value.name = "pct")
mlong[, gene_lab := paste0(gene, " ", stars(q_value))]
mlong[, gene_lab := factor(gene_lab, levels = rev(unique(mlong[order(-cl_mut$overall_pct[match(gene, cl_mut$gene)])]$gene_lab)))]
csize <- as.integer(table(merged$meth_cluster)[cl_levels])
mlong[, cluster := factor(cluster, levels = cl_levels, labels = paste0(cl_levels, "\n(n=", csize, ")"))]
p <- ggplot(mlong, aes(cluster, gene_lab, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.0f%%", pct), colour = pct > 0.6 * max(pct)), size = 3.3) +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +
  scale_fill_gradient(low = "white", high = "#B2182B", name = "% mutated") +
  labs(title = "How often each gene is mutated in each methylation cluster",
       subtitle = "Genes sorted by overall frequency; * q<0.05, ** q<0.01, *** q<0.001 (difference between clusters)",
       x = "Methylation cluster (script 03)", y = NULL) +
  theme_qc + theme(panel.grid = element_blank(), legend.position = "right")
save_plot(p, "05_04_mutations_by_methylation_cluster.png", 9, 6.5)

# ============================================================
# 4c. METHYLATION CLUSTERS vs CLINICAL FEATURES
# ============================================================
step("4c. Clusters vs clinical")
kw_row <- function(var, label) {
  x <- suppressWarnings(as.numeric(merged[[var]]))
  med <- tapply(x, merged$meth_cluster, median, na.rm = TRUE)
  data.frame(variable = label, test = "Kruskal-Wallis",
             summary = paste(sprintf("%s %.1f", names(med), med), collapse = "; "),
             p_value = kruskal.test(x ~ merged$meth_cluster)$p.value)
}
cs_row <- function(var, label) {
  x <- as.character(merged[[var]]); tb <- table(merged$meth_cluster, x)
  top <- names(sort(colSums(tb), decreasing = TRUE))[1]
  data.frame(variable = label, test = "Chi-square (simulated p)",
             summary = paste(sprintf("%s %.0f%% %s", rownames(tb), 100 * tb[, top] / rowSums(tb), top), collapse = "; "),
             p_value = suppressWarnings(chisq.test(tb, simulate.p.value = TRUE, B = 5000)$p.value))
}
cl_clin <- rbind(kw_row("AGE", "Age at diagnosis (median years)"),
                 kw_row("tmb_nonsynonymous", "Tumour mutation burden (median)"),
                 cs_row("stage_simple", "Pathological stage"),
                 cs_row("histology", "Histology (ductal / lobular / other)"),
                 cs_row("pam50_label", "PAM50 subtype"))
cl_clin$q_value <- p.adjust(cl_clin$p_value, method = "BH")
show(transform(cl_clin, p_value = signif(p_value, 3), q_value = signif(q_value, 3)))
save_table(cl_clin, "05_08_clusters_vs_clinical.csv")

# ============================================================
# 4d. SURVIVAL BY METHYLATION CLUSTER
# ============================================================
step("4d. Survival by cluster")
sv <- merged[!is.na(OS_MONTHS) & !is.na(os_event)]
fit    <- survfit(Surv(OS_MONTHS, os_event) ~ meth_cluster, data = sv)
lr     <- survdiff(Surv(OS_MONTHS, os_event) ~ meth_cluster, data = sv)
lr_p   <- pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
s60    <- summary(fit, times = 60, extend = TRUE)
surv_tab <- data.frame(cluster = sub("meth_cluster=", "", names(fit$strata)), patients = fit$n,
                       deaths = as.integer(tapply(sv$os_event, sv$meth_cluster, sum)),
                       five_year_survival_pct = round(100 * s60$surv, 1),
                       ci_low = round(100 * s60$lower, 1), ci_high = round(100 * s60$upper, 1))

# Cox model: cluster effect adjusted for age and stage (reference = largest cluster).
ref_cl <- names(which.max(table(sv$meth_cluster)))
sv[, meth_cluster_ref := relevel(droplevels(meth_cluster), ref = ref_cl)]
sv[, stage_group := ifelse(stage_simple %in% c("Stage I", "Stage II"), "I-II",
                    ifelse(stage_simple %in% c("Stage III", "Stage IV"), "III-IV", NA))]
cox <- tryCatch(coxph(Surv(OS_MONTHS, os_event) ~ meth_cluster_ref + AGE + stage_group, data = sv),
                error = function(e) NULL)
if (!is.null(cox)) {
  cs <- summary(cox)$conf.int; cp <- summary(cox)$coefficients
  cox_tab <- data.frame(term = rownames(cs), hazard_ratio = round(cs[, 1], 2), ci_low = round(cs[, 3], 2),
                        ci_high = round(cs[, 4], 2), p_value = signif(cp[, 5], 3))
  cox_tab$term <- sub("meth_cluster_ref", paste0("cluster (vs ", ref_cl, "): "), cox_tab$term)
  cat("Cox model (overall survival ~ cluster + age + stage):\n"); show(cox_tab)
  save_table(cox_tab, "05_10_cox_model_cluster_age_stage.csv")
}
cat(sprintf("\nLog-rank test across clusters: p = %s | deaths in total: %d of %d\n", fmt_p(lr_p), sum(sv$os_event), nrow(sv)))
show(surv_tab)
save_table(surv_tab, "05_09_survival_by_cluster.csv")

# Fig 5: Kaplan-Meier curves (drawn with ggplot; first 10 years shown).
km <- rbindlist(lapply(seq_along(fit$strata), function(i) {
  idx <- if (i == 1) seq_len(fit$strata[1]) else (sum(fit$strata[1:(i - 1)]) + 1):sum(fit$strata[1:i])
  data.table(cluster = sub("meth_cluster=", "", names(fit$strata)[i]),
             time = c(0, fit$time[idx]), surv = c(1, fit$surv[idx]))
}))
km_lab <- setNames(sprintf("%s (n=%d, deaths=%d)", surv_tab$cluster, surv_tab$patients, surv_tab$deaths), surv_tab$cluster)
p <- ggplot(km, aes(time, surv, colour = cluster)) +
  geom_step(linewidth = 0.9) +
  scale_colour_manual(values = PAL_CLUSTER[names(km_lab)], labels = km_lab, name = NULL) +
  scale_y_continuous(limits = c(0, 1), labels = function(z) paste0(z * 100, "%")) +
  coord_cartesian(xlim = c(0, 120)) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(title = "Overall survival by methylation cluster",
       subtitle = sprintf("Kaplan-Meier estimates; log-rank p = %s; only %d deaths in %d patients, so power is limited",
                          fmt_p(lr_p), sum(sv$os_event), nrow(sv)),
       x = "Months since diagnosis (first 10 years shown)", y = "Patients alive") + theme_qc
save_plot(p, "05_05_survival_by_methylation_cluster.png", 9, 6.5)

# ============================================================
# 4e. GLOBAL METHYLATION vs AGE AND MUTATION BURDEN
# ============================================================
step("4e. Methylation vs age and TMB")
cor_tab <- rbindlist(lapply(c("meth_global_island", "meth_global_opensea", "meth_global_all", "meth_PC1"), function(f) {
  rbindlist(lapply(c(AGE = "AGE", TMB = "tmb_nonsynonymous"), function(v) {
    ct <- suppressWarnings(cor.test(merged[[f]], as.numeric(merged[[v]]), method = "spearman"))
    data.table(methylation_feature = f, versus = v, spearman_rho = round(unname(ct$estimate), 3), p_value = ct$p.value)
  }))
}))
cor_tab[, q_value := p.adjust(p_value, method = "BH")]
show(cor_tab[, .(methylation_feature, versus, spearman_rho, p_value = signif(p_value, 3), q_value = signif(q_value, 3))])
save_table(cor_tab, "05_11_methylation_vs_age_tmb.csv")

# Fig 6: age versus island and open-sea methylation.
rho_lab <- function(f) { r <- cor_tab[methylation_feature == f & versus == "AGE"]; sprintf("rho = %.2f, p = %s", r$spearman_rho, fmt_p(r$p_value)) }
al <- rbind(data.frame(age = merged$AGE, beta = merged$meth_global_island, panel = paste0("CpG islands (", rho_lab("meth_global_island"), ")")),
            data.frame(age = merged$AGE, beta = merged$meth_global_opensea, panel = paste0("Open sea (", rho_lab("meth_global_opensea"), ")")))
n_no_age <- sum(is.na(merged$AGE))
al <- al[!is.na(al$age), ]                 # patients without a recorded age cannot be plotted
if (n_no_age) cat("Patients without a recorded age (left out of the plot):", n_no_age, "\n")
p <- ggplot(al, aes(age, beta)) +
  geom_point(size = 1, alpha = 0.45, colour = "#185FA5") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#EF4444", fill = "#F4B6B6") +
  facet_wrap(~ panel, scales = "free_y") +
  labs(title = "Tumour methylation and age at diagnosis",
       subtitle = "Each point is a patient; red line = linear trend with 95% confidence band; Spearman correlation in panel titles",
       x = "Age at diagnosis (years)", y = "Mean beta") + theme_qc
save_plot(p, "05_06_methylation_vs_age.png", 11, 5.5)

# ============================================================
# 4f. GENE PROMOTER METHYLATION BY PAM50 SUBTYPE
# ============================================================
step("4f. Promoter methylation by subtype")
known <- merged$pam50_label != "Unknown"
prom_sub <- rbindlist(lapply(prom_genes, function(g) {
  v <- merged[[paste0("meth_prom_", g)]]
  med <- tapply(v[known], droplevels(merged$pam50_label[known]), median)
  data.table(gene = g, t(round(med, 3)),
             p_value = kruskal.test(v[known] ~ droplevels(merged$pam50_label[known]))$p.value)
}))
prom_sub[, q_value := p.adjust(p_value, method = "BH")]
# Effect size: how far apart the subtype medians are. With 671 patients even
# tiny differences reach significance, so a gene only counts as
# SUBSTANTIALLY different if its subtype medians also differ by >= 0.10.
sub_cols <- intersect(PAM50_LEVELS, names(prom_sub))
prom_sub[, median_range := round(apply(.SD, 1, max) - apply(.SD, 1, min), 3), .SDcols = sub_cols]
prom_sub[, substantial := q_value < 0.05 & median_range >= 0.10]
setorder(prom_sub, -median_range)
cat("Median core-promoter methylation by PAM50 subtype:\n")
show(prom_sub[, lapply(.SD, function(z) if (is.numeric(z)) signif(z, 3) else z)])
save_table(prom_sub, "05_12_promoter_methylation_by_pam50.csv")

# Fig 2: heatmap - genes x patients, patients grouped by PAM50 subtype.
ord <- order(merged$pam50_label, merged$meth_cluster, -merged$meth_prom_BRCA1)
hm  <- t(as.matrix(merged[ord, ..prom_cols])); rownames(hm) <- prom_genes; colnames(hm) <- merged$patient_id[ord]
ann <- data.frame(PAM50 = merged$pam50_label[ord], Cluster = merged$meth_cluster[ord],
                  BRCA1_mut = factor(ifelse(merged$BRCA1_mutated[ord] == 1, "Mutated", "Wild-type"), levels = c("Wild-type", "Mutated")),
                  TP53_mut  = factor(ifelse(merged$TP53_mutated[ord] == 1, "Mutated", "Wild-type"), levels = c("Wild-type", "Mutated")),
                  row.names = colnames(hm))
ann_col <- list(PAM50 = PAL_PAM50[levels(droplevels(ann$PAM50))], Cluster = PAL_CLUSTER[levels(ann$Cluster)],
                BRCA1_mut = c("Wild-type" = "grey92", Mutated = "black"),
                TP53_mut  = c("Wild-type" = "grey92", Mutated = "black"))
ann$PAM50 <- droplevels(ann$PAM50)
gaps <- cumsum(as.integer(table(ann$PAM50)))[-length(levels(ann$PAM50))]
pheatmap(hm, cluster_rows = TRUE, cluster_cols = FALSE, show_colnames = FALSE, gaps_col = gaps,
         annotation_col = ann, annotation_colors = ann_col, border_color = NA, fontsize = 10,
         color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100), breaks = seq(0, 1, length.out = 101),
         main = "Core-promoter methylation of breast-cancer genes (patients grouped by PAM50)",
         filename = file.path(FIG_DIR, "05_02_promoter_methylation_heatmap.png"), width = 12, height = 5.5)
record(file.path(FIG_DIR, "05_02_promoter_methylation_heatmap.png"), "figure")
cat("  saved figure:", file.path(FIG_DIR, "05_02_promoter_methylation_heatmap.png"), "\n")

# Fig 7: the genes whose promoter methylation differs most between subtypes
# (largest spread of subtype medians among significant genes).
top_genes <- head(prom_sub[q_value < 0.05]$gene, 6)
if (length(top_genes)) {
  gl <- rbindlist(lapply(top_genes, function(g)
    data.table(gene = sprintf("%s (q = %s)", g, fmt_p(prom_sub[gene == g]$q_value)),
               pam50 = merged$pam50_label, beta = merged[[paste0("meth_prom_", g)]])))[pam50 != "Unknown"]
  gl[, gene := factor(gene, levels = unique(gene))]
  p <- ggplot(gl, aes(pam50, beta, fill = pam50)) +
    geom_boxplot(outlier.size = 0.6, width = 0.6, alpha = 0.85) +
    facet_wrap(~ gene, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = PAL_PAM50, guide = "none") +
    labs(title = "Promoter methylation that differs between PAM50 subtypes",
         subtitle = "The 6 genes with the largest spread between subtype medians (Kruskal-Wallis, BH-adjusted q in titles).\nNote: each panel has its own y-axis scale.",
         x = NULL, y = "Core-promoter methylation (beta)") +
    theme_qc + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, "05_07_promoter_methylation_by_pam50.png", 11, 7)
}

# ============================================================
# 5. KEY INSIGHTS AND OUTPUT INDEX
# ============================================================
step("5. Key insights")
sig_mut <- cl_mut[q_value < 0.05]$gene
top_mut <- cl_mut[1]
top_cl  <- cl_levels[which.max(unlist(cl_mut[1, ..cl_levels]))]
cat("\n====================== KEY INSIGHTS ======================\n")
cat(sprintf("* Merged dataset: %d patients x %d columns (%d Month 1 columns + %d methylation features + 2 status columns).\n",
            nrow(merged), n_cols_saved, ncol(m1), ncol(feat) - 1))
cat(sprintf("* BRCA1: promoter hypermethylated in %d tumours (%d Basal); only %d of %d BRCA1-mutated tumours is also methylated - rare overlap, too few cases to test formally (Fisher p = %s).\n",
            sum(merged$meth_BRCA1_hypermethylated), by_sub$brca1_hypermethylated[by_sub$pam50_label == "Basal"],
            tb_brca1["Yes", "Yes"], sum(tb_brca1["Yes", ]), fmt_p(fisher_brca1$p.value)))
cat(sprintf("* Mutations: %d of %d tested genes differ in frequency between methylation clusters (q<0.05); strongest: %s (highest in %s, %.0f%%).\n",
            length(sig_mut), nrow(cl_mut), top_mut$gene, top_cl, max(unlist(cl_mut[1, ..cl_levels]))))
cat(sprintf("* Promoters: %d of %d genes differ significantly between PAM50 subtypes, but only %d substantially (medians differ by >= 0.10): %s.\n",
            sum(prom_sub$q_value < 0.05), nrow(prom_sub), sum(prom_sub$substantial),
            paste(prom_sub$gene[prom_sub$substantial], collapse = ", ")))
cat(sprintf("* Survival: log-rank p = %s across clusters (%d deaths); interpret cautiously.\n", fmt_p(lr_p), sum(sv$os_event)))
r_isl <- cor_tab[methylation_feature == "meth_global_island" & versus == "AGE"]
cat(sprintf("* Age: CpG-island methylation vs age rho = %.2f (p = %s).\n", r_isl$spearman_rho, fmt_p(r_isl$p_value)))
cat("* Next: script 06 builds the publication figures; script 07 the AI-ready feature table.\n")
cat("==========================================================\n")

idx <- do.call(rbind, outputs_written)
idx$size_kb <- round(file.size(idx$file) / 1024, 1)
save_table(idx, "05_99_outputs_index.csv")
cat("\nFiles written:", nrow(idx), "(", sum(idx$type == "figure"), "figures,", sum(idx$type == "table"), "tables )\n")
cat("Total run time:", fmt_time(as.numeric(difftime(Sys.time(), run_start, units = "secs"))), "\n")
cat("Report saved to:", REPORT, "\n")
cat("\nScript 05 complete.\n")

sink()
close(report_con)
