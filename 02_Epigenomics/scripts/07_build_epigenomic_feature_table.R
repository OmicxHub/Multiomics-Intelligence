#!/usr/bin/env Rscript
# ============================================================
# 07_build_epigenomic_feature_table.R
#
# Builds the AI-ready feature table for the Month 2 baseline models.
# (Internship brief: "15_build_epigenomic_feature_table"; Week 4.)
#
# TARGET - exactly as in the REVISED Month 1 pipeline
#   (01R_clinical_genomics/scripts/07_build_ml_feature_table.R)
#   five_year_survival = 1  alive with >= 60 months of follow-up,
#                           OR died after >= 60 months
#   five_year_survival = 0  died within 60 months
#   EXCLUDED                alive but followed for < 60 months - the 5-year
#                           outcome is unknown, so it cannot be labelled
#   The labels are taken directly from Month 1's feature_table_ml.csv.
#
# COHORT: patients with methylation data (671) AND a confirmed 5-year
#   outcome in Month 1 = 211 patients (152 survived / 59 died).
#
# SPLIT: stratified 70 / 30 (seed 42), done FIRST - before any scaling,
#   filling of missing values or fitting of methylation PCs - so nothing
#   about the test patients can leak into training.
#
# FEATURE BLOCKS (script 08 compares models built on each)
#   A  Clinical + genomic (36)  - Month 1's features, unchanged
#   B  Methylation summary (19) - global methylation (8), core-promoter
#                                 methylation of 10 genes, BRCA1 flag
#   C  Methylation patterns (10)- PCs of the 5,000 most variable CpGs,
#                                 FITTED ON TRAINING PATIENTS ONLY and then
#                                 applied to the test patients
#   Deliberately NOT used: survival columns (they contain the answer),
#   script 05's PCs and script 03's clusters (both computed with the test
#   patients included), PAM50 / TCGA clusters (from other data types),
#   plate / QC columns (technical, not biology).
#
# PREPROCESSING - same rules as revised Month 1 script 08
#   Continuous features: z-scored with the TRAINING mean and SD; missing
#   values then set to 0 (= the training mean). Binary features unscaled.
#   Features that are constant in the training set are dropped (they carry
#   no information and break the SVM's internal scaling).
#
# INPUT : ../01R_clinical_genomics/data/processed/TCGA-BRCA/feature_table_ml.csv
#         data/processed/methylation_patient_features.csv     (script 05)
#         data/processed/clinical_genomic_epigenomic_merged.csv (script 05)
#         data/processed/methylation_beta_matrix_patient.rds  (script 04)
# OUTPUT: data/processed/feature_table_epigenomic_ml.csv      (model-ready)
#         data/processed/feature_table_epigenomic_ml_raw.csv  (unscaled)
#         outputs/models/train_test_split_ids.csv
#         outputs/models/feature_sets.csv      (which features form A / B / C)
#         outputs/models/feature_preprocessing.rds (scaler + methylation PCA,
#                                                   to process new patients)
#         outputs/figures/07_features/ - 7 figures; outputs/tables/07_features/
#         outputs/reports/07_feature_table_report.txt
#
# MEMORY: ~2.5 GB (loads the 2 GB patient-level methylation matrix).
#
# Run from the Ubuntu terminal with:
#   Rscript scripts/07_build_epigenomic_feature_table.R
# Test run whose outputs get a suffix (nothing real is overwritten):
#   Rscript scripts/07_build_epigenomic_feature_table.R --tag=test
# ============================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)  # methylation container
  library(matrixStats)           # fast statistics without copying the matrix
  library(data.table)            # tables and CSV writing
  library(ggplot2)               # plots
  library(pheatmap)              # correlation heatmap
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
MONTH1_ML    <- "../01R_clinical_genomics/data/processed/TCGA-BRCA/feature_table_ml.csv"
METH_FEAT    <- "data/processed/methylation_patient_features.csv"
MERGED       <- "data/processed/clinical_genomic_epigenomic_merged.csv"
METH_MATRIX  <- "data/processed/methylation_beta_matrix_patient.rds"

# ---- Outputs ----
FIG_DIR  <- paste0("outputs/figures/07_features", sfx)
TAB_DIR  <- paste0("outputs/tables/07_features", sfx)
MOD_DIR  <- paste0("outputs/models", sfx)
REPORT   <- paste0("outputs/reports/07_feature_table_report", sfx, ".txt")
OUT_ML   <- paste0("data/processed/feature_table_epigenomic_ml", sfx, ".csv")
OUT_RAW  <- paste0("data/processed/feature_table_epigenomic_ml_raw", sfx, ".csv")
for (d in c(FIG_DIR, TAB_DIR, MOD_DIR, dirname(REPORT))) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Settings ----
TRAIN_FRAC <- 0.70       # 70 / 30 split (revised Month 1 used 80 / 20 on 305 patients)
SEED       <- 42         # same seed as Month 1
TOP_VAR_N  <- 5000       # most variable CpGs used for the methylation PCs
N_PC       <- 10         # methylation PCs kept as features

# Month 1's continuous columns (revised script 08, CONT_COLS) - scaled; the
# rest of block A (gene flags, histology one-hot columns) stays 0/1.
M1_CONT <- c("age", "stage_numeric", "tmb_nonsynonymous", "n_truncating", "n_nonsense",
             "n_frameshift_del", "n_frameshift_ins", "n_splice_site")
# Columns that must never be features (they contain the outcome or IDs).
FORBIDDEN <- c("os_event", "os_months", "OS_MONTHS", "OS_STATUS", "VITAL_STATUS", "is_alive",
               "survival_days", "survival_years", "DAYS_TO_DEATH", "DAYS_LAST_FOLLOWUP",
               "dfs_event", "DFS_MONTHS", "DFS_STATUS", "five_year_survival", "patient_id")

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
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-4, formatC(p, format = "e", digits = 1),
                                                   formatC(p, format = "f", digits = 4)))

PAL_OUT   <- c("Survived 5 years (1)" = "#185FA5", "Died within 5 years (0)" = "#EF4444")
PAL_BLOCK <- c("A: Clinical + genomic" = "#0F6E56", "B: Methylation summary" = "#185FA5",
               "C: Methylation PCs" = "#8C510A")
out_label <- function(y) factor(ifelse(y == 1, "Survived 5 years (1)", "Died within 5 years (0)"), levels = names(PAL_OUT))

cat("Script 07 - AI-ready epigenomic feature table\n")
cat("Tag    :", ifelse(nzchar(TAG), TAG, "(none - real run)"), "\n")
cat("Started:", format(run_start), "\n")

# ============================================================
# 1. LOAD INPUTS
# ============================================================
step("1. Load inputs")
for (f in c(MONTH1_ML, METH_FEAT, MERGED, METH_MATRIX))
  if (!file.exists(f)) stop("Input not found: ", f)

m1   <- fread(MONTH1_ML)                    # Month 1: 36 features + target, confirmed outcomes only
mf   <- fread(METH_FEAT)                    # script 05: methylation features (671 patients)
merg <- fread(MERGED, select = c("patient_id", "os_event", "OS_MONTHS"))
cat("Month 1 ML table (revised)  :", nrow(m1), "patients x", ncol(m1) - 2, "features\n")
cat("Methylation features (05)   :", nrow(mf), "patients\n")
stopifnot(!any(duplicated(m1$patient_id)), !any(duplicated(mf$patient_id)))

# ============================================================
# 2. COHORT: METHYLATION + CONFIRMED 5-YEAR OUTCOME
# ============================================================
step("2. Cohort and target")

# Why each of the 671 methylation patients is in or out of the ML cohort.
merg[, outcome_group := fifelse(is.na(OS_MONTHS) | is.na(os_event), "Missing follow-up (excluded)",
                        fifelse(os_event == 0 & OS_MONTHS >= 60, "Survived: alive, followed >= 60 months",
                        fifelse(os_event == 1 & OS_MONTHS >= 60, "Survived: died after 60 months",
                        fifelse(os_event == 1 & OS_MONTHS <  60, "Died within 60 months",
                                "Censored: alive, followed < 60 months (excluded)"))))]
flow_levels <- c("Survived: alive, followed >= 60 months", "Survived: died after 60 months",
                 "Died within 60 months", "Censored: alive, followed < 60 months (excluded)",
                 "Missing follow-up (excluded)")
flow <- merg[, .(patients = .N), by = outcome_group][match(flow_levels, outcome_group)][!is.na(outcome_group)]
flow[, pct := round(100 * patients / nrow(merg), 1)]
flow[, in_ml_cohort := !grepl("excluded", outcome_group)]
cat("The", nrow(merg), "methylation patients by 5-year outcome status:\n"); show(flow)
save_table(flow, "07_01_cohort_flow.csv")

# The cohort = methylation patients present in Month 1's confirmed-outcome table.
cohort <- merge(m1, mf[, .(patient_id)], by = "patient_id")
confirmed_ids <- merg[!grepl("excluded", outcome_group), patient_id]
stopifnot(setequal(cohort$patient_id, confirmed_ids))    # Month 1 labels = rule applied to our data
cat("\nML cohort:", nrow(cohort), "patients (checked: identical to applying the Month 1 rule to our data)\n")
cat("Month 1 confirmed-outcome patients without methylation data:", nrow(m1) - nrow(cohort), "of", nrow(m1), "\n")
tt <- cohort[, .N, by = five_year_survival][order(-five_year_survival)]
tt[, pct := round(100 * N / sum(N), 1)]
cat("\nTarget five_year_survival in the ML cohort:\n"); show(tt)

# Fig 1: where the 671 methylation patients go.
fl <- copy(flow); fl[, outcome_group := factor(outcome_group, levels = rev(flow_levels))]
p <- ggplot(fl, aes(patients, outcome_group, fill = in_ml_cohort)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%d  (%.1f%%)", patients, pct)), hjust = -0.08, size = 3.8) +
  scale_fill_manual(values = c(`TRUE` = "#185FA5", `FALSE` = "#B0B7BF"),
                    labels = c(`TRUE` = "In the ML cohort", `FALSE` = "Excluded"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(title = "Which methylation patients have a confirmed 5-year outcome?",
       subtitle = sprintf("Revised Month 1 target: %d of %d patients are usable for the 5-year survival models",
                          nrow(cohort), nrow(merg)),
       x = "Patients", y = NULL) + theme_qc
save_plot(p, "07_01_cohort_flow.png", 10.5, 4.8)

# ============================================================
# 3. STRATIFIED 70 / 30 SPLIT (done FIRST)
# ============================================================
step("3. Train / test split")
# Same method as Month 1: sample the training share separately within each
# outcome class, so both sets keep the same survived / died balance.
surv_ids <- sort(cohort[five_year_survival == 1, patient_id])
died_ids <- sort(cohort[five_year_survival == 0, patient_id])
train_ids <- c(sample(surv_ids, floor(TRAIN_FRAC * length(surv_ids))),
               sample(died_ids, floor(TRAIN_FRAC * length(died_ids))))
cohort[, split := fifelse(patient_id %in% train_ids, "train", "test")]
stopifnot(!any(cohort[split == "train", patient_id] %in% cohort[split == "test", patient_id]))

split_tab <- cohort[, .(patients = .N, survived = sum(five_year_survival == 1), died = sum(five_year_survival == 0),
                        pct_died = round(100 * mean(five_year_survival == 0), 1)), by = split][order(-patients)]
show(split_tab)
save_table(split_tab, "07_02_split_summary.csv")
split_ids <- cohort[, .(patient_id, split, five_year_survival)]
save_table(split_ids, "train_test_split_ids.csv", dir = MOD_DIR)

# Fig 2: the outcome balance in each set.
sb <- melt(split_tab, id.vars = "split", measure.vars = c("survived", "died"), variable.name = "outcome", value.name = "n")
sb[, outcome := factor(ifelse(outcome == "survived", "Survived 5 years (1)", "Died within 5 years (0)"), levels = names(PAL_OUT))]
sb[, split := factor(ifelse(split == "train", "Training (70%)", "Test (30%)"), levels = c("Training (70%)", "Test (30%)"))]
p <- ggplot(sb, aes(split, n, fill = outcome)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), colour = "white", fontface = "bold", size = 4.2) +
  scale_fill_manual(values = PAL_OUT, name = NULL) +
  labs(title = "Stratified 70 / 30 train-test split",
       subtitle = sprintf("Both sets keep ~%.0f%% deaths; seed %d; the test set is never used for fitting anything",
                          100 * mean(cohort$five_year_survival == 0), SEED),
       x = NULL, y = "Patients") + theme_qc
save_plot(p, "07_02_split_balance.png", 7.5, 5.5)

# ============================================================
# 4. BLOCK A - CLINICAL + GENOMIC (Month 1 features)
# ============================================================
step("4. Block A: clinical + genomic")
A_cols <- setdiff(names(m1), c("patient_id", "five_year_survival"))
cat("Block A features from Month 1:", length(A_cols), "\n")
na_A <- sapply(cohort[, ..A_cols], function(z) sum(is.na(z)))
if (any(na_A > 0)) { cat("Missing values in block A (filled later with the training mean):\n"); print(na_A[na_A > 0]) }

# ============================================================
# 5. BLOCK B - METHYLATION SUMMARY FEATURES (script 05)
# ============================================================
step("5. Block B: methylation summary")
B_cols <- c(grep("^meth_global_", names(mf), value = TRUE),
            grep("^meth_prom_", names(mf), value = TRUE),
            "meth_BRCA1_hypermethylated")
cat("Block B features:", length(B_cols), "\n  ", paste(B_cols, collapse = ", "), "\n")
cat("Not used from script 05 (see header): meth_PC1-10, meth_cluster, pam50_subtype, tcga_meth_cluster, plate/QC columns\n")
cohort <- merge(cohort, mf[, c("patient_id", B_cols), with = FALSE], by = "patient_id")

# ============================================================
# 6. BLOCK C - METHYLATION PCs FITTED ON TRAINING PATIENTS ONLY
# ============================================================
step("6. Block C: methylation PCs (training only)")
se <- readRDS(METH_MATRIX)
B  <- assay(se, withDimnames = FALSE)            # beta matrix, no copy
col_of <- setNames(seq_len(ncol(se)), colnames(se))
tr_col <- col_of[cohort[split == "train", patient_id]]
te_col <- col_of[cohort[split == "test",  patient_id]]
stopifnot(!anyNA(tr_col), !anyNA(te_col))

# 1) Pick the most variable CpGs using the TRAINING patients only.
sd_train <- rowSds(B, cols = tr_col, na.rm = TRUE)
top      <- order(sd_train, decreasing = TRUE)[seq_len(TOP_VAR_N)]
X_tr <- B[top, tr_col, drop = FALSE]
X_te <- B[top, te_col, drop = FALSE]
top_probes <- rowData(se)$probe_id[top]
rm(se, B); invisible(gc())

# 2) Fill gaps with the TRAINING mean of that CpG (applied to both sets).
probe_mean_tr <- rowMeans(X_tr, na.rm = TRUE)
fill <- function(X) { i <- which(is.na(X), arr.ind = TRUE); if (nrow(i)) X[i] <- probe_mean_tr[i[, 1]]; X }
na_tr <- sum(is.na(X_tr)); na_te <- sum(is.na(X_te))
X_tr <- fill(X_tr); X_te <- fill(X_te)
cat(sprintf("CpGs used: %d (chosen on training data) | gaps filled: %d train, %d test\n", TOP_VAR_N, na_tr, na_te))

# 3) Fit the PCA on training patients; project the test patients onto it.
pca <- prcomp(t(X_tr), center = TRUE, scale. = FALSE)
pc_tr <- pca$x[, seq_len(N_PC)]
pc_te <- predict(pca, newdata = t(X_te))[, seq_len(N_PC)]
C_cols <- paste0("methPC", seq_len(N_PC))
colnames(pc_tr) <- colnames(pc_te) <- C_cols
pcs <- rbind(data.table(patient_id = cohort[split == "train", patient_id], pc_tr),
             data.table(patient_id = cohort[split == "test",  patient_id], pc_te))
cohort <- merge(cohort, pcs, by = "patient_id")
var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pc_var <- data.table(PC = paste0("methPC", seq_along(var_exp))[1:20], pct_variance = round(var_exp[1:20], 2),
                     cumulative_pct = round(cumsum(var_exp)[1:20], 2))
cat(sprintf("Training PCA: PC1 %.1f%%, PC2 %.1f%%; the %d kept PCs explain %.1f%% of training variance\n",
            var_exp[1], var_exp[2], N_PC, sum(var_exp[seq_len(N_PC)])))
save_table(pc_var, "07_05_methylation_pca_variance_train.csv")

# Fig 4: training patients (dots) and projected test patients (triangles).
pv <- data.table(PC1 = cohort$methPC1, PC2 = cohort$methPC2,
                 set = factor(ifelse(cohort$split == "train", "Training", "Test (projected)"), levels = c("Training", "Test (projected)")),
                 outcome = out_label(cohort$five_year_survival))
p <- ggplot(pv, aes(PC1, PC2, colour = outcome, shape = set)) +
  geom_point(size = 2.3, alpha = 0.8, stroke = 1) +
  scale_colour_manual(values = PAL_OUT, name = NULL) +
  scale_shape_manual(values = c("Training" = 16, "Test (projected)" = 2), name = NULL) +
  guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  labs(title = "Methylation PCs fitted on training patients only",
       subtitle = "Test patients were projected onto the training PCA afterwards (no leakage)",
       x = sprintf("methPC1 (%.1f%% of training variance)", var_exp[1]),
       y = sprintf("methPC2 (%.1f%%)", var_exp[2])) + theme_qc + theme(legend.box = "vertical")
save_plot(p, "07_04_methylation_pca_train_test.png", 8.5, 7)

# Fig 5: scree plot of the training PCA.
sc <- pc_var[1:15]; sc[, PC := factor(PC, levels = PC)]
sc[, kept := ifelse(seq_len(.N) <= N_PC, "Kept as feature", "Not kept")]
p <- ggplot(sc, aes(PC, pct_variance, fill = kept)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", pct_variance)), vjust = -0.4, size = 3.1) +
  scale_fill_manual(values = c("Kept as feature" = "#8C510A", "Not kept" = "grey75"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Variance explained by the training-set methylation PCs",
       subtitle = sprintf("%s most variable CpGs; first %d PCs kept (%.1f%% of variance)",
                          format(TOP_VAR_N, big.mark = ","), N_PC, sum(var_exp[seq_len(N_PC)])),
       x = NULL, y = "Variance explained (%)") +
  theme_qc + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p, "07_05_methylation_pca_scree_train.png", 9, 5)

# ============================================================
# 7. ASSEMBLE, DROP CONSTANT FEATURES, SCALE WITH TRAINING STATISTICS
# ============================================================
step("7. Assemble and preprocess")
feat_cols <- c(A_cols, B_cols, C_cols)
bad <- intersect(feat_cols, FORBIDDEN)
if (length(bad)) stop("Outcome/ID columns found among features: ", paste(bad, collapse = ", "))
block_of <- c(setNames(rep(names(PAL_BLOCK)[1], length(A_cols)), A_cols),
              setNames(rep(names(PAL_BLOCK)[2], length(B_cols)), B_cols),
              setNames(rep(names(PAL_BLOCK)[3], length(C_cols)), C_cols))

# Raw (unscaled) table - written before any preprocessing.
setcolorder(cohort, c("patient_id", "split", "five_year_survival", feat_cols))
raw_tab <- cohort[, c("patient_id", "split", "five_year_survival", feat_cols), with = FALSE]
fwrite(raw_tab, OUT_RAW); record(OUT_RAW, "processed data")
cat("  saved raw feature table:", OUT_RAW, "(", nrow(raw_tab), "x", length(feat_cols), "features )\n")

# Drop features that do not vary in the TRAINING set.
tr <- cohort[split == "train"]
const <- feat_cols[sapply(feat_cols, function(v) { x <- tr[[v]]; length(unique(x[!is.na(x)])) <= 1 })]
if (length(const)) {
  cat("Dropped - constant in the training set:", paste(const, collapse = ", "), "\n")
  save_table(data.table(feature = const, block = block_of[const],
                        reason = "constant (no variation) in the training set"), "07_07_dropped_constant_features.csv")
}
feat_cols <- setdiff(feat_cols, const)

# Which features are continuous (scaled) vs binary (left as 0/1)?
is_binary <- sapply(feat_cols, function(v) all(cohort[[v]][!is.na(cohort[[v]])] %in% c(0, 1)))
cont_cols <- feat_cols[!is_binary | feat_cols %in% M1_CONT]
bin_cols  <- setdiff(feat_cols, cont_cols)
cat(sprintf("Features kept: %d | continuous (scaled): %d | binary (unscaled): %d\n",
            length(feat_cols), length(cont_cols), length(bin_cols)))

# Training mean / SD, applied to BOTH sets; then missing -> 0 (= training mean).
tr_mean <- sapply(cont_cols, function(v) mean(tr[[v]], na.rm = TRUE))
tr_sd   <- sapply(cont_cols, function(v) sd(tr[[v]], na.rm = TRUE))
ml <- copy(cohort[, c("patient_id", "split", "five_year_survival", feat_cols), with = FALSE])
n_filled <- integer(0)
for (v in cont_cols) {
  set(ml, j = v, value = (ml[[v]] - tr_mean[[v]]) / tr_sd[[v]])
  k <- is.na(ml[[v]]); if (any(k)) { n_filled[v] <- sum(k); set(ml, i = which(k), j = v, value = 0) }
}
for (v in bin_cols) {                                    # binary: fill any gap with the training majority value
  k <- is.na(ml[[v]]); if (any(k)) { n_filled[v] <- sum(k); set(ml, i = which(k), j = v, value = as.numeric(mean(tr[[v]], na.rm = TRUE) >= 0.5)) }
}
if (length(n_filled)) { cat("Missing values filled (training mean for continuous features):\n"); print(n_filled) }

stopifnot(!anyNA(ml), nrow(ml) == nrow(cohort))
fwrite(ml, OUT_ML); record(OUT_ML, "processed data")
cat("  saved model-ready feature table:", OUT_ML, "(", nrow(ml), "x", length(feat_cols), "features )\n")

# Everything needed to process a NEW patient in exactly the same way.
prep <- list(cont_cols = cont_cols, bin_cols = bin_cols, train_mean = tr_mean, train_sd = tr_sd,
             dropped_constant = const, meth_top_probes = top_probes, meth_probe_train_mean = probe_mean_tr,
             meth_pca = list(rotation = pca$rotation[, seq_len(N_PC)], center = pca$center),
             seed = SEED, train_frac = TRAIN_FRAC, created = format(Sys.time()))
saveRDS(prep, file.path(MOD_DIR, "feature_preprocessing.rds")); record(file.path(MOD_DIR, "feature_preprocessing.rds"), "model object")
cat("  saved preprocessing parameters:", file.path(MOD_DIR, "feature_preprocessing.rds"), "\n")

# Feature sets for script 08: A = Month 1 baseline, B+C = methylation only, A+B+C = combined.
fsets <- data.table(feature = feat_cols, block = block_of[feat_cols],
                    type = ifelse(feat_cols %in% cont_cols, "continuous (scaled)", "binary (0/1)"))
fsets[, set_A_clinical_genomic := block == names(PAL_BLOCK)[1]]
fsets[, set_BC_methylation     := block != names(PAL_BLOCK)[1]]
fsets[, set_ABC_combined       := TRUE]
save_table(fsets, "feature_sets.csv", dir = MOD_DIR)
cat("\nFeature sets for script 08:\n")
show(fsets[, .(features = .N), by = block])
cat(sprintf("  Set A (clinical + genomic): %d | Set B+C (methylation): %d | Set A+B+C (combined): %d\n",
            sum(fsets$set_A_clinical_genomic), sum(fsets$set_BC_methylation), nrow(fsets)))

# Feature dictionary (training-set summaries only).
describe <- function(v) {
  if (v %in% M1_CONT) return("Month 1 clinical / genomic measure")
  if (v %in% A_cols && grepl("^diag_", v)) return(paste("Histology one-hot:", sub("^diag_", "", v)))
  if (v %in% A_cols) return(paste(v, "gene mutated (1) or not (0) - Month 1 gene flag"))
  if (grepl("^meth_global_", v)) return(paste("Mean methylation,", sub("^meth_global_", "", v), "CpGs"))
  if (grepl("^meth_prom_", v)) return(paste("Core-promoter methylation of", sub("^meth_prom_", "", v)))
  if (v == "meth_BRCA1_hypermethylated") return("BRCA1 promoter silenced by methylation (1/0)")
  if (grepl("^methPC", v)) return(paste("Methylation principal component", sub("^methPC", "", v), "(training-set PCA)"))
  ""
}
dict <- data.table(feature = feat_cols, block = block_of[feat_cols],
                   type = ifelse(feat_cols %in% cont_cols, "continuous (scaled)", "binary (0/1)"),
                   description = vapply(feat_cols, describe, character(1)),
                   train_mean = round(sapply(feat_cols, function(v) mean(tr[[v]], na.rm = TRUE)), 4),
                   train_sd = round(sapply(feat_cols, function(v) sd(tr[[v]], na.rm = TRUE)), 4),
                   n_missing = sapply(feat_cols, function(v) sum(is.na(cohort[[v]]))))
save_table(dict, "07_03_feature_dictionary.csv")
save_table(data.table(feature = cont_cols, train_mean = round(tr_mean, 6), train_sd = round(tr_sd, 6)),
           "07_04_preprocessing_parameters.csv")

# ============================================================
# 8. EXPLORATION ON THE TRAINING SET ONLY
#    Nothing below changes the feature table - it is for understanding
#    only, and uses training patients so the test set stays untouched.
# ============================================================
step("8. Training-set exploration")
trs <- ml[split == "train"]                               # scaled training data
y   <- trs$five_year_survival

# Univariate screen: how well does each feature ALONE separate the classes?
# AUC = probability that a random survivor has a higher value than a random
# non-survivor (0.5 = no signal; far from 0.5 in either direction = signal).
uni_auc <- function(x, y) { r <- rank(x); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
screen <- rbindlist(lapply(feat_cols, function(v) {
  x <- trs[[v]]
  data.table(feature = v, block = block_of[v], auc = uni_auc(x, y),
             p_value = tryCatch(suppressWarnings(wilcox.test(x[y == 1], x[y == 0])$p.value), error = function(e) NA_real_))
}))
screen[, strength := abs(auc - 0.5)]
screen[, auc_either_direction := pmax(auc, 1 - auc)]      # e.g. AUC 0.32 = 0.68 read the other way round
screen[, direction := ifelse(auc >= 0.5, "higher in survivors", "higher in non-survivors")]
screen[, q_value := p.adjust(p_value, method = "BH")]
setorder(screen, -strength)
cat("Top 15 single features (training set, exploratory - not used to select features):\n")
show(head(screen[, .(feature, block, auc = round(auc_either_direction, 3), direction,
                     p_value = signif(p_value, 3), q_value = signif(q_value, 3))], 15))
cat("(auc shown in its stronger direction: 0.5 = no signal, 1 = perfect separation)\n")
save_table(screen, "07_06_univariate_feature_screen_train.csv")

# Fig 7: univariate signal of the top 25 features, coloured by block.
top25 <- head(screen, 25)
top25[, feature := factor(feature, levels = rev(feature))]
p <- ggplot(top25, aes(auc - 0.5, feature, fill = block)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, colour = "grey40") +
  scale_fill_manual(values = PAL_BLOCK, name = NULL) +
  scale_x_continuous(labels = function(z) sprintf("%.2f", z + 0.5)) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(title = "How well does each feature alone separate the outcomes?",
       subtitle = "Training set only; single-feature AUC (0.5 = no signal); right = higher in survivors, left = higher in non-survivors",
       x = "Single-feature AUC", y = NULL) + theme_qc
save_plot(p, "07_07_univariate_feature_screen.png", 10, 8)

# Fig 3: block B methylation features by outcome (unscaled values, training set).
b_cont <- intersect(setdiff(B_cols, "meth_BRCA1_hypermethylated"), feat_cols)
bl <- rbindlist(lapply(b_cont, function(v) {
  pv <- screen[feature == v, p_value]
  data.table(feature = sprintf("%s (p = %s)", sub("^meth_", "", v), fmt_p(pv)),
             value = cohort[split == "train"][[v]], outcome = out_label(cohort[split == "train", five_year_survival]))
}))
bl[, feature := factor(feature, levels = unique(feature))]
p <- ggplot(bl, aes(outcome, value, fill = outcome)) +
  geom_boxplot(outlier.size = 0.6, width = 0.6, alpha = 0.85) +
  facet_wrap(~ feature, ncol = 6, scales = "free_y") +
  scale_fill_manual(values = PAL_OUT, name = NULL) +
  labs(title = "Methylation summary features by 5-year outcome (training set)",
       subtitle = "Wilcoxon p-values in panel titles (not corrected for multiple testing); each panel has its own y-axis",
       x = NULL, y = "Beta (methylation level)") +
  theme_qc + theme(axis.text.x = element_blank(), strip.text = element_text(size = 8))
save_plot(p, "07_03_methylation_features_by_outcome.png", 13, 8)

# Fig 6: correlation between all features (training set) - shows redundancy.
cm <- suppressWarnings(cor(as.matrix(trs[, ..feat_cols]), method = "spearman", use = "pairwise.complete.obs"))
cm[is.na(cm)] <- 0
ann <- data.frame(Block = factor(block_of[feat_cols], levels = names(PAL_BLOCK)), row.names = feat_cols)
pheatmap(cm, annotation_row = ann, annotation_col = ann, annotation_colors = list(Block = PAL_BLOCK),
         color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100), breaks = seq(-1, 1, length.out = 101),
         border_color = NA, fontsize = 9, fontsize_row = 6.5, fontsize_col = 6.5, treeheight_row = 20, treeheight_col = 20,
         annotation_names_row = FALSE, annotation_names_col = FALSE,
         main = "Correlation between features (training set, Spearman)",
         filename = file.path(FIG_DIR, "07_06_feature_correlation_heatmap.png"), width = 12, height = 11.5)
record(file.path(FIG_DIR, "07_06_feature_correlation_heatmap.png"), "figure")
cat("  saved figure:", file.path(FIG_DIR, "07_06_feature_correlation_heatmap.png"), "\n")
hi <- which(abs(cm) >= 0.9 & upper.tri(cm), arr.ind = TRUE)
if (nrow(hi)) {
  cat(sprintf("Feature pairs with |correlation| >= 0.9: %d (e.g. %s)\n", nrow(hi),
              paste(head(sprintf("%s ~ %s", rownames(cm)[hi[, 1]], colnames(cm)[hi[, 2]]), 3), collapse = "; ")))
}

# ============================================================
# 9. KEY POINTS AND OUTPUT INDEX
# ============================================================
step("9. Summary")
top_meth <- screen[block != names(PAL_BLOCK)[1]][1]
top_all  <- screen[1]
cat("\n====================== SUMMARY ======================\n")
cat("* Target: five_year_survival, revised Month 1 definition (confirmed 5-year outcomes only).\n")
cat(sprintf("* Cohort: %d patients (%d survived / %d died, %.1f%% died); %d methylation patients excluded as censored or missing.\n",
            nrow(cohort), sum(cohort$five_year_survival == 1), sum(cohort$five_year_survival == 0),
            100 * mean(cohort$five_year_survival == 0), nrow(merg) - nrow(cohort)))
cat(sprintf("* Split: %d training / %d test patients (%d / %d deaths), stratified, seed %d.\n",
            sum(cohort$split == "train"), sum(cohort$split == "test"),
            sum(cohort$split == "train" & cohort$five_year_survival == 0), sum(cohort$split == "test" & cohort$five_year_survival == 0), SEED))
cat(sprintf("* Features: %d in total - A clinical+genomic %d, B methylation summary %d, C methylation PCs %d%s.\n",
            length(feat_cols), sum(block_of[feat_cols] == names(PAL_BLOCK)[1]), sum(block_of[feat_cols] == names(PAL_BLOCK)[2]),
            sum(block_of[feat_cols] == names(PAL_BLOCK)[3]),
            if (length(const)) sprintf(" (%d constant features dropped)", length(const)) else ""))
cat("* Leakage checks passed: split before preprocessing; scaler, gap filling and PCA fitted on training only; no outcome columns among features.\n")
cat(sprintf("* Strongest single feature (training): %s (AUC %.2f, %s); strongest methylation feature: %s (AUC %.2f, %s).\n",
            top_all$feature, top_all$auc_either_direction, top_all$direction,
            top_meth$feature, top_meth$auc_either_direction, top_meth$direction))
cat(sprintf("* Only %d feature(s) pass q < 0.05 on their own - expected with 41 training deaths; script 08 tests them together.\n",
            sum(screen$q_value < 0.05, na.rm = TRUE)))
cat("* Next: script 08 trains LR, Elastic Net, RF, XGBoost and SVM on sets A, B+C and A+B+C.\n")
cat("=====================================================\n")

idx <- do.call(rbind, outputs_written)
idx$size_kb <- round(file.size(idx$file) / 1024, 1)
save_table(idx, "07_99_outputs_index.csv")
cat("\nFiles written:", nrow(idx), "(", sum(idx$type == "figure"), "figures,", sum(idx$type == "table"), "tables )\n")
cat("Total run time:", fmt_time(as.numeric(difftime(Sys.time(), run_start, units = "secs"))), "\n")
cat("Report saved to:", REPORT, "\n")
cat("\nScript 07 complete.\n")

sink()
close(report_con)
