#!/usr/bin/env Rscript
# ============================================================
# 04_match_harmonise_ids.R
#
# Harmonises the methylation sample IDs with the Month 1 cohort and
# collapses the data to ONE tumour sample per patient.
# (Internship brief: "12_match_methylation_ids".)
#
# WHAT THIS SCRIPT DOES, IN ORDER
#   1. Loads the QC-cleaned methylation object from script 03
#   2. Loads the Month 1 cohort (964 patients) and the tumour sample IDs
#      Month 1 used for its mutation data
#   3. ID checks - barcode format, patient/sample IDs consistent, primary
#      tumours only, every methylation patient present in Month 1
#   4. Replicates - 6 patients have 2-3 tumour samples; keeps exactly one per
#      patient using a fixed, documented rule and records WHY
#   5. Builds the patient-level object: columns renamed from barcodes to
#      patient IDs, so it can be merged directly with Month 1 in script 05
#   6. Cohort flow - explains what happened to every Month 1 patient
#      (matched / removed at QC / only non-tumour samples / no data on GDC)
#   7. Representativeness - are the patients WITH methylation similar to
#      those without? (age, stage, survival, mutations, PAM50)
#   8. Saves the patient-level object, an ID map and a master ID table
#
# HOW ONE SAMPLE IS CHOSEN FOR A REPLICATE PATIENT (in this order)
#   1. Same tumour vial as Month 1's mutation data (e.g. both "-01A"), so
#      methylation and mutations describe the same piece of tumour
#   2. QC status "Pass" before "Flagged" (from script 03)
#   3. Vial A before later vials (B, C ...)
#   4. Fewer sample-specific missing values
#   5. Alphabetical barcode (only as a final, reproducible tie-break)
#
# INPUT : data/processed/methylation_beta_matrix_cleaned.rds   (script 03)
#         ../01_clinical_genomics/data/processed/clinical_genomics_merged.csv
#         ../01_clinical_genomics/data/processed/clinical_sample_cleaned.csv
#         data/processed/methylation_matched_manifest.csv       (script 01)
#         data/processed/methylation_query_object.rds  (optional; adds detail
#                                                       to the cohort flow)
# OUTPUT: data/processed/methylation_beta_matrix_patient.rds
#           SummarizedExperiment, probes x patients (columns = patient IDs)
#         data/processed/methylation_patient_id_map.csv   (one row per patient)
#         data/processed/matched_patient_ids.csv          (all 964 Month 1 patients)
#         outputs/figures/04_ids/  - 5 figures (PNG, 300 dpi)
#         outputs/tables/04_ids/   - summary tables (CSV)
#         outputs/reports/04_id_matching_report.txt - everything printed
#
# MEMORY: loads the ~2 GB cleaned matrix and makes one ~2 GB copy when
# dropping the replicate columns, so peak is ~4-5 GB. Close web browsers.
#
# Run from the Ubuntu terminal with:
#   Rscript scripts/04_match_harmonise_ids.R
# Test run whose outputs get a suffix (nothing real is overwritten):
#   Rscript scripts/04_match_harmonise_ids.R --tag=test
# ============================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)  # the methylation container (matrix + annotation)
  library(data.table)            # fast tables and CSV writing
  library(ggplot2)               # all plots
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
INPUT <- getflag("input", "data/processed/methylation_beta_matrix_cleaned.rds")
TAG   <- getflag("tag", "")
sfx   <- if (nzchar(TAG)) paste0("_", TAG) else ""

# ---- Input files ----
MONTH1_DIR     <- "../01_clinical_genomics/data/processed"
MONTH1_MERGED  <- file.path(MONTH1_DIR, "clinical_genomics_merged.csv")   # same file script 01 used
MONTH1_SAMPLES <- file.path(MONTH1_DIR, "clinical_sample_cleaned.csv")    # tumour sample IDs
MONTH1_SAMPLES_ALT <- file.path(MONTH1_DIR, "TCGA-BRCA", "sample_cleaned.csv")
MANIFEST       <- "data/processed/methylation_matched_manifest.csv"
QUERY_CACHE    <- "data/processed/methylation_query_object.rds"

# ---- Output locations ----
FIG_DIR  <- paste0("outputs/figures/04_ids", sfx)
TAB_DIR  <- paste0("outputs/tables/04_ids", sfx)
REPORT   <- paste0("outputs/reports/04_id_matching_report", sfx, ".txt")
OUT_RDS  <- paste0("data/processed/methylation_beta_matrix_patient", sfx, ".rds")
ID_MAP   <- paste0("data/processed/methylation_patient_id_map", sfx, ".csv")
MASTER   <- paste0("data/processed/matched_patient_ids", sfx, ".csv")
for (d in c(FIG_DIR, TAB_DIR, dirname(REPORT), dirname(OUT_RDS)))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

SEED <- 42          # same seed as Month 1 (used by the simulated chi-square p-values)
set.seed(SEED)

# ---- Everything printed goes to the screen AND the report file ----
report_con <- file(REPORT, open = "wt")
sink(report_con, split = TRUE)

# ---- Timing / section headers ----
run_start <- Sys.time()
fmt_time  <- function(secs) if (secs < 60) sprintf("%.0f s", secs) else sprintf("%.1f min", secs / 60)
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

# ---- Consistent plot style (same as script 03) ----
theme_qc <- theme_minimal(base_size = 12) +
  theme(plot.title       = element_text(face = "bold", size = 14),
        plot.subtitle    = element_text(colour = "grey30", size = 10.5),
        legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        strip.text       = element_text(face = "bold"),
        plot.margin      = margin(10, 15, 10, 10))

save_plot <- function(p, name, width, height) {
  path <- file.path(FIG_DIR, name)
  ggsave(path, p, width = width, height = height, dpi = 300, bg = "white")
  record(path, "figure"); cat("  saved figure:", path, "\n")
}
save_table <- function(df, name, dir = TAB_DIR) {
  path <- file.path(dir, name)
  fwrite(df, path); record(path, "table"); cat("  saved table :", path, "\n")
}
show  <- function(df) print(as.data.frame(df), row.names = FALSE)
fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-4, formatC(p, format = "e", digits = 1),
                                                   formatC(p, format = "f", digits = 4)))
norm_id <- function(x) toupper(trimws(as.character(x)))   # guards against stray spaces / case

PAM50_LEVELS <- c("LumA", "LumB", "Her2", "Basal", "Normal", "Unknown")
PAL_PAM50  <- c(LumA = "#2C7BB6", LumB = "#74ADD1", Her2 = "#F46D43",
                Basal = "#D73027", Normal = "#1A9850", Unknown = "grey70")
PAL_GROUP  <- c("With methylation" = "#185FA5", "Without methylation" = "#B0B7BF")

cat("Script 04 - ID harmonisation and one sample per patient\n")
cat("Input  :", INPUT, "\n")
cat("Tag    :", ifelse(nzchar(TAG), TAG, "(none - real run)"), "\n")
cat("Started:", format(run_start), "\n")

# ============================================================
# 1. LOAD THE CLEANED METHYLATION OBJECT
# ============================================================
step("1. Load methylation data")
if (!file.exists(INPUT)) stop("Input not found: ", INPUT, "\nRun scripts/03_clean_methylation_data.R first.")
se <- readRDS(INPUT)                                   # probes x samples
cd <- as.data.frame(colData(se))                       # one row per sample
cd$barcode    <- norm_id(cd$barcode)
cd$patient_id <- norm_id(cd$patient_id)
cd$sample_id  <- norm_id(cd$sample_id)
cd$col_index  <- seq_len(nrow(cd))                     # position in the matrix

cat("Probes              :", fmt_n(nrow(se)), "\n")
cat("Samples (columns)   :", ncol(se), "\n")
cat("Patients            :", length(unique(cd$patient_id)), "\n")
cat("Patients with >1 sample:", sum(table(cd$patient_id) > 1), "\n")

# ============================================================
# 2. LOAD THE MONTH 1 COHORT AND ITS TUMOUR SAMPLE IDS
# ============================================================
step("2. Load Month 1 cohort")
if (!file.exists(MONTH1_MERGED)) stop("Month 1 merged dataset not found: ", MONTH1_MERGED)
m1 <- fread(MONTH1_MERGED)
m1$patient_id <- norm_id(m1$patient_id)
cat("Month 1 patients    :", nrow(m1), "(unique:", length(unique(m1$patient_id)), ")\n")

# The tumour sample Month 1's mutation data came from, e.g. TCGA-A7-A26J-01A.
# Used so that, for replicate patients, methylation comes from the same vial.
samp_file <- if (file.exists(MONTH1_SAMPLES)) MONTH1_SAMPLES else MONTH1_SAMPLES_ALT
if (file.exists(samp_file)) {
  m1s <- fread(samp_file)
  m1s$patient_id <- norm_id(m1s$patient_id)
  m1s$sample_id  <- norm_id(m1s$sample_id)
  if ("sample_type" %in% names(m1s)) m1s <- m1s[sample_type == "Primary Tumor"]
  month1_samples <- unique(m1s$sample_id)
  cat("Month 1 primary-tumour sample IDs:", length(month1_samples), "(from", samp_file, ")\n")
} else {
  month1_samples <- character(0)
  cat("Month 1 sample file not found - the 'same vial as Month 1' rule will be skipped.\n")
}

# ============================================================
# 3. ID CONSISTENCY CHECKS
#    Every check must pass; a failure here would mean samples could be
#    matched to the wrong patient later on, so the script stops.
# ============================================================
step("3. ID checks")
bc_pattern <- "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]-[0-9]{2}[A-Z]-[A-Z0-9]{4}-[0-9]{2}$"
checks <- data.frame(
  check = c("Barcode has the standard TCGA format",
            "patient_id = first 12 characters of the barcode",
            "sample_id = first 16 characters of the barcode",
            "Sample type code is 01 (primary tumour)",
            "No duplicated barcodes",
            "Month 1 patient IDs are unique",
            "Every methylation patient is in the Month 1 cohort"),
  n_pass  = c(sum(grepl(bc_pattern, cd$barcode)),
              sum(cd$patient_id == substr(cd$barcode, 1, 12)),
              sum(cd$sample_id  == substr(cd$barcode, 1, 16)),
              sum(substr(cd$barcode, 14, 15) == "01"),
              sum(!duplicated(cd$barcode)),
              sum(!duplicated(m1$patient_id)),
              sum(unique(cd$patient_id) %in% m1$patient_id)),
  n_total = c(rep(nrow(cd), 5), nrow(m1), length(unique(cd$patient_id))))
checks$result <- ifelse(checks$n_pass == checks$n_total, "PASS", "FAIL")
show(checks)
save_table(checks, "04_01_id_checks.csv")
if (any(checks$result == "FAIL"))
  stop("One or more ID checks failed (see table above). Fix the IDs before matching.")
cat("All ID checks passed.\n")

# ============================================================
# 4. RESOLVE REPLICATES: EXACTLY ONE SAMPLE PER PATIENT
# ============================================================
step("4. One sample per patient")
cd$vial_matches_month1 <- cd$sample_id %in% month1_samples
cd$vial_rank <- match(cd$vial, LETTERS)                 # A = 1, B = 2 ...

# Sort so the preferred sample comes first within each patient (rules 1-5).
ord <- order(cd$patient_id,
             -as.integer(cd$vial_matches_month1),       # 1. same vial as Month 1 first
             cd$qc_status != "Pass",                    # 2. Pass before Flagged
             cd$vial_rank,                              # 3. vial A before B
             cd$pct_missing,                            # 4. fewer missing values
             cd$barcode)                                # 5. alphabetical tie-break
cd_sorted <- cd[ord, ]
cd_sorted$rank_in_patient <- ave(seq_len(nrow(cd_sorted)), cd_sorted$patient_id, FUN = seq_along)
chosen <- cd_sorted[cd_sorted$rank_in_patient == 1, ]

# Explains, in words, which rule separated the kept sample from the runner-up.
why_chosen <- function(pt) {
  s <- cd_sorted[cd_sorted$patient_id == pt, ]
  if (nrow(s) == 1) return("only sample")
  a <- s[1, ]; b <- s[2, ]
  if (a$vial_matches_month1 != b$vial_matches_month1)
    return(sprintf("same tumour vial as Month 1 mutation data (%s)", substr(a$sample_id, 14, 16)))
  if ((a$qc_status == "Pass") != (b$qc_status == "Pass")) return("better QC status (Pass vs Flagged)")
  if (a$vial_rank != b$vial_rank) return(sprintf("vial %s preferred over vial %s", a$vial, b$vial))
  if (a$pct_missing != b$pct_missing)
    return(sprintf("fewer missing values (%.2f%% vs %.2f%%)", a$pct_missing, b$pct_missing))
  "alphabetical tie-break"
}
chosen$selection_reason <- vapply(chosen$patient_id, why_chosen, character(1))

# Detailed table for the replicate patients only, including how closely each
# dropped sample agrees with the kept one (Pearson r over all probes).
rep_pat <- unique(cd$patient_id[duplicated(cd$patient_id)])
rep_tab <- cd_sorted[cd_sorted$patient_id %in% rep_pat,
                     c("patient_id", "barcode", "vial", "plate", "qc_status", "pct_missing",
                       "cor_to_median", "vial_matches_month1", "rank_in_patient", "col_index")]
rep_tab$kept <- rep_tab$rank_in_patient == 1
# withDimnames = FALSE: plain assay(se) re-attaches row/column names and can
# silently copy the whole 2 GB matrix on EVERY call; this handle does not.
B <- assay(se, withDimnames = FALSE)
rep_tab$r_with_kept_sample <- vapply(seq_len(nrow(rep_tab)), function(k) {
  keep_col <- chosen$col_index[chosen$patient_id == rep_tab$patient_id[k]]
  if (rep_tab$col_index[k] == keep_col) return(1)
  a <- B[, keep_col]; b <- B[, rep_tab$col_index[k]]; ok <- !is.na(a) & !is.na(b)
  cor(a[ok], b[ok])
}, numeric(1))
rm(B)
rep_tab$decision <- ifelse(rep_tab$kept,
                           paste("KEPT -", chosen$selection_reason[match(rep_tab$patient_id, chosen$patient_id)]),
                           "dropped")
rep_out <- rep_tab[, c("patient_id", "barcode", "vial", "plate", "qc_status", "pct_missing",
                       "vial_matches_month1", "r_with_kept_sample", "decision")]
rep_out$pct_missing        <- round(rep_out$pct_missing, 3)
rep_out$r_with_kept_sample <- round(rep_out$r_with_kept_sample, 4)

cat("Replicate patients:", length(rep_pat), "| samples involved:", nrow(rep_tab),
    "| samples dropped:", sum(!rep_tab$kept), "\n\n")
show(rep_out)
save_table(rep_out, "04_02_replicate_resolution.csv")
cat("\nWhy each replicate patient's sample was chosen:\n")
print(table(chosen$selection_reason[chosen$patient_id %in% rep_pat]))

stopifnot(!any(duplicated(chosen$patient_id)))            # exactly one per patient
cat("\nSamples kept:", nrow(chosen), "=", length(unique(cd$patient_id)), "patients (one each)\n")

# ============================================================
# 5. BUILD THE PATIENT-LEVEL OBJECT (columns = patient IDs)
# ============================================================
step("5. Patient-level object")
keep_cols <- sort(chosen$col_index)                       # keep the original column order
cat("Dropping", ncol(se) - length(keep_cols), "replicate columns and renaming columns to patient IDs...\n")
se_p <- se[, keep_cols]                                   # the one ~2 GB copy
rm(se); invisible(gc())

ch <- chosen[match(keep_cols, chosen$col_index), ]
colData(se_p)$methylation_barcode <- colnames(se_p)       # keep the original barcode
colData(se_p)$selection_reason    <- ch$selection_reason
colData(se_p)$vial_matches_month1 <- ch$vial_matches_month1
colnames(se_p) <- ch$patient_id                           # columns now named by patient
stopifnot(!any(duplicated(colnames(se_p))),
          identical(colnames(se_p), norm_id(colData(se_p)$patient_id)))
S4Vectors::metadata(se_p)$script04 <- list(
  description = "One primary-tumour sample per patient; columns named by patient_id",
  selection_rules = c("same vial as Month 1 mutation data", "QC Pass before Flagged",
                      "vial A before later vials", "fewer missing values", "alphabetical barcode"),
  created = format(Sys.time()))
cat("Patient-level object:", fmt_n(nrow(se_p)), "probes x", ncol(se_p), "patients\n")

id_map <- data.frame(
  patient_id          = colnames(se_p),
  methylation_barcode = colData(se_p)$methylation_barcode,
  sample_id           = colData(se_p)$sample_id,
  plate               = colData(se_p)$plate,
  qc_status           = colData(se_p)$qc_status,
  pam50               = colData(se_p)$pam50,
  n_samples_available = colData(se_p)$n_samples_for_patient,
  vial_matches_month1 = colData(se_p)$vial_matches_month1,
  selection_reason    = colData(se_p)$selection_reason)
save_table(id_map, basename(ID_MAP), dir = dirname(ID_MAP))

# ============================================================
# 6. COHORT FLOW: WHAT HAPPENED TO EVERY MONTH 1 PATIENT?
# ============================================================
step("6. Cohort flow")
man <- fread(MANIFEST)
man_patients <- unique(norm_id(man$patient_id))           # had a primary-tumour file (script 01)

gdc_any <- character(0)                                   # had ANY 450K file on the GDC
if (file.exists(QUERY_CACHE)) {
  q   <- readRDS(QUERY_CACHE)
  res <- as.data.frame(q$results[[1]])
  pid <- if ("cases.submitter_id" %in% names(res)) res$cases.submitter_id else substr(res$cases, 1, 12)
  gdc_any <- unique(norm_id(pid))
  cat("GDC query cache found:", length(gdc_any), "patients with any TCGA-BRCA 450K file\n")
} else {
  cat("GDC query cache not found - patients without methylation are grouped together.\n")
}

final_patients <- colnames(se_p)
status <- ifelse(m1$patient_id %in% final_patients, "Matched: methylation available",
          ifelse(m1$patient_id %in% man_patients, "Removed at methylation QC (script 03)",
          ifelse(m1$patient_id %in% gdc_any, "Methylation only from non-primary-tumour samples",
          ifelse(length(gdc_any) > 0, "No 450K methylation data on the GDC",
                 "No usable methylation sample"))))
status_levels <- c("Matched: methylation available", "Removed at methylation QC (script 03)",
                   "Methylation only from non-primary-tumour samples", "No 450K methylation data on the GDC",
                   "No usable methylation sample")

master <- data.frame(
  patient_id          = m1$patient_id,
  has_methylation     = m1$patient_id %in% final_patients,
  methylation_status  = status,
  methylation_barcode = id_map$methylation_barcode[match(m1$patient_id, id_map$patient_id)],
  month1_sample_id    = if (length(month1_samples)) {
    s <- m1s[, .(sample_id = sample_id[1]), by = patient_id]
    s$sample_id[match(m1$patient_id, s$patient_id)]
  } else NA_character_)
save_table(master, basename(MASTER), dir = dirname(MASTER))

flow <- as.data.frame(table(status = factor(status, levels = status_levels)))
flow <- flow[flow$Freq > 0, ]; names(flow)[2] <- "patients"
flow$pct_of_month1 <- round(100 * flow$patients / nrow(m1), 1)
cat("\nMonth 1 cohort (", nrow(m1), "patients ) by methylation status:\n"); show(flow)
save_table(flow, "04_03_cohort_flow.csv")
if (any(status == "Removed at methylation QC (script 03)")) {
  cat("\nRemoved at QC:", paste(master$patient_id[status == "Removed at methylation QC (script 03)"], collapse = ", "), "\n")
}

# Fig 1: where the 964 Month 1 patients ended up.
fl <- flow; fl$status <- factor(fl$status, levels = rev(status_levels))
p <- ggplot(fl, aes(patients, status)) +
  geom_col(aes(fill = status == "Matched: methylation available"), width = 0.65) +
  geom_text(aes(label = sprintf("%d  (%.1f%%)", patients, pct_of_month1)), hjust = -0.08, size = 3.8) +
  scale_fill_manual(values = c(`TRUE` = "#185FA5", `FALSE` = "#B0B7BF"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(title = "What happened to each Month 1 patient",
       subtitle = sprintf("%d Month 1 patients; %d go forward with methylation data", nrow(m1), length(final_patients)),
       x = "Patients", y = NULL) + theme_qc
save_plot(p, "04_01_cohort_flow.png", 10, 4.5)

# Fig 2: the replicate decisions - each dot is a sample, filled = kept.
rp <- rep_tab
rp$status <- ifelse(rp$kept, "Kept", "Dropped")
rp$vial_lab <- paste("Vial", rp$vial)
p <- ggplot(rp, aes(pct_missing, patient_id)) +
  geom_line(aes(group = patient_id), colour = "grey75") +
  geom_point(aes(colour = status, shape = vial_lab), size = 4, stroke = 1.1) +
  scale_colour_manual(values = c(Kept = "#185FA5", Dropped = "#EF4444"), name = NULL) +
  scale_shape_manual(values = c("Vial A" = 16, "Vial B" = 17, "Vial C" = 15), name = NULL) +
  labs(title = "Choosing one sample for each replicate patient",
       subtitle = "Rules: same vial as Month 1 > QC Pass > vial A > fewer missing values\nFull details: 04_02_replicate_resolution.csv",
       x = "Sample-specific missing values (%)", y = NULL) + theme_qc
save_plot(p, "04_02_replicate_resolution.png", 9, 5)

# ============================================================
# 7. REPRESENTATIVENESS: WITH vs WITHOUT METHYLATION
#    If patients with methylation differ from those without (older, later
#    stage, worse survival...), results from methylation may not apply to
#    the whole cohort. This is a selection-bias check.
# ============================================================
step("7. With vs without methylation")
m1$group <- factor(ifelse(m1$patient_id %in% final_patients, "With methylation", "Without methylation"),
                   levels = names(PAL_GROUP))

# PAM50 for all Month 1 patients (the table ships inside TCGAbiolinks).
m1$pam50 <- NA_character_
if (requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  st <- tryCatch(as.data.frame(suppressMessages(TCGAbiolinks::TCGAquery_subtype(tumor = "brca"))),
                 error = function(e) NULL)
  if (!is.null(st)) m1$pam50 <- st$BRCA_Subtype_PAM50[match(m1$patient_id, norm_id(st$patient))]
}
m1$pam50[m1$pam50 %in% c("NA", "")] <- NA
m1$pam50_label <- factor(ifelse(is.na(m1$pam50), "Unknown", m1$pam50), levels = PAM50_LEVELS)

num_row <- function(var, label, digits = 1) {           # numeric: median [IQR] + Wilcoxon test
  if (!var %in% names(m1)) return(NULL)
  x <- suppressWarnings(as.numeric(m1[[var]]))
  a <- x[m1$group == "With methylation"]; b <- x[m1$group == "Without methylation"]
  f <- function(v) sprintf(paste0("%.", digits, "f [%.", digits, "f-%.", digits, "f]"),
                           median(v, na.rm = TRUE), quantile(v, 0.25, na.rm = TRUE), quantile(v, 0.75, na.rm = TRUE))
  p <- tryCatch(suppressWarnings(wilcox.test(a, b)$p.value), error = function(e) NA)
  data.frame(variable = label, with_methylation = f(a), without_methylation = f(b),
             test = "Wilcoxon (median [IQR])", p_value = p)
}
pct_row <- function(var, label, positive = 1) {         # yes/no: % + chi-square test
  if (!var %in% names(m1)) return(NULL)
  x <- m1[[var]] %in% c(positive, TRUE, "TRUE", "1")
  x[is.na(m1[[var]])] <- NA
  a <- x[m1$group == "With methylation"]; b <- x[m1$group == "Without methylation"]
  p <- tryCatch(suppressWarnings(chisq.test(table(m1$group, x))$p.value), error = function(e) NA)
  data.frame(variable = label,
             with_methylation = sprintf("%.1f%% (%d/%d)", 100 * mean(a, na.rm = TRUE), sum(a, na.rm = TRUE), sum(!is.na(a))),
             without_methylation = sprintf("%.1f%% (%d/%d)", 100 * mean(b, na.rm = TRUE), sum(b, na.rm = TRUE), sum(!is.na(b))),
             test = "Chi-square (% yes)", p_value = p)
}
cat_row <- function(var, label) {                        # categories: distribution + chi-square
  if (!var %in% names(m1)) return(NULL)
  x <- as.character(m1[[var]])
  tb <- table(m1$group, x)
  p <- tryCatch(suppressWarnings(chisq.test(tb, simulate.p.value = TRUE, B = 5000)$p.value), error = function(e) NA)
  fmt <- function(g) { v <- tb[g, ]; paste(sprintf("%s %.0f%%", names(v), 100 * v / sum(v)), collapse = "; ") }
  data.frame(variable = label, with_methylation = fmt("With methylation"),
             without_methylation = fmt("Without methylation"),
             test = "Chi-square, simulated p (distribution)", p_value = p)
}
comp <- rbind(
  data.frame(variable = "Patients", with_methylation = as.character(sum(m1$group == "With methylation")),
             without_methylation = as.character(sum(m1$group == "Without methylation")), test = "", p_value = NA),
  num_row("AGE", "Age at diagnosis (years)"),
  cat_row("stage_simple", "Pathological stage"),
  pct_row("os_event", "Died during follow-up"),
  num_row("OS_MONTHS", "Follow-up time (months)"),
  num_row("tmb_nonsynonymous", "Tumour mutation burden (TMB)", digits = 2),
  pct_row("TP53_mutated", "TP53 mutated"),
  pct_row("PIK3CA_mutated", "PIK3CA mutated"),
  cat_row("pam50_label", "PAM50 subtype"))
comp$p_value    <- signif(comp$p_value, 3)
comp$difference <- ifelse(is.na(comp$p_value), "", ifelse(comp$p_value < 0.05, "DIFFERS (p<0.05)", "similar"))
cat("Patients with vs without methylation data:\n")
show(comp[, c("variable", "with_methylation", "without_methylation", "p_value", "difference")])
save_table(comp, "04_04_with_vs_without_methylation.csv")

# Fig 3: age and follow-up time, with vs without methylation.
long_num <- rbind(
  data.frame(group = m1$group, value = suppressWarnings(as.numeric(m1$AGE)),
             metric = sprintf("Age at diagnosis (p = %s)", fmt_p(comp$p_value[comp$variable == "Age at diagnosis (years)"]))),
  data.frame(group = m1$group, value = suppressWarnings(as.numeric(m1$OS_MONTHS)),
             metric = sprintf("Follow-up months (p = %s)", fmt_p(comp$p_value[comp$variable == "Follow-up time (months)"]))))
long_num <- long_num[!is.na(long_num$value), ]
p <- ggplot(long_num, aes(group, value, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.85) +
  geom_jitter(width = 0.15, height = 0, size = 0.6, alpha = 0.35, colour = "grey20") +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = PAL_GROUP, guide = "none") +
  labs(title = "Are patients with methylation data representative? (1) age and follow-up",
       subtitle = "Wilcoxon test p-values in the panel titles; p < 0.05 would suggest a difference",
       x = NULL, y = NULL) + theme_qc
save_plot(p, "04_03_age_followup_by_group.png", 10, 5.5)

# Fig 4 / 5: stage and PAM50 make-up, with vs without methylation.
comp_bar <- function(var, title, pal = NULL, pname) {
  d <- as.data.frame(table(group = m1$group, category = as.character(m1[[var]])))
  d <- d[d$Freq > 0, ]
  d$share <- ave(d$Freq, d$group, FUN = function(z) z / sum(z))
  p <- ggplot(d, aes(group, Freq, fill = category)) +
    geom_col(position = "fill", width = 0.6) +
    geom_text(aes(label = ifelse(share >= 0.06, Freq, "")), position = position_fill(vjust = 0.5),
              size = 3.3, colour = "white", fontface = "bold") +
    scale_y_continuous(labels = function(z) paste0(z * 100, "%")) +
    guides(fill = guide_legend(nrow = 2)) +
    labs(title = title, subtitle = sprintf("Numbers = patients; chi-square p = %s",
                                           fmt_p(comp$p_value[comp$variable == pname])),
         x = NULL, y = "Share of patients", fill = NULL) + theme_qc
  if (!is.null(pal)) p <- p + scale_fill_manual(values = pal)
  p
}
if ("stage_simple" %in% names(m1)) {
  p <- comp_bar("stage_simple", "Are patients with methylation data representative? (2) stage",
                pname = "Pathological stage")
  save_plot(p, "04_04_stage_by_group.png", 8, 6)
}
p <- comp_bar("pam50_label", "Are patients with methylation data representative? (3) PAM50 subtype",
              pal = PAL_PAM50, pname = "PAM50 subtype")
save_plot(p, "04_05_pam50_by_group.png", 8, 6)

# ============================================================
# 8. SAVE THE PATIENT-LEVEL OBJECT AND CHECK IT RELOADS
# ============================================================
step("8. Save patient-level object")
fp <- sum(assay(se_p, withDimnames = FALSE)[, 1], na.rm = TRUE)   # fingerprint of the first column
dims_saved <- dim(se_p)
saveRDS(se_p, OUT_RDS, compress = FALSE); record(OUT_RDS, "processed data")
cat(sprintf("Saved %s (%.1f MB)\n", OUT_RDS, file.size(OUT_RDS) / 1024^2))
rm(se_p); invisible(gc())
chk <- readRDS(OUT_RDS)
ok_reload <- identical(dim(chk), dims_saved) &&
             isTRUE(all.equal(sum(assay(chk, withDimnames = FALSE)[, 1], na.rm = TRUE), fp)) &&
             !any(duplicated(colnames(chk)))
cat("Reload check - dimensions", paste(dim(chk), collapse = " x "), "| values match and IDs unique:", ok_reload, "\n")
if (!ok_reload) stop("Saved patient-level object failed the reload check.")
rm(chk); invisible(gc())

# ============================================================
# 9. KEY INSIGHTS AND OUTPUT INDEX
# ============================================================
step("9. Key insights")
differs <- comp$variable[comp$difference == "DIFFERS (p<0.05)"]
cat("\n====================== KEY INSIGHTS ======================\n")
cat(sprintf("* ID checks: all %d passed - barcodes, patient IDs and sample IDs are consistent.\n", nrow(checks)))
cat(sprintf("* One sample per patient: %d samples -> %d patients (%d replicate samples dropped).\n",
            nrow(cd), length(final_patients), nrow(cd) - length(final_patients)))
cat(sprintf("* Replicates: %d of %d patients kept the same tumour vial as Month 1's mutation data.\n",
            sum(chosen$vial_matches_month1[chosen$patient_id %in% rep_pat]), length(rep_pat)))
cat(sprintf("* Cohort: %d of %d Month 1 patients (%.1f%%) have methylation data for integration.\n",
            length(final_patients), nrow(m1), 100 * length(final_patients) / nrow(m1)))
if (length(differs)) {
  cat("* Representativeness: patients with methylation DIFFER from those without in:",
      paste(differs, collapse = "; "), "- state this as a limitation.\n")
} else {
  cat("* Representativeness: no significant differences between patients with and without methylation.\n")
}
cat("* Next: script 05 merges this patient-level methylation data with the Month 1 clinical-genomic table.\n")
cat("==========================================================\n")

idx <- do.call(rbind, outputs_written)
idx$size_kb <- round(file.size(idx$file) / 1024, 1)
save_table(idx, "04_99_outputs_index.csv")
cat("\nFiles written:", nrow(idx), "(", sum(idx$type == "figure"), "figures,",
    sum(idx$type == "table"), "tables )\n")
cat("Total run time:", fmt_time(as.numeric(difftime(Sys.time(), run_start, units = "secs"))), "\n")
cat("Report saved to:", REPORT, "\n")
cat("\nScript 04 complete.\n")

sink()
close(report_con)
