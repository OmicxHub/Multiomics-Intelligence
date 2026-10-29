#!/usr/bin/env Rscript
# ============================================================
# 02b_assemble_beta_matrix.R
#
# Assembles the 680 downloaded TCGA-BRCA 450K beta-value files into ONE
# annotated object (probes x samples) and saves it.
#
# Script 02 downloads and md5-verifies the raw files; this script only
# READS them. It is kept separate so the download script that was run to
# fetch the data stays exactly as it was.
#
# WHY NOT TCGAbiolinks::GDCprepare()?
#   1. Memory. This machine has 7.7 GB RAM. GDCprepare() was measured at
#      725 MB peak for 10 files and 1,259 MB for 40 (~18 MB extra per
#      file), projecting to ~12.6 GB for 680 files. This script allocates
#      the final matrix once (486,427 x 680 doubles = 2.65 GB) and fills
#      it one file at a time, so peak memory stays around 3 GB.
#   2. Annotation. GDCprepare(summarizedExperiment = TRUE) annotates probes
#      via the sesame / sesameData / ExperimentHub packages (not installed;
#      it stops with "sesameData package is needed"). Those packages wrap
#      the Zhou lab InfiniumAnnotation HM450 hg38 manifests, so we read that
#      SAME source directly as plain files.
#   Files are read exactly as GDCprepare() reads them (data.table::fread,
#   no header, character probe ID + numeric beta). Pass --gdcprepare to use
#   GDCprepare() for the matrix on a machine with >= 16 GB RAM.
#
# OUTPUT: data/processed/methylation_beta_matrix_raw.rds
#   A SummarizedExperiment:
#     assay "beta"  probes x samples (486,427 x 680), values 0-1, NA = masked
#     rowData       probe annotation: type (cg/ch/rs/ct), Infinium design
#                   (I/II), hg38 chr + position, gene, distance to TSS,
#                   CpG-island context, genomic region, quality mask
#     colData       sample annotation: barcode fields (tissue source site,
#                   vial, plate = processing batch), replicate flags,
#                   PAM50 subtype, pathology, TCGA methylation cluster
#
# Run from the Ubuntu terminal with:
#   Rscript scripts/02b_assemble_beta_matrix.R
#
# Useful flags:
#   --pilot=40      first 40 files + all replicate-patient files, saved to a
#                   separate *_pilot40.rds (never overwrites the real output)
#   --verify-md5    re-hash every raw file before assembling (slow over USB;
#                   script 02 has already verified them)
#   --gdcprepare    build the matrix with GDCprepare() (needs >= 16 GB RAM)
#   --purge-raw     delete data/raw/ after the .rds is saved and reloads
#                   correctly (OFF by default)
#
# Prerequisite: 02_download_beta_matrix.R has downloaded all 680 files.
# ============================================================

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(data.table)
  library(matrixStats)
})

# ---- 0. Working directory, flags, logging ----
setwd("/mnt/e/R_Studio/Projects/MOI/Month_2")

args    <- commandArgs(trailingOnly = TRUE)
getflag <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}
PILOT          <- suppressWarnings(as.integer(getflag("pilot", NA)))
VERIFY_MD5     <- "--verify-md5" %in% args
USE_GDCPREPARE <- "--gdcprepare" %in% args
PURGE_RAW      <- "--purge-raw" %in% args

RAW_DIR  <- "data/raw"
MANIFEST <- "data/processed/methylation_matched_manifest.csv"
OUT_RDS  <- "data/processed/methylation_beta_matrix_raw.rds"
REF_DIR  <- "data/reference"
LOGFILE  <- "outputs/reports/02b_assembly_log.txt"
if (!is.na(PILOT)) {
  OUT_RDS <- sub("\\.rds$", paste0("_pilot", PILOT, ".rds"), OUT_RDS)
  LOGFILE <- sub("\\.txt$", paste0("_pilot", PILOT, ".txt"), LOGFILE)
}

# Everything printed is also written to the log file.
dir.create(dirname(LOGFILE), recursive = TRUE, showWarnings = FALSE)
log_con <- file(LOGFILE, open = "wt")
sink(log_con, split = TRUE)

run_start <- Sys.time()
fmt_time <- function(secs) {
  if (!is.finite(secs)) return("--")
  if (secs < 60) return(sprintf("%.0f s", secs))
  if (secs < 3600) return(sprintf("%.1f min", secs / 60))
  sprintf("%.1f h", secs / 3600)
}
step_start <- Sys.time()
step <- function(title) {
  now <- Sys.time()
  cat(sprintf("\n[%s | total %s] ==== %s ====\n", format(now, "%H:%M:%S"),
              fmt_time(as.numeric(difftime(now, run_start, units = "secs"))), title))
  step_start <<- now
}

cat("Working directory :", getwd(), "\n")
cat("Started           :", format(run_start), "\n")
cat("Settings          : pilot =", ifelse(is.na(PILOT), "no (all files)", PILOT),
    "| verify md5 =", VERIFY_MD5, "| assembler =",
    ifelse(USE_GDCPREPARE, "GDCprepare()", "lean"), "| purge raw =", PURGE_RAW, "\n")
cat("Output            :", OUT_RDS, "\n")
cat("Log               :", LOGFILE, "\n")

# ---- 1. Manifest ----
step("1. Load manifest")
if (!file.exists(MANIFEST)) stop("Manifest not found: ", MANIFEST, " - run script 01 first.")
manifest <- read.csv(MANIFEST, stringsAsFactors = FALSE)
need <- c("id", "file_id", "file_name", "md5sum", "file_size", "cases", "project",
          "data_category", "data_type", "patient_id", "sample.submitter_id", "sample_type")
if (length(setdiff(need, colnames(manifest))))
  stop("Manifest missing column(s): ", paste(setdiff(need, colnames(manifest)), collapse = ", "))

rep_patients <- names(which(table(manifest$patient_id) > 1))
cat("Files in manifest   :", nrow(manifest), "\n")
cat("Unique patients     :", length(unique(manifest$patient_id)), "\n")
cat("Replicate patients  :", length(rep_patients), "(",
    sum(manifest$patient_id %in% rep_patients), "files ) - kept here; collapsed to one",
    "sample per patient in a later script\n")

if (!is.na(PILOT)) {
  rows <- sort(unique(c(seq_len(min(PILOT, nrow(manifest))),
                        which(manifest$patient_id %in% rep_patients))))
  manifest <- manifest[rows, , drop = FALSE]
  cat(">>> PILOT MODE: first", PILOT, "files + replicate patients =", nrow(manifest), "files\n")
}
n_s <- nrow(manifest)

# ---- 2. Check the raw files ----
step("2. Check raw files")
dest <- file.path(RAW_DIR, manifest$project,
                  gsub(" ", "_", manifest$data_category),
                  gsub(" ", "_", manifest$data_type),
                  gsub(" ", "_", manifest$file_id),
                  gsub(" ", "_", manifest$file_name))
present <- file.exists(dest)
size_ok <- present & file.size(dest) == as.numeric(manifest$file_size)
cat("Files present       :", sum(present), "of", n_s, "\n")
cat("Exact size match    :", sum(size_ok), "of", n_s, "\n")
if (!all(size_ok)) {
  print(utils::head(manifest[!size_ok, c("file_id", "cases")], 10))
  stop(sum(!size_ok), " raw file(s) missing or wrong size. ",
       "Run: Rscript scripts/02_download_beta_matrix.R --download-only")
}
if (VERIFY_MD5) {
  cat("Re-hashing all files (--verify-md5)...\n")
  md5_ok <- logical(n_s)
  for (i in seq_len(n_s)) {
    md5_ok[i] <- unname(tools::md5sum(dest[i])) == manifest$md5sum[i]
    if (i %% 50 == 0 || i == n_s) cat("   hashed", i, "/", n_s, "| passing:", sum(md5_ok[1:i]), "\n")
  }
  if (!all(md5_ok)) stop(sum(!md5_ok), " file(s) failed md5 - re-run script 02 --download-only.")
} else {
  cat("md5 not re-checked (script 02 verified all 680; use --verify-md5 to re-hash).\n")
}

# ---- 3. Probe annotation reference files ----
step("3. Annotation reference files")
REF_BASE <- "https://raw.githubusercontent.com/zhou-lab/InfiniumAnnotationV1/main/Anno/HM450/"
REF <- c(gencode = "HM450.hg38.manifest.gencode.v36.tsv.gz",
         design  = "HM450.hg38.manifest.tsv.gz",
         mask    = "HM450.hg38.mask.tsv.gz")
dir.create(REF_DIR, recursive = TRUE, showWarnings = FALSE)
for (f in REF) {
  p <- file.path(REF_DIR, f)
  if (!file.exists(p) || file.size(p) == 0) {
    cat("Downloading", f, "...\n")
    st <- system2("curl", c("-sS", "-L", "--fail", "--retry", "5", "--max-time", "900",
                            "-o", shQuote(p), shQuote(paste0(REF_BASE, f))))
    if (st != 0 || !file.exists(p)) stop("Could not download annotation file ", f)
  }
  cat(sprintf("  %-42s %6.1f MB\n", f, file.size(p) / 1024^2))
}

# ---- 4. Probe list + probe annotation (rowData) ----
# Built BEFORE the big matrix is allocated: the annotation step briefly
# needs ~1.5 GB of working memory, and doing it first keeps the peak at
# ~2.9 GB instead of ~4.2 GB (matrix + annotation at the same time).
step("4. Probe annotation")
read_betas <- function(f) fread(f, header = FALSE, sep = "\t",
                                colClasses = c("character", "numeric"),
                                showProgress = FALSE)
ref_probes <- read_betas(dest[1])$V1
cat("Probes per file     :", length(ref_probes), "(taken from the first file)\n")

gen <- fread(file.path(REF_DIR, REF[["gencode"]]),
             select = c("probeID", "CpG_chrm", "CpG_beg", "genesUniq", "distToTSS", "CGIposition"))
des <- fread(file.path(REF_DIR, REF[["design"]]),
             select = c("Probe_ID", "CpG_chrm", "CpG_beg", "type"))
msk <- fread(file.path(REF_DIR, REF[["mask"]]),
             select = c("Probe_ID", "maskUniq", "M_general"))

pid <- ref_probes
ig  <- match(pid, gen$probeID)
id  <- match(pid, des$Probe_ID)
im  <- match(pid, msk$Probe_ID)

# distToTSS lists one distance per overlapping transcript ("-1390;10744");
# keep the distance nearest the TSS.
nearest_tss <- function(x) {
  out <- rep(NA_integer_, length(x))
  has <- !is.na(x) & nzchar(x)
  out[has] <- vapply(strsplit(x[has], ";", fixed = TRUE), function(v) {
    v <- suppressWarnings(as.integer(v)); v <- v[!is.na(v)]
    if (length(v)) v[which.min(abs(v))] else NA_integer_
  }, integer(1))
  out
}

chr     <- ifelse(is.na(ig), des$CpG_chrm[id], gen$CpG_chrm[ig])
beg0    <- ifelse(is.na(ig), des$CpG_beg[id],  gen$CpG_beg[ig])
cgi_raw <- gen$CGIposition[ig]
gene    <- gen$genesUniq[ig]
dtss    <- nearest_tss(gen$distToTSS[ig])

rd <- data.frame(
  probe_id     = pid,
  probe_type   = substr(pid, 1, 2),                    # cg / ch / rs / ct
  design_type  = des$type[id],                         # Infinium I / II
  chr          = chr,
  pos_hg38     = ifelse(is.na(beg0), NA_integer_, as.integer(beg0) + 1L),
  gene         = gene,
  dist_to_tss  = dtss,
  cgi_position = cgi_raw,                              # Island / N_Shore / S_Shelf ...
  cgi_context  = ifelse(is.na(chr), NA_character_,
                 ifelse(is.na(cgi_raw), "OpenSea",
                 ifelse(cgi_raw == "Island", "Island",
                 ifelse(grepl("Shore", cgi_raw), "Shore", "Shelf")))),
  region       = ifelse(is.na(chr), NA_character_,
                 ifelse(!is.na(dtss) & abs(dtss) <= 1500, "Promoter (TSS +/-1500)",
                 ifelse(!is.na(gene) & nzchar(gene), "Gene body", "Intergenic"))),
  mask_general = ifelse(is.na(im), FALSE, msk$M_general[im]),
  mask_reason  = msk$maskUniq[im],
  stringsAsFactors = FALSE,
  row.names = pid
)
rm(gen, des, msk, pid, ig, id, im, chr, beg0, cgi_raw, gene, dtss); invisible(gc())

cat("Probes with hg38 coordinates  :", sum(!is.na(rd$chr)), "of", nrow(rd), "\n")
cat("Flagged by quality mask       :", sum(rd$mask_general), "\n")
cat("\nProbe type x design:\n"); print(table(type = rd$probe_type, design = rd$design_type, useNA = "ifany"))
cat("\nCpG-island context:\n");  print(table(rd$cgi_context, useNA = "ifany"))
cat("\nGenomic region:\n");      print(table(rd$region, useNA = "ifany"))

# ---- 5. Sample annotation (colData) ----
step("5. Sample annotation")
bc    <- manifest$cases
parts <- do.call(rbind, strsplit(bc, "-", fixed = TRUE))
n_per_patient <- as.integer(table(manifest$patient_id)[manifest$patient_id])
cd <- data.frame(
  barcode       = bc,
  patient_id    = manifest$patient_id,
  sample_id     = manifest$sample.submitter_id,
  file_id       = manifest$file_id,
  sample_type   = manifest$sample_type,
  tss           = parts[, 2],                  # tissue source site
  vial          = substr(parts[, 4], 3, 3),    # A / B ...
  portion       = substr(parts[, 5], 1, 2),
  analyte       = substr(parts[, 5], 3, 3),
  plate         = parts[, 6],                  # processing batch
  center        = parts[, 7],
  n_samples_for_patient = n_per_patient,
  is_replicate_patient  = n_per_patient > 1,
  stringsAsFactors = FALSE,
  row.names = bc
)
# Published TCGA-BRCA subtype calls ship inside TCGAbiolinks (no network needed).
st <- tryCatch(as.data.frame(suppressMessages(TCGAquery_subtype(tumor = "brca"))),
               error = function(e) NULL)
if (!is.null(st)) {
  j <- match(cd$patient_id, st$patient)
  cd$pam50             <- st$BRCA_Subtype_PAM50[j]
  cd$pathology         <- st$BRCA_Pathology[j]
  cd$meth_cluster_tcga <- as.character(st[["DNA.Methylation Clusters"]][j])
  cat("PAM50 available for", sum(!is.na(cd$pam50)), "of", nrow(cd), "samples\n")
  print(table(PAM50 = cd$pam50, useNA = "ifany"))
} else {
  cat("TCGAquery_subtype() unavailable - subtype columns not added.\n")
}
rm(st, parts); invisible(gc())
cat("\nProcessing plates (batches):", length(unique(cd$plate)),
    "| tissue source sites:", length(unique(cd$tss)), "\n")

# ---- 6. Assemble the beta matrix ----
step("6. Assemble beta matrix")
cat("Close other programs (especially web browsers): this step needs ~3 GB free RAM.\n")
if (USE_GDCPREPARE) {
  QUERY_CACHE <- "data/processed/methylation_query_object.rds"
  if (!file.exists(QUERY_CACHE)) stop("GDCquery cache not found: ", QUERY_CACHE)
  query <- readRDS(QUERY_CACHE)
  res   <- query$results[[1]]
  query$results[[1]] <- res[match(manifest$id, res$file_id), , drop = FALSE]
  stopifnot(identical(query$results[[1]]$file_id, manifest$id))
  beta <- GDCprepare(query, directory = RAW_DIR, summarizedExperiment = FALSE)
  if (!identical(colnames(beta), manifest$cases)) beta <- beta[, manifest$cases, drop = FALSE]
  if (!identical(rownames(beta), rownames(rd)))   beta <- beta[rownames(rd), , drop = FALSE]
  ASSEMBLY <- "TCGAbiolinks::GDCprepare(summarizedExperiment = FALSE)"
} else {
  cat(sprintf("Allocating matrix   : %d x %d = %.2f GB\n",
              length(ref_probes), n_s, length(ref_probes) * n_s * 8 / 1024^3))
  beta <- matrix(NA_real_, nrow = length(ref_probes), ncol = n_s,
                 dimnames = list(ref_probes, manifest$cases))

  n_realigned <- 0
  t_asm <- Sys.time()
  for (i in seq_len(n_s)) {
    d <- read_betas(dest[i])
    if (identical(d$V1, ref_probes)) {
      beta[, i] <- d$V2
    } else {
      # Different probe order/content: align by probe ID, not by position.
      beta[, i] <- d$V2[match(ref_probes, d$V1)]
      n_realigned <- n_realigned + 1
    }
    if (i %% 10 == 0 || i == n_s) {
      el   <- as.numeric(difftime(Sys.time(), t_asm, units = "secs"))
      frac <- i / n_s
      fill <- round(frac * 30)
      cat(sprintf("\r  [%s%s] %5.1f%% | %d/%d files | elapsed %-8s | ETA %-8s",
                  strrep("=", fill), strrep(".", 30 - fill), 100 * frac, i, n_s,
                  fmt_time(el), fmt_time(el / frac - el)))
      flush(stdout())
    }
  }
  cat("\n")
  rm(d); invisible(gc())
  cat("Files needing probe re-alignment:", n_realigned, "\n")
  ASSEMBLY <- "lean assembler (pre-allocated matrix, filled file by file)"
}

# ---- 7. Verify the matrix against the source files ----
step("7. Verify matrix")
cat("Dimensions (probes x samples) :", paste(dim(beta), collapse = " x "), "\n")
cat("Columns = manifest barcodes   :", identical(colnames(beta), manifest$cases), "\n")
rng <- range(beta, na.rm = TRUE)
cat(sprintf("Beta value range              : %.4f to %.4f\n", rng[1], rng[2]))
if (rng[1] < 0 || rng[2] > 1) stop("Beta values outside [0,1] - assembly is wrong.")
na_pct <- colCounts(beta, value = NA) / nrow(beta) * 100
cat(sprintf("Missing values per sample     : median %.1f%%, range %.1f%% - %.1f%%\n",
            median(na_pct), min(na_pct), max(na_pct)))
cat(sprintf("Missing values overall        : %.2f%%\n", mean(na_pct)))

set.seed(1)
spot <- sort(sample(seq_len(n_s), min(5, n_s)))
spot_ok <- vapply(spot, function(i) {
  d <- fread(dest[i], header = FALSE, sep = "\t",
             colClasses = c("character", "numeric"), showProgress = FALSE)
  isTRUE(all.equal(unname(beta[, i]), d$V2[match(rownames(beta), d$V1)]))
}, logical(1))
cat("Spot-check (re-read files)    :", sum(spot_ok), "of", length(spot), "exact matches\n")
if (!all(spot_ok)) stop("Spot-check failed: matrix columns do not match their source files.")

# ---- 8. Build, save and reload-check the object ----
step("8. Save and verify")
se <- SummarizedExperiment(assays = list(beta = beta), rowData = rd, colData = cd)
S4Vectors::metadata(se) <- list(
  description = "TCGA-BRCA Illumina HM450 beta values (GDC SeSAMe level-3), Primary Tumor, matched to Month 1 cohort",
  assembly    = ASSEMBLY,
  annotation  = "Zhou lab InfiniumAnnotation HM450 hg38 (GENCODE v36) + general quality mask",
  manifest    = MANIFEST,
  n_files     = n_s,
  created     = format(Sys.time()),
  R_version   = R.version.string
)
print(se)

fingerprint <- colSums(beta[, spot, drop = FALSE], na.rm = TRUE)
dir.create(dirname(OUT_RDS), recursive = TRUE, showWarnings = FALSE)
cat("\nSaving", OUT_RDS, "(uncompressed: beta values barely compress, and\n",
    "uncompressed saves and loads faster with less memory)...\n")
t_save <- Sys.time()
saveRDS(se, OUT_RDS, compress = FALSE)
cat(sprintf("Saved: %.1f MB in %s\n", file.size(OUT_RDS) / 1024^2,
            fmt_time(as.numeric(difftime(Sys.time(), t_save, units = "secs")))))

cat("Reloading to confirm the file is complete and readable...\n")
rm(se, beta); invisible(gc())
t_load <- Sys.time()
check  <- readRDS(OUT_RDS)
reload_dims <- dim(check)
same <- isTRUE(all.equal(colSums(assay(check)[, spot, drop = FALSE], na.rm = TRUE), fingerprint))
cat(sprintf("Reloaded in %s | dimensions %s | values match pre-save fingerprint: %s\n",
            fmt_time(as.numeric(difftime(Sys.time(), t_load, units = "secs"))),
            paste(reload_dims, collapse = " x "), same))
if (!same || reload_dims[2] != n_s) stop("Reloaded object does not match what was saved.")
rm(check); invisible(gc())

# ---- 9. Optional: remove the raw download cache ----
if (PURGE_RAW) {
  if (is.na(PILOT)) {
    cat("\n--purge-raw: .rds verified, deleting", RAW_DIR, "...\n")
    unlink(RAW_DIR, recursive = TRUE, force = TRUE)
    dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
    cat("Raw files removed; only the assembled .rds is kept.\n")
  } else {
    cat("\n--purge-raw ignored in pilot mode.\n")
  }
} else {
  cat("\ndata/raw/ kept. Delete it later with --purge-raw once script 03 has run successfully.\n")
}

# ---- 10. Summary ----
step("Done")
mem <- gc()
cat("\n===== SCRIPT 02b SUMMARY =====\n")
cat("Samples assembled  :", n_s, "| patients:", length(unique(manifest$patient_id)), "\n")
cat("Probes             :", reload_dims[1], "\n")
cat("Assembly method    :", ASSEMBLY, "\n")
cat("Peak R memory      :", round(sum(mem[, ncol(mem)]) / 1024, 2), "GB\n")
cat("Total run time     :", fmt_time(as.numeric(difftime(Sys.time(), run_start, units = "secs"))), "\n")
cat("Output             :", OUT_RDS, "\n")
cat("Log                :", LOGFILE, "\n")
cat("==============================\n")
cat("\nScript 02b complete.\n")

sink()
close(log_con)
