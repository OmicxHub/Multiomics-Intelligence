#!/usr/bin/env Rscript
# ============================================================
# 02_download_beta_matrix.R
#
# Downloads the 680 matched TCGA-BRCA 450K methylation beta-value
# files from the GDC and assembles them into a single object.
#
# WHY THIS DOES NOT USE GDCdownload():
#   The default TCGAbiolinks method POSTs a batch of UUIDs to the GDC
#   /data endpoint and receives ONE continuous streamed .tar.gz per
#   chunk. That stream has no resume, no per-file checksum, and no
#   integrity check until the very end - so a single garbled segment
#   anywhere kills the whole chunk ("truncated gzip input" /
#   "tar: Unexpected EOF in archive"), and retrying restarts from
#   byte zero. method="client" fails separately because no
#   gdc-client binary is on PATH ("argument is of length zero").
#
#   Instead we fetch each file individually from
#       https://api.gdc.cancer.gov/data/{file_id}
#   which is verified to support HTTP range requests (resume) and to
#   return plain UNCOMPRESSED text - no tar, no gzip, nothing that
#   can be "truncated". Every file is checked against its md5sum from
#   the manifest, and only files that fail are retried.
#
# Run from the Ubuntu terminal with:
#   Rscript scripts/02_download_beta_matrix.R
#
# Useful flags:
#   --pilot=5        only do the first 5 files (end-to-end smoke test)
#   --jobs=4         parallel curl transfers (default 4)
#   --batch=20       files per reporting batch (default 20)
#   --passes=4       max retry passes over failed files (default 4)
#   --download-only  fetch and verify, but stop before GDCprepare()
#   --purge-raw      delete data/raw/ after the .rds is verified
#                    (OFF by default - see section 10)
#
# EXPECT THIS TO TAKE A LONG TIME. Measured throughput from this machine
# is ~30-70 KB/s to the GDC (and to AWS us-east-1 generally), versus
# ~430 KB/s to a CDN with a local edge node - so the limit is the network
# path to US-hosted research data, not this script. More parallelism does
# not help: 16 connections measured SLOWER than 4. At that rate the full
# 8.27 GiB is on the order of 1.5-3 days of transfer.
#
# The script is built for exactly that: it is fully resumable. Stop it
# whenever you like (Ctrl-C, reboot, unplugging the drive) and re-run the
# same command - verified files are kept, partial files resume mid-byte,
# and only what is missing is fetched.
#
# Prerequisite: 01_prepare_methylation_data.R must have completed.
# ============================================================

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(tools)
})

# ---- 0. Working directory + argument parsing ----
setwd("/mnt/e/R_Studio/Projects/MOI/Month_2")
cat("Working directory set to:", getwd(), "\n")

args    <- commandArgs(trailingOnly = TRUE)
getflag <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}
PILOT     <- suppressWarnings(as.integer(getflag("pilot", NA)))
JOBS      <- suppressWarnings(as.integer(getflag("jobs", "4")))
BATCH     <- suppressWarnings(as.integer(getflag("batch", "20")))
MAX_PASS  <- suppressWarnings(as.integer(getflag("passes", "4")))
PURGE_RAW <- "--purge-raw" %in% args
DL_ONLY   <- "--download-only" %in% args

RAW_DIR  <- "data/raw"
OUT_RDS  <- "data/processed/methylation_beta_matrix_raw.rds"
MANIFEST <- "data/processed/methylation_matched_manifest.csv"
API      <- "https://api.gdc.cancer.gov/data/"

cat("Settings: jobs =", JOBS, "| batch =", BATCH, "| max passes =", MAX_PASS,
    "| pilot =", ifelse(is.na(PILOT), "no (all files)", PILOT),
    "| purge raw =", PURGE_RAW, "| download only =", DL_ONLY, "\n\n")

# ---- 1. Load and validate the matched manifest from Script 01 ----
if (!file.exists(MANIFEST)) {
  stop("Matched manifest not found at ", MANIFEST,
       "\nRun 01_prepare_methylation_data.R first.")
}
manifest <- read.csv(MANIFEST, stringsAsFactors = FALSE)

need_cols <- c("id", "file_id", "file_name", "md5sum", "file_size",
               "cases", "project", "data_category", "data_type",
               "platform", "patient_id")
missing_cols <- setdiff(need_cols, colnames(manifest))
if (length(missing_cols)) {
  stop("Manifest is missing required column(s): ",
       paste(missing_cols, collapse = ", "))
}

cat("===== MANIFEST LOADED =====\n")
cat("Files to download   :", nrow(manifest), "\n")
cat("Unique file UUIDs   :", length(unique(manifest$id)), "\n")
cat("Unique barcodes     :", length(unique(manifest$cases)), "\n")
cat("Unique PATIENTS     :", length(unique(manifest$patient_id)), "\n")
cat("Total download size :",
    round(sum(as.numeric(manifest$file_size)) / 1024^3, 2), "GiB\n")

# NOTE: 680 files but only 672 patients - 6 patients contributed more
# than one Primary Tumour aliquot/vial. All 680 barcodes are distinct,
# so the GDCprepare() duplicate guard (which tests duplicated(cases))
# will NOT trip. Collapsing to one file per patient is deliberately
# left to the cleaning step so this raw object stays faithful to GDC.
dupe_patients <- sum(table(manifest$patient_id) > 1)
if (dupe_patients > 0) {
  cat("NOTE:", dupe_patients, "patients have >1 file (replicate vials/aliquots).\n")
  cat("      Raw object will have", nrow(manifest), "columns for",
      length(unique(manifest$patient_id)), "patients.\n")
  cat("      De-duplicate to one sample per patient in script 03 BEFORE\n")
  cat("      merging with the Month 1 cohort, or those patients will\n")
  cat("      fan out into duplicate rows.\n")
}
cat("===========================\n\n")

if (!is.na(PILOT)) {
  manifest <- manifest[seq_len(min(PILOT, nrow(manifest))), , drop = FALSE]
  cat(">>> PILOT MODE: restricted to", nrow(manifest), "files.\n\n")
}

# ---- 2. Build destination paths in the exact layout GDCprepare() expects ----
# GDCprepare() reconstructs its file paths as:
#   {directory}/{project}/{data_category}/{data_type}/{file_id}/{file_name}
# with spaces replaced by underscores in every component except project.
# Deriving the paths from the same manifest columns guarantees a match.
dest_dir <- file.path(RAW_DIR,
                      manifest$project,
                      gsub(" ", "_", manifest$data_category),
                      gsub(" ", "_", manifest$data_type),
                      gsub(" ", "_", manifest$file_id))
dest <- file.path(dest_dir, gsub(" ", "_", manifest$file_name))

cat("Files will be written under:\n  ",
    dirname(dest_dir[1]), "/{file_id}/\n\n", sep = "")
invisible(lapply(unique(dest_dir), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

# ---- 2b. Single-instance lock ----
# Two copies of this script running at once will have their curls write
# the SAME files concurrently, interleaving bytes and producing corrupt
# (even oversized) output. Accidentally launching it twice is easy, so
# refuse to start if another copy is already running.
LOCK <- file.path(RAW_DIR, ".download.lock")
if (file.exists(LOCK)) {
  old_pid <- suppressWarnings(readLines(LOCK, warn = FALSE)[1])
  alive   <- !is.na(old_pid) && nzchar(old_pid) &&
             dir.exists(file.path("/proc", old_pid))
  if (alive) {
    stop("ANOTHER DOWNLOAD IS ALREADY RUNNING (process ", old_pid, ").\n\n",
         "Do not start a second copy - both would write the same files and\n",
         "corrupt them. Watch the running one instead:\n",
         "  tail -f outputs/reports/02_overnight.log\n\n",
         "If you are certain it is dead, delete ", LOCK, " and re-run.")
  }
  cat("Found a stale lock from process", old_pid,
      "(no longer running) - clearing it.\n")
  unlink(LOCK)
}
writeLines(as.character(Sys.getpid()), LOCK)
on.exit(unlink(LOCK), add = TRUE)
cat("Lock acquired (process", Sys.getpid(), ").\n\n")

# ---- 3. Helpers: drive check + verification ----
drive_alive <- function() {
  # E: is an external USB drive - fail loudly rather than thrash if it drops.
  ok <- dir.exists(RAW_DIR) && file.access(RAW_DIR, mode = 2) == 0
  if (!ok) {
    stop("Cannot write to ", normalizePath(RAW_DIR, mustWork = FALSE),
         " - has the external E: drive disconnected?\n",
         "Reconnect it and re-run this script; verified files are kept ",
         "and only the missing ones will be fetched.")
  }
  invisible(TRUE)
}

verify_one <- function(i) {
  # TRUE only if the file exists, is exactly the expected size, and its
  # md5 matches the manifest.
  if (!file.exists(dest[i])) return(FALSE)
  if (file.size(dest[i]) != as.numeric(manifest$file_size[i])) return(FALSE)
  unname(tools::md5sum(dest[i])) == manifest$md5sum[i]
}

verify_all <- function(idx, label) {
  cat(label, "- verifying", length(idx), "file(s) against manifest md5sums...\n")
  ok <- logical(length(idx))
  for (k in seq_along(idx)) {
    ok[k] <- verify_one(idx[k])
    if (k %% 50 == 0 || k == length(idx)) {
      cat("   verified", k, "/", length(idx),
          "|  passing:", sum(ok[seq_len(k)]), "\n")
    }
  }
  ok
}

# ---- 4. Work out what still needs downloading (idempotent / resumable) ----
drive_alive()
cat("Checking for files already present from a previous run...\n")
already_ok <- verify_all(seq_len(nrow(manifest)), "Pre-flight")
pending    <- which(!already_ok)

cat("\nAlready complete & verified:", sum(already_ok), "\n")
cat("Still to download          :", length(pending), "\n\n")

# ---- 5. Download loop: batched passes, retrying ONLY the failures ----
# Measured throughput to the GDC from this machine is roughly 30-70 KB/s
# and does NOT improve with more parallel connections (16 connections
# measured slower than 4), so this is a long transfer. Files are fetched
# in small batches rather than one huge curl invocation so that progress,
# cumulative totals and an ETA are printed regularly, and so that every
# batch is checksum-verified as it lands. Stopping the script at any
# point is safe - re-running picks up exactly where it left off.
total_bytes <- sum(as.numeric(manifest$file_size))
run_start   <- Sys.time()
fail_count  <- integer(nrow(manifest))   # consecutive md5 failures per file

fmt_eta <- function(secs) {
  if (!is.finite(secs)) return("--")
  if (secs < 3600) return(sprintf("%.0f min", secs / 60))
  if (secs < 86400) return(sprintf("%.1f h", secs / 3600))
  sprintf("%.1f days", secs / 86400)
}

# Live single-line progress bar, redrawn in place with a carriage return.
draw_bar <- function(done_bytes, done_n, rate_bps) {
  frac <- if (total_bytes > 0) done_bytes / total_bytes else 0
  w    <- 28
  fill <- max(0, min(w, round(frac * w)))
  bar  <- paste0("[", strrep("=", fill), strrep(".", w - fill), "]")
  eta  <- if (rate_bps > 0) (total_bytes - done_bytes) / rate_bps else Inf
  cat(sprintf(
    "\r  %s %5.1f%% | %d/%d files | %.2f/%.2f GiB | %6.1f KB/s | ETA %-9s",
    bar, 100 * frac, done_n, nrow(manifest),
    done_bytes / 1024^3, total_bytes / 1024^3,
    rate_bps / 1024, fmt_eta(eta)))
  flush(stdout())
}

# Appends a durable one-line record per batch, so a multi-day run leaves
# a readable log even though the live bar itself is not log-friendly.
LOGFILE <- "outputs/reports/02_download_log.txt"
dir.create(dirname(LOGFILE), recursive = TRUE, showWarnings = FALSE)
logmsg <- function(...) {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
  cat(line, "\n", file = LOGFILE, sep = "", append = TRUE)
}

# Runs one batch of transfers. curl is launched in the BACKGROUND and R
# polls the growing files on disk so it can redraw a live progress bar;
# a sentinel file tells us when curl has exited.
run_batch <- function(idx, base_bytes, base_n) {
  # A file already at/over its expected size but failing md5 is
  # unrecoverable, so remove it and start clean. Short files are left in
  # place so that curl --continue-at - can resume them.
  for (i in idx) {
    if (file.exists(dest[i]) &&
        file.size(dest[i]) >= as.numeric(manifest$file_size[i])) {
      unlink(dest[i])
    }
  }

  cfg <- tempfile(fileext = ".curlcfg")
  writeLines(
    as.vector(rbind(
      paste0("url = \"",    API, manifest$id[idx], "\""),
      paste0("output = \"", normalizePath(dest[idx], mustWork = FALSE), "\"")
    )),
    cfg
  )

  curl_opts <- paste(
    "--parallel", "--parallel-max", JOBS,
    "--continue-at", "-",          # resume partial transfers
    "--location",                  # follow redirects
    "--fail",                      # non-zero exit on HTTP errors
    "--silent", "--show-error",    # quiet: R draws the progress bar
    "--retry", "5",
    "--retry-delay", "5",
    "--retry-connrefused",
    "--connect-timeout", "30",
    "--speed-limit", "1000",       # abort a transfer stalled under
    "--speed-time", "120",         #   1 KB/s for 2 min, so it gets retried
    "--max-time", "3600",
    "--config", shQuote(cfg))

  done_flag <- tempfile(fileext = ".done")
  sh        <- tempfile(fileext = ".sh")
  writeLines(c("#!/usr/bin/env bash",
               paste("curl", curl_opts),
               paste("echo $? >", shQuote(done_flag))), sh)

  system2("bash", shQuote(sh), wait = FALSE)

  # Poll until curl signals completion, redrawing the bar as bytes land.
  repeat {
    cur_bytes  <- sum(file.size(dest[idx]), na.rm = TRUE)
    done_bytes <- base_bytes + cur_bytes
    elapsed    <- as.numeric(difftime(Sys.time(), run_start, units = "secs"))
    rate       <- if (elapsed > 0) (done_bytes - start_bytes) / elapsed else 0
    draw_bar(done_bytes, base_n, rate)
    if (file.exists(done_flag)) break
    Sys.sleep(3)
  }
  status <- suppressWarnings(as.integer(readLines(done_flag, warn = FALSE)[1]))
  unlink(c(cfg, sh, done_flag))
  cat("\n")
  status
}

pass <- 1
while (length(pending) > 0 && pass <= MAX_PASS) {
  drive_alive()
  cat("========== DOWNLOAD PASS", pass, "of", MAX_PASS, "==========\n")
  cat("Outstanding:", length(pending), "file(s),",
      round(sum(as.numeric(manifest$file_size[pending])) / 1024^3, 2), "GiB\n")
  cat("Fetching in batches of", BATCH, "with", JOBS, "parallel transfers.\n\n")

  batches <- split(pending, ceiling(seq_along(pending) / BATCH))
  still   <- integer(0)

  # Bytes/files already banked before this pass, so the bar shows overall
  # progress across the whole 680, not just the current batch.
  verified_now <- vapply(seq_len(nrow(manifest)), verify_one, logical(1))
  base_bytes   <- sum(as.numeric(manifest$file_size[verified_now]))
  base_n       <- sum(verified_now)
  if (!exists("start_bytes")) start_bytes <- base_bytes

  for (b in seq_along(batches)) {
    drive_alive()
    idx <- batches[[b]]
    t0  <- Sys.time()

    cat(sprintf("  batch %d/%d (%d files)\n", b, length(batches), length(idx)))
    run_batch(idx, base_bytes, base_n)
    mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

    ok_now <- vapply(idx, verify_one, logical(1))
    still  <- c(still, idx[!ok_now])

    # A file that keeps failing may have a corrupt PREFIX, in which case
    # resuming it with --continue-at appends good bytes onto bad ones and
    # it can never pass md5. After a second failure, delete it so the next
    # pass starts that file from byte zero.
    for (i in idx[!ok_now]) {
      fail_count[i] <- fail_count[i] + 1
      if (fail_count[i] >= 2 && file.exists(dest[i])) {
        unlink(dest[i])
        cat("  (discarding repeatedly-failing partial:",
            manifest$file_id[i], "- will restart it cleanly)\n")
        fail_count[i] <- 0
      }
    }

    base_n     <- base_n + sum(ok_now)
    base_bytes <- base_bytes + sum(as.numeric(manifest$file_size[idx[ok_now]]))

    msg <- sprintf(
      "batch %d/%d | %d/%d verified | %.1f min | TOTAL %d/%d files (%.1f%%) | %.2f/%.2f GiB",
      b, length(batches), sum(ok_now), length(idx), mins,
      base_n, nrow(manifest), 100 * base_n / nrow(manifest),
      base_bytes / 1024^3, total_bytes / 1024^3)
    cat("  ", msg, "\n\n", sep = "")
    logmsg(msg)
  }

  pending <- still
  cat("\nPass", pass, "finished. Still outstanding:", length(pending), "\n\n")

  if (length(pending) > 0 && pass < MAX_PASS) {
    cat("Waiting 15s before retrying the failures...\n\n")
    Sys.sleep(15)
  }
  pass <- pass + 1
}

# ---- 6. Hard gate: refuse to assemble from an incomplete download ----
final_ok <- verify_all(seq_len(nrow(manifest)), "Final")
cat("\n===== DOWNLOAD SUMMARY =====\n")
cat("Expected files :", nrow(manifest), "\n")
cat("Verified files :", sum(final_ok), "\n")
cat("Failed files   :", sum(!final_ok), "\n")
cat("============================\n\n")

if (!all(final_ok)) {
  bad <- manifest[!final_ok, c("id", "cases", "file_name")]
  cat("The following files did not verify:\n")
  print(utils::head(bad, 20))
  if (nrow(bad) > 20) cat("... and", nrow(bad) - 20, "more.\n")
  stop("Download incomplete: ", sum(!final_ok), " of ", nrow(manifest),
       " files failed md5 verification.\n",
       "Re-run this script - verified files are kept and only these ",
       "will be retried.")
}
cat("All", nrow(manifest), "files present and md5-verified.\n\n")

if (DL_ONLY) {
  cat("--download-only given: stopping before assembly.\n")
  cat("Re-run without that flag to build the matrix.\n")
  quit(save = "no", status = 0)
}

# ---- 7. Rebuild the GDCquery object that GDCprepare() requires ----
# Rather than re-querying with barcode filters (which can return a
# different set than the manifest), query the full cohort once and then
# subset the results to exactly our manifest file_ids. That makes the
# query object provably identical to what we downloaded.
# GDCquery() needs a live connection to the GDC, and on a flaky link it
# fails with "SSL connect error" - which would be a nasty surprise after
# a multi-day download. So: retry it, and CACHE the successful object to
# disk. Later runs reuse the cache and need no network for this step.
QUERY_CACHE <- "data/processed/methylation_query_object.rds"

query <- NULL
if (file.exists(QUERY_CACHE)) {
  cat("Loading cached GDCquery object from", QUERY_CACHE, "\n")
  query <- tryCatch(readRDS(QUERY_CACHE), error = function(e) NULL)
  if (!is.null(query) && !all(manifest$id %in% query$results[[1]]$file_id)) {
    cat("Cached query does not cover the manifest - discarding it.\n")
    query <- NULL
  }
}

if (is.null(query)) {
  cat("Rebuilding GDCquery object (needs GDC connectivity)...\n")
  for (attempt in 1:5) {
    cat("  GDCquery attempt", attempt, "of 5...\n")
    query <- tryCatch(
      GDCquery(
        project       = "TCGA-BRCA",
        data.category = "DNA Methylation",
        platform      = "Illumina Human Methylation 450",
        data.type     = "Methylation Beta Value"
      ),
      error = function(e) {
        cat("  failed:", conditionMessage(e), "\n")
        NULL
      })
    if (!is.null(query)) break
    if (attempt < 5) { cat("  waiting 30s before retrying...\n"); Sys.sleep(30) }
  }
  if (is.null(query)) {
    stop("GDCquery() failed 5 times - the GDC is unreachable right now.\n",
         "Your DOWNLOADED FILES ARE SAFE and still verified; nothing is lost.\n",
         "Just re-run this script when the connection is better.")
  }
  saveRDS(query, QUERY_CACHE)
  cat("Cached the query object to", QUERY_CACHE,
      "- future runs will not need the network for this step.\n")
}

res <- query$results[[1]]
cat("Query covers", nrow(res), "files for the full cohort.\n")

keep <- match(manifest$id, res$file_id)
if (any(is.na(keep))) {
  stop(sum(is.na(keep)), " manifest file_id(s) were not found in the ",
       "fresh GDC query. GDC may have issued a new data release since ",
       "script 01 was run - re-run script 01 to refresh the manifest.")
}
query$results[[1]] <- res[keep, , drop = FALSE]

cat("Query subset to", nrow(query$results[[1]]), "files.\n")
stopifnot(identical(query$results[[1]]$file_id, manifest$id))
stopifnot(!any(duplicated(query$results[[1]]$cases)))
cat("Query matches the manifest exactly; no duplicated barcodes.\n\n")

# ---- 8. Assemble into a single object ----
cat("Assembling with GDCprepare() - reading and merging",
    nrow(manifest), "files. This takes a while and is memory-hungry.\n")
meth_data <- GDCprepare(query, directory = RAW_DIR, summarizedExperiment = TRUE)

cat("\n===== ASSEMBLED OBJECT =====\n")
cat("Class      :", paste(class(meth_data), collapse = ", "), "\n")
cat("Dimensions :", paste(dim(meth_data), collapse = " x "),
    "(probes x samples)\n")
cat("Probes     :", nrow(meth_data), "\n")
cat("Samples    :", ncol(meth_data), "\n")
if (ncol(meth_data) != nrow(manifest)) {
  cat("WARNING: sample count", ncol(meth_data),
      "does not equal manifest rows", nrow(manifest), "\n")
}
cat("============================\n\n")

# ---- 9. Save, then prove the saved file actually reloads ----
dir.create(dirname(OUT_RDS), recursive = TRUE, showWarnings = FALSE)
cat("Saving to", OUT_RDS, "...\n")
saveRDS(meth_data, OUT_RDS)
cat("Saved. File size:", round(file.size(OUT_RDS) / 1024^2, 1), "MB\n")

cat("Reloading the saved object to confirm it is not corrupt...\n")
rm(meth_data); invisible(gc())
check <- readRDS(OUT_RDS)
cat("Reloaded OK. Dimensions:", paste(dim(check), collapse = " x "), "\n")
reload_dims <- dim(check)
rm(check); invisible(gc())

# ---- 10. Optional: remove the bulky raw cache ----
# OFF by default. data/raw/ is ~8.3 GiB and re-downloading it over an
# external USB drive is expensive, so deletion must be asked for
# explicitly with --purge-raw once you are happy with the .rds.
if (PURGE_RAW) {
  if (file.exists(OUT_RDS) && reload_dims[2] == nrow(manifest)) {
    cat("\n--purge-raw given and .rds verified: deleting", RAW_DIR, "...\n")
    unlink(RAW_DIR, recursive = TRUE, force = TRUE)
    dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
    cat("Raw cache removed. Only the assembled .rds is kept.\n")
  } else {
    cat("\nRefusing to purge: the saved object did not verify as expected.\n")
  }
} else {
  cat("\ndata/raw/ kept (",
      round(sum(as.numeric(manifest$file_size)) / 1024^3, 2),
      " GiB). Re-run with --purge-raw to delete it once you are happy.\n",
      sep = "")
}

cat("\nScript 02 complete.\n")
