#!/usr/bin/env Rscript
# ============================================================
# 08b_train_epigenomic_baseline_models.R
#
# COPY OF SCRIPT 08 WITH ONLY THREE CHANGES (script 08 is not modified):
#   1. ROC curves: 3 separate figures, one per feature set (A, B+C, A+B+C),
#      each showing the curves of all 5 models (script 08 had one combined figure).
#   2. Metrics heatmap: Specificity and MCC columns removed.
#   3. Outputs are saved in separate "08b" folders, so script 08's results
#      are not overwritten:
#        outputs/figures/08b_models/, outputs/tables/08b_models/,
#        outputs/models/08b_models/, outputs/reports/08b_models_report.txt
#   Everything else (models, settings, tables, other figures) is identical.
#
# ---------------- original script 08 description ----------------
#
# Trains the Month 2 baseline machine-learning models and answers the
# central question: does adding DNA methylation improve 5-year survival
# prediction over Month 1's clinical + genomic features?
# (Internship brief: "16_train_epigenomic_baseline_model"; Week 4.)
#
# MODELS - the same five as revised Month 1, same settings:
#   1. Logistic Regression (glm)
#   2. Elastic Net        (glmnet, alpha = 0.5, lambda chosen by 5-fold CV, AUC)
#   3. Random Forest      (500 trees, mtry = sqrt(p), class weights)
#   4. XGBoost            (depth 4, eta 0.05, subsample / colsample 0.8,
#                          min_child_weight 5, class weight - see fix 2 below)
#   5. SVM                (radial kernel, cost 1, gamma 1/p, class weights)
#
# FEATURE SETS (built by script 07; same 147 training / 64 test patients):
#   A     clinical + genomic (34)  - the Month 1 baseline, re-trained here
#   B+C   methylation only   (29)
#   A+B+C combined           (63)
#
# TARGET: five_year_survival (revised Month 1): 1 = survived 5 years,
#   0 = died within 5 years. As in Month 1, the "positive" class is 1, so
#   Recall = survivors correctly identified and Specificity = non-survivors
#   correctly identified (the clinically important one). Threshold = 0.5.
#
# THREE DELIBERATE CHANGES FROM MONTH 1 (bug fixes):
#   1. XGBoost rounds: Month 1 chose the number of boosting rounds by early
#      stopping on the TEST set (evals = list(train, test)) and then scored
#      that same test set - a form of data leakage. Here the number of rounds
#      is chosen by 5-fold cross-validation on the TRAINING set only.
#   2. XGBoost class weight: Month 1 set scale_pos_weight = n_survived /
#      n_died, which up-weights the MAJORITY class (1 = survived) and makes the
#      imbalance worse. Here it is n_died / n_survived, which balances it
#      (the same direction Month 1 correctly used for Random Forest and SVM).
#   3. AUC direction is fixed (controls = 0, cases = 1). pROC's automatic
#      direction would silently turn a worse-than-chance model into a
#      better-than-chance one.
#
# ROBUSTNESS: the test set has only 18 deaths, so test AUCs are uncertain.
#   The script therefore also reports 95% confidence intervals, paired
#   DeLong tests between feature sets, and repeated 5-fold cross-validation
#   (5 repeats) on the training set.
#
# INPUT : data/processed/feature_table_epigenomic_ml.csv   (script 07)
#         outputs/models/feature_sets.csv                   (script 07)
#         outputs/models/train_test_split_ids.csv           (script 07)
# OUTPUT: outputs/models/08_trained_models.rds   - all 15 fitted models
#         outputs/tables/08_models/ - metrics, comparisons, importances,
#                                     test-set predictions
#         outputs/figures/08_models/ - 9 figures
#         outputs/reports/08_models_report.txt - everything printed
#
# Run from the Ubuntu terminal with:
#   Rscript scripts/08_train_epigenomic_baseline_models.R
# Test run whose outputs get a suffix (nothing real is overwritten):
#   Rscript scripts/08_train_epigenomic_baseline_models.R --tag=test
# ============================================================

# ---- Required packages (stop with a clear message if any is missing) ----
needed  <- c("data.table", "ggplot2", "pROC", "glmnet", "randomForest", "xgboost", "e1071")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R package(s): ", paste(missing, collapse = ", "), "\nInstall them once with:\n",
       "  Rscript -e 'install.packages(c(", paste0('"', missing, '"', collapse = ", "),
       "), repos = \"https://cloud.r-project.org\")'")
}
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(pROC)
  library(glmnet); library(randomForest); library(xgboost); library(e1071)
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

FEATURES  <- "data/processed/feature_table_epigenomic_ml.csv"
SETS_FILE <- "outputs/models/feature_sets.csv"
SPLIT     <- "outputs/models/train_test_split_ids.csv"
# 08b: separate output folders so script 08's results are not overwritten.
FIG_DIR   <- paste0("outputs/figures/08b_models", sfx)
TAB_DIR   <- paste0("outputs/tables/08b_models", sfx)
MOD_DIR   <- paste0("outputs/models/08b_models", sfx)
REPORT    <- paste0("outputs/reports/08b_models_report", sfx, ".txt")
for (d in c(FIG_DIR, TAB_DIR, MOD_DIR, dirname(REPORT))) dir.create(d, recursive = TRUE, showWarnings = FALSE)

SEED       <- 42
THRESHOLD  <- 0.5        # same as Month 1
CV_FOLDS   <- 5
CV_REPEATS <- 5
XGB_MAX_ROUNDS <- 300    # Month 1's nrounds; the actual number is chosen by training-set CV
MODELS <- c("Logistic Regression", "Elastic Net", "Random Forest", "XGBoost", "SVM")
SETS   <- c("A: Clinical + genomic", "B+C: Methylation", "A+B+C: Combined")
set.seed(SEED)

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
        strip.text = element_text(face = "bold"),
        plot.margin = ggplot2::margin(10, 15, 10, 10))   # ggplot2:: because randomForest also has a margin() function
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
                                                   formatC(p, format = "f", digits = 3)))

PAL_SET   <- setNames(c("#0F6E56", "#185FA5", "#8C510A"), SETS)
PAL_MODEL <- setNames(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"), MODELS)
PAL_BLOCK <- c("A: Clinical + genomic" = "#0F6E56", "B: Methylation summary" = "#185FA5",
               "C: Methylation PCs" = "#8C510A")

cat("Script 08 - baseline models: clinical+genomic vs methylation vs combined\n")
cat("Tag    :", ifelse(nzchar(TAG), TAG, "(none - real run)"), "\n")
cat("Packages: glmnet", as.character(packageVersion("glmnet")), "| randomForest", as.character(packageVersion("randomForest")),
    "| xgboost", as.character(packageVersion("xgboost")), "| e1071", as.character(packageVersion("e1071")),
    "| pROC", as.character(packageVersion("pROC")), "\n")
cat("Started:", format(run_start), "\n")

# ============================================================
# 1. LOAD THE FEATURE TABLE AND FEATURE SETS
# ============================================================
step("1. Load data")
for (f in c(FEATURES, SETS_FILE, SPLIT)) if (!file.exists(f)) stop("Input not found: ", f, "\nRun script 07 first.")
ft    <- fread(FEATURES)
fsets <- fread(SETS_FILE)
spl   <- fread(SPLIT)

stopifnot(!anyNA(ft), identical(sort(ft$patient_id), sort(spl$patient_id)),
          all(ft$split == spl$split[match(ft$patient_id, spl$patient_id)]))
forbidden <- c("os_event", "os_months", "OS_MONTHS", "OS_STATUS", "survival_days", "dfs_event", "DFS_MONTHS")
stopifnot(!any(fsets$feature %in% forbidden))

set_features <- list(fsets[set_A_clinical_genomic == TRUE, feature],
                     fsets[set_BC_methylation == TRUE, feature],
                     fsets[set_ABC_combined == TRUE, feature])
names(set_features) <- SETS
block_of <- setNames(fsets$block, fsets$feature)

train <- ft[split == "train"]; test <- ft[split == "test"]
y_tr <- train$five_year_survival; y_te <- test$five_year_survival
cat(sprintf("Training: %d patients (%d survived / %d died) | Test: %d patients (%d survived / %d died)\n",
            nrow(train), sum(y_tr == 1), sum(y_tr == 0), nrow(test), sum(y_te == 1), sum(y_te == 0)))
for (s in SETS) cat(sprintf("  %-24s %2d features\n", s, length(set_features[[s]])))
cat("Checks passed: no missing values; split identical to script 07; no survival columns among features.\n")

# ============================================================
# 2. MODEL FUNCTIONS (fit on training data -> P(survived) for new data)
# ============================================================
w0 <- function(y) sum(y == 1) / sum(y == 0)            # weight for the minority class (died), as in Month 1

fit_lr <- function(X, y) {
  warn <- character(0)
  m <- withCallingHandlers(glm(y ~ ., data = data.frame(y = y, X, check.names = FALSE), family = binomial()),
                           warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(model = m, warnings = unique(warn),
       predict = function(Xn) as.numeric(predict(m, newdata = data.frame(Xn, check.names = FALSE), type = "response")))
}
fit_en <- function(X, y) {
  set.seed(SEED)
  m <- cv.glmnet(as.matrix(X), y, alpha = 0.5, family = "binomial", nfolds = 5, type.measure = "auc")
  list(model = m, predict = function(Xn) as.numeric(predict(m, as.matrix(Xn), s = "lambda.min", type = "response")))
}
fit_rf <- function(X, y) {
  set.seed(SEED)
  m <- randomForest(x = X, y = factor(y, levels = c(0, 1)), ntree = 500, mtry = floor(sqrt(ncol(X))),
                    classwt = c("0" = w0(y), "1" = 1), importance = TRUE)
  list(model = m, predict = function(Xn) as.numeric(predict(m, Xn, type = "prob")[, "1"]))
}
# scale_pos_weight multiplies the weight of the POSITIVE class (1 = survived,
# the majority). Month 1 used n_survived / n_died (~2.6), which makes the
# imbalance WORSE - the reason its XGBoost predicted every patient as a
# survivor. The balancing value is n_died / n_survived (~0.39).
xgb_params <- function(y) list(objective = "binary:logistic", eval_metric = "auc", max_depth = 4, eta = 0.05,
                               subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 5,
                               scale_pos_weight = sum(y == 0) / sum(y == 1), nthread = 1)
# Number of boosting rounds chosen by 5-fold CV on the TRAINING data only.
xgb_best_rounds <- function(X, y) {
  set.seed(SEED)
  cv <- xgb.cv(params = xgb_params(y), data = xgb.DMatrix(as.matrix(X), label = y), nrounds = XGB_MAX_ROUNDS,
               nfold = 5, stratified = TRUE, early_stopping_rounds = 30, verbose = FALSE)
  bi <- cv$best_iteration
  if (is.null(bi)) bi <- cv$early_stop$best_iteration
  if (is.null(bi)) bi <- which.max(cv$evaluation_log[[grep("test.*auc.*mean", names(cv$evaluation_log))[1]]])
  as.integer(bi)
}
fit_xgb <- function(X, y, nrounds) {
  set.seed(SEED)
  m <- xgb.train(params = xgb_params(y), data = xgb.DMatrix(as.matrix(X), label = y), nrounds = nrounds, verbose = 0)
  list(model = m, nrounds = nrounds, predict = function(Xn) as.numeric(predict(m, xgb.DMatrix(as.matrix(Xn)))))
}
fit_svm <- function(X, y) {
  keep <- names(X)[vapply(X, function(z) length(unique(z)) > 1, logical(1))]   # SVM cannot scale constant columns
  set.seed(SEED)
  m <- svm(x = as.matrix(X[, ..keep]), y = factor(y, levels = c(0, 1)), kernel = "radial", probability = TRUE,
           class.weights = c("0" = w0(y), "1" = 1), cost = 1, gamma = 1 / length(keep))
  list(model = m, predict = function(Xn) as.numeric(attr(predict(m, as.matrix(Xn[, ..keep]), probability = TRUE),
                                                          "probabilities")[, "1"]))
}
fit_model <- function(model, X, y, xgb_rounds = NULL) {
  switch(model, "Logistic Regression" = fit_lr(X, y), "Elastic Net" = fit_en(X, y), "Random Forest" = fit_rf(X, y),
         "XGBoost" = fit_xgb(X, y, xgb_rounds), "SVM" = fit_svm(X, y))
}

# Month 1's metrics (positive class = 1) plus AUC 95% CI and balanced accuracy.
roc_fixed <- function(y, p) pROC::roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)
metrics <- function(y, p, thr = THRESHOLD) {
  pc <- as.integer(p >= thr)
  tp <- sum(y == 1 & pc == 1); tn <- sum(y == 0 & pc == 0); fp <- sum(y == 0 & pc == 1); fn <- sum(y == 1 & pc == 0)
  prec <- if ((tp + fp) == 0) 0 else tp / (tp + fp); rec <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
  spec <- if ((tn + fp) == 0) 0 else tn / (tn + fp); f1 <- if ((prec + rec) == 0) 0 else 2 * prec * rec / (prec + rec)
  md <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)); mcc <- if (md == 0) 0 else (tp * tn - fp * fn) / md
  r <- roc_fixed(y, p); ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
  data.table(AUC = as.numeric(pROC::auc(r)), AUC_low = ci[1], AUC_high = ci[3], Accuracy = (tp + tn) / length(y),
             Balanced_accuracy = (rec + spec) / 2, Recall = rec, Precision = prec, Specificity = spec,
             F1 = f1, MCC = mcc, TP = tp, TN = tn, FP = fp, FN = fn)
}

# ============================================================
# 3. TRAIN ON 70%, EVALUATE ON THE UNTOUCHED 30%
# ============================================================
step("3. Train and test (15 models)")
fits <- list(); res <- list(); rocs <- list(); preds <- data.table(patient_id = test$patient_id, truth = y_te)
xgb_rounds <- setNames(integer(length(SETS)), SETS)
lr_warn <- list()

for (s in SETS) {
  cols <- set_features[[s]]
  Xtr <- train[, ..cols]; Xte <- test[, ..cols]
  xgb_rounds[s] <- xgb_best_rounds(Xtr, y_tr)
  cat(sprintf("\n--- %s (%d features) | XGBoost rounds chosen by training CV: %d ---\n", s, length(cols), xgb_rounds[s]))
  for (mname in MODELS) {
    f <- fit_model(mname, Xtr, y_tr, xgb_rounds[s])
    p <- f$predict(Xte)
    key <- paste(mname, s, sep = " | ")
    fits[[key]] <- f; rocs[[key]] <- roc_fixed(y_te, p)
    preds[[key]] <- round(p, 5)
    res[[key]] <- cbind(data.table(Model = mname, Feature_set = s, n_features = length(cols)), metrics(y_te, p))
    if (mname == "Logistic Regression" && length(f$warnings)) lr_warn[[s]] <- f$warnings
  }
  show(rbindlist(res)[Feature_set == s, .(Model, AUC = round(AUC, 3), CI = sprintf("%.2f-%.2f", AUC_low, AUC_high),
                                            BalAcc = round(Balanced_accuracy, 3), Recall = round(Recall, 3),
                                            Specificity = round(Specificity, 3), F1 = round(F1, 3), MCC = round(MCC, 3))])
}
if (length(lr_warn)) {
  cat("\nLogistic Regression warnings (expected with many features and few deaths - see Elastic Net):\n")
  for (s in names(lr_warn)) cat("  ", s, ":", paste(lr_warn[[s]], collapse = " | "), "\n")
}

results <- rbindlist(res)
num_cols <- setdiff(names(results), c("Model", "Feature_set", "n_features", "TP", "TN", "FP", "FN"))
results[, (num_cols) := lapply(.SD, round, 4), .SDcols = num_cols]
save_table(results, "08_01_test_metrics_all_models.csv")
save_table(preds, "08_09_test_set_predictions.csv")

# ============================================================
# 4. DOES METHYLATION HELP? PAIRED COMPARISONS ON THE SAME TEST PATIENTS
# ============================================================
step("4. Compare feature sets (DeLong)")
cmp <- rbindlist(lapply(MODELS, function(mname) {
  rbindlist(lapply(c("A+B+C: Combined", "B+C: Methylation"), function(s2) {
    r1 <- rocs[[paste(mname, "A: Clinical + genomic", sep = " | ")]]; r2 <- rocs[[paste(mname, s2, sep = " | ")]]
    t <- tryCatch(pROC::roc.test(r1, r2, method = "delong", paired = TRUE), error = function(e) NULL)
    data.table(Model = mname, comparison = paste(s2, "vs A"), AUC_A = as.numeric(auc(r1)), AUC_other = as.numeric(auc(r2)),
               delta_AUC = as.numeric(auc(r2)) - as.numeric(auc(r1)), p_value = if (is.null(t)) NA_real_ else t$p.value)
  }))
}))
cmp[, (c("AUC_A", "AUC_other", "delta_AUC")) := lapply(.SD, round, 3), .SDcols = c("AUC_A", "AUC_other", "delta_AUC")]
cmp[, p_value := signif(p_value, 3)]
show(cmp)
cat(sprintf("(Only %d deaths in the test set, so only large AUC differences can reach p < 0.05.)\n", sum(y_te == 0)))
save_table(cmp, "08_02_feature_set_comparison_delong.csv")

# ============================================================
# 5. ROBUSTNESS: REPEATED 5-FOLD CROSS-VALIDATION ON THE TRAINING SET
#    Every model x feature set is refitted 25 times on 4/5 of the training
#    patients and scored on the remaining 1/5 - a much larger evaluation
#    than the 18 test deaths. (Note: script 07 fitted the methylation PCs on
#    all training patients, which makes the CV estimates for sets with PCs
#    slightly optimistic; the test-set results above are not affected.)
# ============================================================
step(sprintf("5. Cross-validation (%d x %d-fold)", CV_REPEATS, CV_FOLDS))
strat_folds <- function(y, k, seed) {
  set.seed(seed); f <- integer(length(y))
  for (cl in unique(y)) { i <- which(y == cl); f[i] <- sample(rep(seq_len(k), length.out = length(i))) }
  f
}
cv <- list(); t_cv <- Sys.time()
for (r in seq_len(CV_REPEATS)) {
  folds <- strat_folds(y_tr, CV_FOLDS, SEED + r)
  for (k in seq_len(CV_FOLDS)) {
    tr_i <- which(folds != k); va_i <- which(folds == k)
    for (s in SETS) {
      cols <- set_features[[s]]
      Xa <- train[tr_i, ..cols]; Xv <- train[va_i, ..cols]
      for (mname in MODELS) {
        p <- tryCatch(suppressWarnings(fit_model(mname, Xa, y_tr[tr_i], xgb_rounds[s])$predict(Xv)), error = function(e) NULL)
        auc_v <- if (is.null(p)) NA_real_ else as.numeric(auc(roc_fixed(y_tr[va_i], p)))
        cv[[length(cv) + 1]] <- data.table(repeat_no = r, fold = k, Model = mname, Feature_set = s, AUC = auc_v)
      }
    }
  }
  cat(sprintf("  repeat %d of %d done (%s)\n", r, CV_REPEATS, fmt_time(as.numeric(difftime(Sys.time(), t_cv, units = "secs")))))
}
cv <- rbindlist(cv)
cv_sum <- cv[, .(cv_AUC_mean = round(mean(AUC, na.rm = TRUE), 3), cv_AUC_sd = round(sd(AUC, na.rm = TRUE), 3),
                 failed_fits = sum(is.na(AUC))), by = .(Model, Feature_set)]
# Paired fold-by-fold gain from adding methylation to set A.
gain <- merge(cv[Feature_set == "A: Clinical + genomic", .(repeat_no, fold, Model, AUC_A = AUC)],
              cv[Feature_set == "A+B+C: Combined", .(repeat_no, fold, Model, AUC_ABC = AUC)], by = c("repeat_no", "fold", "Model"))
gain_sum <- gain[, .(mean_gain = round(mean(AUC_ABC - AUC_A, na.rm = TRUE), 3),
                     pct_folds_improved = round(100 * mean(AUC_ABC > AUC_A, na.rm = TRUE), 0)), by = Model]
cat("\nCross-validated AUC (mean +/- SD over 25 folds):\n")
show(dcast(cv_sum, Model ~ Feature_set, value.var = "cv_AUC_mean"))
cat("\nAdding methylation to clinical + genomic (A+B+C minus A), fold by fold:\n"); show(gain_sum)
save_table(cv_sum, "08_03_cross_validation_auc.csv")
save_table(gain_sum, "08_04_cv_gain_from_methylation.csv")

# ============================================================
# 6. WHICH FEATURES DO THE MODELS USE? (combined set A+B+C)
# ============================================================
step("6. Feature importance (combined set)")
ABC <- "A+B+C: Combined"
# Elastic Net: non-zero coefficients (all three sets).
en_coef <- rbindlist(lapply(SETS, function(s) {
  cf <- as.matrix(coef(fits[[paste("Elastic Net", s, sep = " | ")]]$model, s = "lambda.min"))
  d <- data.table(feature = rownames(cf), coefficient = cf[, 1])[feature != "(Intercept)" & coefficient != 0]
  d[, Feature_set := s]; d[, block := block_of[feature]]; d
}))
cat("Elastic Net - features kept (non-zero coefficient) per set:\n"); show(en_coef[, .N, by = Feature_set])
save_table(en_coef[order(Feature_set, -abs(coefficient))], "08_05_elastic_net_coefficients.csv")

rf_imp <- {
  im <- importance(fits[[paste("Random Forest", ABC, sep = " | ")]]$model, type = 1)
  data.table(feature = rownames(im), importance = im[, 1], block = block_of[rownames(im)])[order(-importance)]
}
xgb_imp <- {
  im <- as.data.table(xgb.importance(model = fits[[paste("XGBoost", ABC, sep = " | ")]]$model))
  im[, .(feature = Feature, importance = Gain, block = block_of[Feature])][order(-importance)]
}
save_table(rf_imp, "08_06_random_forest_importance_combined.csv")
save_table(xgb_imp, "08_07_xgboost_gain_combined.csv")
cat("\nTop 10 Random Forest features (combined set):\n"); show(head(rf_imp[, .(feature, block, importance = round(importance, 2))], 10))
cat("\nTop 10 XGBoost features by gain (combined set):\n"); show(head(xgb_imp[, .(feature, block, importance = round(importance, 3))], 10))

# Share of each model's importance coming from each block.
share <- rbind(
  rf_imp[importance > 0, .(importance = sum(importance)), by = block][, method := "Random Forest (accuracy decrease)"],
  xgb_imp[, .(importance = sum(importance)), by = block][, method := "XGBoost (gain)"],
  en_coef[Feature_set == ABC, .(importance = sum(abs(coefficient))), by = block][, method := "Elastic Net (|coefficient|)"])
share[, pct := round(100 * importance / sum(importance), 1), by = method]
cat("\nShare of importance by feature block (combined set):\n"); show(dcast(share, method ~ block, value.var = "pct", fill = 0))
save_table(share, "08_08_importance_share_by_block.csv")

# ============================================================
# 7. FIGURES
# ============================================================
step("7. Figures")
# Fig 1 (08b change): ROC curves - one figure per feature set, one line per model.
roc_df <- rbindlist(lapply(names(rocs), function(k) {
  r <- rocs[[k]]; parts <- strsplit(k, " | ", fixed = TRUE)[[1]]
  data.table(Model = parts[1], Feature_set = parts[2], fpr = 1 - r$specificities, tpr = r$sensitivities)[order(fpr, tpr)]
}))
roc_file <- c("A: Clinical + genomic" = "08_01a_roc_curves_A_clinical_genomic.png",
              "B+C: Methylation"      = "08_01b_roc_curves_BC_methylation.png",
              "A+B+C: Combined"       = "08_01c_roc_curves_ABC_combined.png")
for (s in SETS) {
  d <- roc_df[Feature_set == s]
  d[, Model := factor(Model, levels = MODELS)]
  # AUC values go in the legend, where they cannot collide with the curves.
  auc_lab <- sapply(MODELS, function(m) sprintf("%s (AUC %.2f)", m, results[Model == m & Feature_set == s, AUC]))
  p <- ggplot(d, aes(fpr, tpr, colour = Model)) +
    geom_abline(linetype = "dashed", colour = "grey60") +
    geom_step(linewidth = 0.9) +
    scale_colour_manual(values = PAL_MODEL, labels = auc_lab, name = NULL) +
    guides(colour = guide_legend(ncol = 2)) +
    coord_equal() +
    labs(title = sprintf("ROC curves - %s (%d features)", s, length(set_features[[s]])),
         subtitle = sprintf("Untouched test set (%d patients, %d deaths)\nPositive class = survived 5 years; dashed line = chance",
                            length(y_te), sum(y_te == 0)),
         x = "False positive rate (1 - specificity)", y = "True positive rate (recall)") + theme_qc
  save_plot(p, roc_file[[s]], 8.5, 9)
}

# Fig 2: test AUC with 95% CI - the headline comparison.
rd <- copy(results); rd[, Model := factor(Model, levels = rev(MODELS))]; rd[, Feature_set := factor(Feature_set, levels = SETS)]
p <- ggplot(rd, aes(AUC, Model, colour = Feature_set)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = AUC_low, xmax = AUC_high), width = 0.25, orientation = "y",
                position = position_dodge(width = 0.6)) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = PAL_SET, name = NULL) +
  scale_x_continuous(limits = c(min(0.3, min(rd$AUC_low)), 1), breaks = seq(0.3, 1, 0.1)) +
  labs(title = "Test-set AUC with 95% confidence interval",
       subtitle = "Wide intervals = only 18 test deaths; overlapping intervals mean the sets are not clearly different",
       x = "AUC (0.5 = chance, 1 = perfect)", y = NULL) + theme_qc
save_plot(p, "08_02_test_auc_with_ci.png", 10, 6)

# Fig 3: cross-validated AUC distribution.
cvp <- copy(cv[!is.na(AUC)]); cvp[, Model := factor(Model, levels = MODELS)]; cvp[, Feature_set := factor(Feature_set, levels = SETS)]
p <- ggplot(cvp, aes(Model, AUC, fill = Feature_set)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
  geom_boxplot(outlier.size = 0.7, width = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = PAL_SET, name = NULL) +
  labs(title = sprintf("Cross-validated AUC on the training set (%d repeats x %d folds)", CV_REPEATS, CV_FOLDS),
       subtitle = "Each box = 25 held-out folds; a larger, steadier estimate than the single test split",
       x = NULL, y = "AUC in held-out fold") + theme_qc
save_plot(p, "08_03_cross_validation_auc.png", 11, 6)

# Fig 4: metrics as a heatmap (08b change: Specificity and MCC removed).
mh <- melt(results[, .(Model, Feature_set, AUC, Balanced_accuracy, Recall, F1)],
           id.vars = c("Model", "Feature_set"), variable.name = "metric", value.name = "value")
mh[, row := factor(paste(Model, "|", sub(":.*", "", Feature_set)),
                   levels = rev(paste(rep(MODELS, each = 3), "|", sub(":.*", "", rep(SETS, times = 5)))))]
mh[, metric := factor(gsub("_", " ", metric), levels = gsub("_", " ", c("AUC", "Balanced_accuracy", "Recall", "F1")))]
p <- ggplot(mh, aes(metric, row, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", value), colour = value > 0.8), size = 3.2) +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +
  # One white-to-blue scale for every metric (darker = higher); negative MCC shown as white.
  scale_fill_gradient(low = "white", high = "#185FA5", limits = c(0, 1), oob = scales::squish, name = NULL) +
  labs(title = "All test-set metrics (threshold 0.5)",
       subtitle = "Rows = model | feature set (A = clinical+genomic, B+C = methylation, A+B+C = combined)",
       x = NULL, y = NULL) +
  theme_qc + theme(panel.grid = element_blank(), legend.position = "right")
save_plot(p, "08_04_metrics_heatmap.png", 10, 8)

# Fig 5: confusion matrices of the five models on the combined set.
cm <- results[Feature_set == ABC, .(Model, TP, FN, FP, TN)]
cm <- rbindlist(lapply(seq_len(nrow(cm)), function(i) data.table(Model = cm$Model[i],
  truth = c("Survived", "Survived", "Died", "Died"), predicted = c("Survived", "Died", "Survived", "Died"),
  n = c(cm$TP[i], cm$FN[i], cm$FP[i], cm$TN[i]))))
cm[, correct := truth == predicted]; cm[, Model := factor(Model, levels = MODELS)]
cm[, truth := factor(truth, levels = c("Died", "Survived"))]; cm[, predicted := factor(predicted, levels = c("Survived", "Died"))]
p <- ggplot(cm, aes(predicted, truth, fill = correct)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = n), size = 5, fontface = "bold", colour = "white") +
  facet_wrap(~ Model, nrow = 1) +
  scale_fill_manual(values = c(`TRUE` = "#185FA5", `FALSE` = "#EF4444"), labels = c(`TRUE` = "Correct", `FALSE` = "Wrong"), name = NULL) +
  labs(title = "Confusion matrices - combined feature set (test set, threshold 0.5)",
       subtitle = sprintf("Test set: %d survived, %d died; the bottom-right cell of each panel = deaths correctly flagged",
                          sum(y_te == 1), sum(y_te == 0)),
       x = "Predicted", y = "Actual") + theme_qc + theme(panel.grid = element_blank())
save_plot(p, "08_05_confusion_matrices_combined.png", 13, 4.5)

# Figs 6-8: importance of features in the combined models, coloured by block.
imp_plot <- function(d, title, xlab, n = 20) {
  d <- head(d[order(-abs(importance))], n); d[, feature := factor(feature, levels = rev(feature))]
  ggplot(d, aes(importance, feature, fill = block)) +
    geom_col(width = 0.7) + geom_vline(xintercept = 0, colour = "grey40") +
    scale_fill_manual(values = PAL_BLOCK, name = NULL, drop = FALSE) +
    guides(fill = guide_legend(nrow = 1)) +
    labs(title = title, x = xlab, y = NULL) + theme_qc
}
en_abc <- en_coef[Feature_set == ABC, .(feature, importance = coefficient, block)]
if (nrow(en_abc)) {
  p <- imp_plot(en_abc, sprintf("Elastic Net - %d features kept (combined set)", nrow(en_abc)),
                "Coefficient (positive = more likely to survive 5 years)") +
    labs(subtitle = "Features with a zero coefficient were dropped by the penalty")
  save_plot(p, "08_06_elastic_net_coefficients.png", 10, 7)
}
p <- imp_plot(rf_imp, "Random Forest - top 20 features (combined set)", "Mean decrease in accuracy when the feature is scrambled")
save_plot(p, "08_07_random_forest_importance.png", 10, 7.5)
p <- imp_plot(xgb_imp, "XGBoost - top 20 features (combined set)", "Gain (share of the improvement the feature provides)")
save_plot(p, "08_08_xgboost_importance.png", 10, 7.5)

# Fig 9: how much of each model's importance comes from each block.
sh <- copy(share); sh[, block := factor(block, levels = names(PAL_BLOCK))]
p <- ggplot(sh, aes(pct, method, fill = block)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = ifelse(pct >= 6, paste0(pct, "%"), "")), position = position_stack(vjust = 0.5),
            colour = "white", fontface = "bold", size = 3.6) +
  scale_fill_manual(values = PAL_BLOCK, name = NULL) +
  labs(title = "Where the combined models get their information",
       subtitle = "Share of total importance by feature block (combined set)", x = "% of importance", y = NULL) + theme_qc
save_plot(p, "08_09_importance_share_by_block.png", 10, 4.5)

# ============================================================
# 8. SAVE MODELS AND A MONTH 1 REFERENCE TABLE
# ============================================================
step("8. Save models")
saveRDS(list(models = lapply(fits, `[[`, "model"), xgb_rounds = xgb_rounds, feature_sets = set_features,
             threshold = THRESHOLD, seed = SEED, created = format(Sys.time())),
        file.path(MOD_DIR, "08_trained_models.rds"))
record(file.path(MOD_DIR, "08_trained_models.rds"), "model object")
cat("  saved models:", file.path(MOD_DIR, "08_trained_models.rds"), "\n")

# Revised Month 1 results are shown for context only: different cohort (305
# patients, no methylation requirement), 80/20 split, 62 test patients.
m1_ref <- data.table(Model = c("Elastic Net", "XGBoost", "SVM", "Logistic Regression", "Random Forest"),
                     Month1_AUC_305pts = c(0.815, 0.792, 0.780, 0.726, 0.725))
m2A <- results[Feature_set == "A: Clinical + genomic", .(Model, Month2_setA_AUC = round(AUC, 3))]
m2C <- results[Feature_set == ABC, .(Model, Month2_combined_AUC = round(AUC, 3))]
ref <- Reduce(function(a, b) merge(a, b, by = "Model"), list(m1_ref, m2A, m2C))
ref[, note := ifelse(Model == "XGBoost", "Month 1: rounds tuned on its test set (optimistic) and inverted class weight", "")]
cat("\nContext - revised Month 1 (different cohort, not directly comparable) vs Month 2:\n"); show(ref)
save_table(ref, "08_10_month1_reference.csv")

# ============================================================
# 9. KEY INSIGHTS AND OUTPUT INDEX
# ============================================================
step("9. Key insights")
short <- c("Logistic Regression" = "LR", "Elastic Net" = "EN", "Random Forest" = "RF", "XGBoost" = "XGB", "SVM" = "SVM")
best <- results[, .SD[which.max(AUC)], by = Feature_set]
cv_best <- cv_sum[, .SD[which.max(cv_AUC_mean)], by = Feature_set]
sig <- cmp[comparison == "A+B+C: Combined vs A" & !is.na(p_value) & p_value < 0.05]
cat("\n====================== KEY INSIGHTS ======================\n")
for (i in seq_len(nrow(best)))
  cat(sprintf("* Best on the test set for %-22s: %s, AUC %.2f (95%% CI %.2f-%.2f); best in CV: %s, %.2f.\n",
              best$Feature_set[i], best$Model[i], best$AUC[i], best$AUC_low[i], best$AUC_high[i],
              cv_best[Feature_set == best$Feature_set[i], Model], cv_best[Feature_set == best$Feature_set[i], cv_AUC_mean]))
cat(sprintf("* Adding methylation to clinical + genomic changed test AUC by: %s; significant (DeLong p<0.05) for %d of 5 models.\n",
            paste(sprintf("%s %+.3f", short[cmp[comparison == "A+B+C: Combined vs A", Model]],
                          cmp[comparison == "A+B+C: Combined vs A", delta_AUC]), collapse = ", "), nrow(sig)))
cat(sprintf("* Cross-validation, adding methylation (share of 25 folds improved, mean AUC change): %s.\n",
            paste(sprintf("%s %d%% (%+.3f)", short[gain_sum$Model], gain_sum$pct_folds_improved, gain_sum$mean_gain), collapse = ", ")))
cat("  (LR = logistic regression, EN = elastic net, RF = random forest, XGB = XGBoost)\n")
cat(sprintf("* Top Random Forest feature: %s; top XGBoost feature: %s.\n", rf_imp$feature[1], xgb_imp$feature[1]))
cat(sprintf("* Caution: only %d deaths in the test set - AUC intervals are wide; treat differences as indicative.\n", sum(y_te == 0)))
cat("* Month 1 XGBoost had two issues (rounds tuned on its test set; class weight inverted) - both fixed here; state this in the report.\n")
cat("==========================================================\n")

idx <- do.call(rbind, outputs_written)
idx$size_kb <- round(file.size(idx$file) / 1024, 1)
save_table(idx, "08_99_outputs_index.csv")
cat("\nFiles written:", nrow(idx), "(", sum(idx$type == "figure"), "figures,", sum(idx$type == "table"), "tables )\n")
cat("Total run time:", fmt_time(as.numeric(difftime(Sys.time(), run_start, units = "secs"))), "\n")
cat("Report saved to:", REPORT, "\n")
cat("\nScript 08b complete.\n")

sink()
close(report_con)
