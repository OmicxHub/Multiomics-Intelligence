# ============================================================
# Script 08: ML Models + Cox Regression
# Project : Clinical-Genomics Harmonisation Pipeline
# Author  : Muteeba Azhar | MOI Internship | July 2026
#
# PURPOSE:
#   SECTION A — Binary ML Classification (305 patients)
#     Train 5 models on confirmed-outcome cohort.
#     Evaluate each on held-out test set.
#     Print individual plots (ROC, CM, importance) + compare.
#
#   SECTION B — Cox Proportional Hazards (964 patients)
#     Survival regression using all patients incl censored.
#     Censoring handled mathematically — no exclusion.
#     Outputs: C-index, hazard ratios, KM curves.
#
# DATA LEAKAGE PREVENTION:
#   1. Split 80/20 FIRST — then scale using train only
#   2. Scaler (mean/SD) fitted on training set
#   3. Same scaler applied to test set
#   4. Gene flags and one-hot columns not scaled (binary)
#
# MODELS:
#   ML-1 : Logistic Regression (glm — interpretable)
#   ML-2 : Elastic Net (glmnet alpha=0.5 — penalised LR)
#   ML-3 : Random Forest (randomForest — ensemble trees)
#   ML-4 : XGBoost (xgboost — gradient boosting)
#   ML-5 : Support Vector Machine (e1071 — kernel SVM)
#   Cox   : Cox Proportional Hazards (survival::coxph)
#
# EACH ML MODEL OUTPUTS:
#   - Performance metrics (AUC, Accuracy, Recall,
#     Precision, F1, Specificity, MCC)
#   - ROC curve plot (printed + saved)
#   - Confusion matrix plot (printed + saved)
#   - Feature importance / coefficients (printed + saved)
#
# FINAL COMPARISON:
#   - ROC overlay (all 5 ML models)
#   - Metrics bar chart (all models)
#   - Model performance summary table
#
# INPUTS:
#   data/processed/TCGA-BRCA/feature_table_ml.csv  (305 pts)
#   data/processed/TCGA-BRCA/feature_table_cox.csv (964 pts)
#
# OUTPUTS:
#   outputs/TCGA-BRCA/models/  — 5 model .rds files
#   outputs/TCGA-BRCA/models/train_test_split_ids.csv
#   outputs/TCGA-BRCA/tables/08_model_performance_all.csv
#   outputs/TCGA-BRCA/tables/08_cox_hazard_ratios.csv
#   outputs/TCGA-BRCA/figures/fig08_[model]_roc.png
#   outputs/TCGA-BRCA/figures/fig08_[model]_confusion.png
#   outputs/TCGA-BRCA/figures/fig08_[model]_importance.png
#   outputs/TCGA-BRCA/figures/fig08_roc_overlay.png
#   outputs/TCGA-BRCA/figures/fig08_model_comparison.png
#   outputs/TCGA-BRCA/figures/fig08_cox_forest_plot.png
#   outputs/TCGA-BRCA/figures/fig08_km_curves.png
# ============================================================

#packageVersion("rlang")
#packageVersion("vctrs")

library(tidyverse)
library(pROC)
library(randomForest)
library(xgboost)
library(glmnet)
library(survival)
library(survminer)
library(broom)
library(scales)



# Install SVM package if needed
if (!requireNamespace("e1071", quietly=TRUE))
  install.packages("e1071")
library(e1071)

setwd("E:/R_Studio/Projects/MOI/01R_clinical_genomics")

PATHS <- list(
  proc_brca        = "data/processed/TCGA-BRCA",
  out_brca_tables  = "outputs/TCGA-BRCA/tables",
  out_brca_figures = "outputs/TCGA-BRCA/figures",
  out_brca_models  = "outputs/TCGA-BRCA/models"
)

set.seed(42)  # Fixed for reproducibility

cat("====================================================\n")
cat(" SCRIPT 08: ML Models + Cox Regression\n")
cat(" SECTION A: 5 ML models on 305 confirmed patients\n")
cat(" SECTION B: Cox regression on all 964 patients\n")
cat("====================================================\n\n")


# ═══════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════

# Compute all classification metrics from predictions
compute_metrics <- function(actual, predicted_prob,
                             threshold=0.5, model_name="") {
  pred_class <- as.integer(predicted_prob >= threshold)
  tp <- sum(actual==1 & pred_class==1)
  tn <- sum(actual==0 & pred_class==0)
  fp <- sum(actual==0 & pred_class==1)
  fn <- sum(actual==1 & pred_class==0)

  acc   <- (tp+tn)/length(actual)
  prec  <- if((tp+fp)==0) 0 else tp/(tp+fp)
  rec   <- if((tp+fn)==0) 0 else tp/(tp+fn)
  spec  <- if((tn+fp)==0) 0 else tn/(tn+fp)
  f1    <- if((prec+rec)==0) 0 else 2*prec*rec/(prec+rec)
  mcc_d <- sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  mcc   <- if(mcc_d==0) 0 else (tp*tn-fp*fn)/mcc_d
  roc_obj <- pROC::roc(actual, predicted_prob, quiet=TRUE)
  auc   <- as.numeric(pROC::auc(roc_obj))

  list(model=model_name, threshold=threshold,
       AUC=round(auc,4), Accuracy=round(acc,4),
       Recall=round(rec,4), Precision=round(prec,4),
       Specificity=round(spec,4), F1=round(f1,4),
       MCC=round(mcc,4), TP=tp, TN=tn, FP=fp, FN=fn,
       roc_obj=roc_obj)
}

# Print confusion matrix to console
print_cm <- function(m) {
  cat(sprintf("  Confusion Matrix (threshold=%.2f):\n",
              m$threshold))
  cat("                 Predicted\n")
  cat("              Died  Survived\n")
  cat(sprintf("  Actual Died  | %3d  | %3d  |\n",m$TP,m$FN))
  cat(sprintf("  Actual Surv  | %3d  | %3d  |\n\n",m$FP,m$TN))
  cat(sprintf("  AUC=%.4f | Acc=%.4f | Recall=%.4f\n",
              m$AUC, m$Accuracy, m$Recall))
  cat(sprintf("  Spec=%.4f | F1=%.4f | MCC=%.4f\n\n",
              m$Specificity, m$F1, m$MCC))
}

# Plot and save ROC curve
plot_roc <- function(roc_obj, model_name, auc_val, path) {
  roc_df <- tibble(
    fpr=1-roc_obj$specificities,
    tpr=roc_obj$sensitivities)
  p <- ggplot(roc_df,aes(x=fpr,y=tpr)) +
    geom_ribbon(aes(ymin=0,ymax=tpr),
                fill="#185FA5",alpha=0.08) +
    geom_line(colour="#185FA5",linewidth=1.4) +
    geom_abline(linetype="dashed",colour="grey50") +
    annotate("text",x=0.6,y=0.15,
             label=sprintf("AUC = %.4f",auc_val),
             size=5,fontface="bold",colour="#185FA5") +
    scale_x_continuous(labels=percent_format(),
                       limits=c(0,1)) +
    scale_y_continuous(labels=percent_format(),
                       limits=c(0,1)) +
    labs(title=sprintf("%s — ROC Curve",model_name),
         subtitle=sprintf("Test set n=%d",
                          length(roc_obj$response)),
         x="False Positive Rate (1-Specificity)",
         y="True Positive Rate (Sensitivity)",
         caption="TCGA-BRCA | 5-year survival prediction") +
    theme_minimal(base_size=13) +
    theme(plot.background=element_rect(fill="white",
                                       colour=NA))
  print(p)
  ggsave(path,p,width=7,height=6,dpi=300)
}

# Plot and save confusion matrix
plot_cm <- function(m, model_name, path) {
  cm_df <- tibble(
    Actual    =c("Died","Died","Survived","Survived"),
    Predicted =c("Died","Survived","Died","Survived"),
    N   =c(m$TP,m$FN,m$FP,m$TN),
    pct =round(c(m$TP,m$FN,m$FP,m$TN)/
                 (m$TP+m$FN+m$FP+m$TN)*100,1),
    correct=c(TRUE,FALSE,FALSE,TRUE))
  p <- ggplot(cm_df,
              aes(x=Predicted,y=Actual,fill=correct)) +
    geom_tile(colour="white",linewidth=2) +
    geom_text(aes(label=sprintf("%d\n(%.1f%%)",N,pct)),
              size=7,fontface="bold",colour="white") +
    scale_fill_manual(
      values=c("TRUE"="#185FA5","FALSE"="#EF4444")) +
    scale_x_discrete(position="top") +
    labs(
      title   =sprintf("%s — Confusion Matrix",model_name),
      subtitle=sprintf(
        "Recall=%.3f | Spec=%.3f | F1=%.3f | MCC=%.3f",
        m$Recall,m$Specificity,m$F1,m$MCC),
      x="Predicted Class",y="Actual Class") +
    theme_minimal(base_size=13) +
    theme(plot.background=element_rect(fill="white",NA),
          legend.position="none",
          panel.grid=element_blank())
  print(p)
  ggsave(path,p,width=6,height=5,dpi=300)
}


# ═══════════════════════════════════════════════════════════════
# LOAD DATA
# ═══════════════════════════════════════════════════════════════

cat("LOADING DATA\n")
cat("─────────────\n")

feat_ml <- read_csv(
  file.path(PATHS$proc_brca,"feature_table_ml.csv"),
  show_col_types=FALSE)

feat_cox <- read_csv(
  file.path(PATHS$proc_brca,"feature_table_cox.csv"),
  show_col_types=FALSE)

cat(sprintf("ML dataset  : %d patients, %d columns\n",
            nrow(feat_ml), ncol(feat_ml)))
cat(sprintf("Cox dataset : %d patients, %d columns\n",
            nrow(feat_cox), ncol(feat_cox)))

cat("\n--- ML dataset: str() ---\n")
str(feat_ml)

cat("\n--- ML dataset: head(3) ---\n")
print(head(feat_ml,3))

# Feature columns (no patient_id, no target)
FEATURE_COLS <- setdiff(colnames(feat_ml),
                        c("patient_id","five_year_survival"))

# Continuous features (will be scaled)
CONT_COLS <- c("age","stage_numeric","tmb_nonsynonymous",
               "n_truncating","n_nonsense","n_frameshift_del",
               "n_frameshift_ins","n_splice_site")

cat(sprintf("\nFeature columns: %d\n",length(FEATURE_COLS)))
cat(sprintf("Continuous (will scale): %d\n",length(CONT_COLS)))
cat(sprintf("Binary (not scaled): %d\n",
            length(FEATURE_COLS)-length(CONT_COLS)))

# ══════════════════════════════════════════════════════════════
# ███████████████  SECTION A — ML MODELS  ████████████████████
# ══════════════════════════════════════════════════════════════

cat("\n")
cat("╔══════════════════════════════════════════════╗\n")
cat("║        SECTION A — ML MODELS                ║\n")
cat("║  5 models | 305 patients | 80/20 split      ║\n")
cat("╚══════════════════════════════════════════════╝\n\n")


# ── A1: STRATIFIED 80/20 TRAIN/TEST SPLIT ─────────────────────
cat("A1: Stratified 80/20 train/test split\n")
cat("──────────────────────────────────────\n")

# Split FIRST — scaling happens AFTER
set.seed(42)
survived_ids <- feat_ml$patient_id[feat_ml$five_year_survival==1]
died_ids     <- feat_ml$patient_id[feat_ml$five_year_survival==0]

train_surv <- sample(survived_ids, floor(0.8*length(survived_ids)))
train_died <- sample(died_ids,     floor(0.8*length(died_ids)))
train_ids  <- c(train_surv, train_died)
test_ids   <- feat_ml$patient_id[!feat_ml$patient_id %in% train_ids]

train_raw <- feat_ml |> filter(patient_id %in% train_ids)
test_raw  <- feat_ml |> filter(!patient_id %in% train_ids)

cat(sprintf("Train: %d patients (Survived=%d, Died=%d)\n",
            nrow(train_raw),
            sum(train_raw$five_year_survival==1),
            sum(train_raw$five_year_survival==0)))
cat(sprintf("Test : %d patients (Survived=%d, Died=%d)\n",
            nrow(test_raw),
            sum(test_raw$five_year_survival==1),
            sum(test_raw$five_year_survival==0)))

# Save split IDs for reproducibility
write_csv(
  tibble(patient_id=c(train_ids,test_ids),
         split=c(rep("train",length(train_ids)),
                 rep("test",length(test_ids)))),
  file.path(PATHS$out_brca_models,"train_test_split_ids.csv"))
cat("Saved: train_test_split_ids.csv\n")


# ── A2: SCALE CONTINUOUS FEATURES (train fit → both apply) ────
cat("\nA2: Feature scaling (train statistics only)\n")
cat("────────────────────────────────────────────\n")
cat("Fitting scaler on TRAINING set only — no data leakage.\n")

# Fit on training
col_means <- colMeans(train_raw[,CONT_COLS], na.rm=TRUE)
col_sds   <- apply(train_raw[,CONT_COLS], 2, sd, na.rm=TRUE)
col_sds[col_sds==0] <- 1  # avoid /0

cat("Training set means (continuous features):\n")
print(round(col_means,3))
cat("\nTraining set SDs:\n")
print(round(col_sds,3))

# Scale function
scale_df <- function(df, means, sds, cont_cols) {
  df_sc <- df
  for (col in cont_cols) {
    df_sc[[col]] <- (df[[col]] - means[col]) / sds[col]
  }
  df_sc
}

train_sc <- scale_df(train_raw, col_means, col_sds, CONT_COLS)
test_sc  <- scale_df(test_raw,  col_means, col_sds, CONT_COLS)

# Impute any remaining NA with 0 (post-scaling = mean of train)
for (col in FEATURE_COLS) {
  na_tr <- is.na(train_sc[[col]])
  na_te <- is.na(test_sc[[col]])
  if (any(na_tr)) {
    train_sc[[col]][na_tr] <- 0
    cat(sprintf("Imputed %d NA in train %s\n", sum(na_tr), col))
  }
  if (any(na_te)) {
    test_sc[[col]][na_te] <- 0
    cat(sprintf("Imputed %d NA in test %s\n", sum(na_te), col))
  }
}

# Matrices for model training
X_train <- as.matrix(train_sc[,FEATURE_COLS])
X_test  <- as.matrix(test_sc[,FEATURE_COLS])
y_train <- train_sc$five_year_survival
y_test  <- test_sc$five_year_survival

cat(sprintf("\nX_train: %d x %d\n",nrow(X_train),ncol(X_train)))
cat(sprintf("X_test : %d x %d\n",nrow(X_test),ncol(X_test)))

# Save scaler for future use
saveRDS(list(means=col_means, sds=col_sds,
             cont_cols=CONT_COLS,
             feature_cols=FEATURE_COLS),
        file.path(PATHS$out_brca_models,"feature_scaler.rds"))
cat("Saved: feature_scaler.rds\n")


# Store all model results
all_metrics <- list()
roc_objects <- list()


# ═══════════════════════════════════════════════════════════════
# ML-1: LOGISTIC REGRESSION
# ═══════════════════════════════════════════════════════════════

cat("\n──────────────────────────────────────────────\n")
cat(" ML-1: Logistic Regression\n")
cat("──────────────────────────────────────────────\n")

lr_df   <- data.frame(X_train, y=factor(y_train))
lr_model <- glm(y ~ ., data=lr_df, family=binomial())

cat("--- LR model summary ---\n")
print(summary(lr_model))

# Predict on test
lr_prob <- predict(lr_model,
                   newdata=data.frame(X_test),
                   type="response")

m_lr <- compute_metrics(y_test, lr_prob, 0.5, "Logistic Regression")
cat("\n--- LR Test Set Performance ---\n")
print_cm(m_lr)

# LR Odds Ratio plot
lr_coef <- broom::tidy(lr_model) |>
  filter(term != "(Intercept)") |>
  mutate(
    OR        = exp(estimate),
    CI_lower  = exp(estimate - 1.96*std.error),
    CI_upper  = exp(estimate + 1.96*std.error),
    direction = if_else(estimate>0,
                        "Reduces survival","Improves survival")
  ) |>
  arrange(desc(abs(estimate))) |>
  head(20)

cat("\n--- LR Top 20 coefficients ---\n")
print(lr_coef |> select(term,estimate,OR,CI_lower,CI_upper,p.value))

p_lr_imp <- ggplot(lr_coef |>
                     mutate(term=fct_reorder(term,estimate)),
                   aes(x=estimate, y=term, fill=direction)) +
  geom_col(width=0.7, colour="white") +
  geom_vline(xintercept=0, colour="grey40", linewidth=0.8) +
  scale_fill_manual(values=c("Improves survival"="#0F6E56",
                              "Reduces survival"="#EF4444")) +
  labs(title="Logistic Regression — Feature Coefficients",
       subtitle="Positive = increases risk of dying | Negative = improves survival",
       x="Log Odds Ratio", y=NULL, fill=NULL,
       caption="Top 20 features by absolute coefficient") +
  theme_minimal(base_size=11) +
  theme(plot.background=element_rect(fill="white",NA),
        legend.position="bottom")

print(p_lr_imp)
ggsave(file.path(PATHS$out_brca_figures,"fig08_lr_importance.png"),
       p_lr_imp, width=10, height=8, dpi=300)
cat("Saved: fig08_lr_importance.png\n")

plot_roc(m_lr$roc_obj, "Logistic Regression", m_lr$AUC,
         file.path(PATHS$out_brca_figures,"fig08_lr_roc.png"))
cat("Saved: fig08_lr_roc.png\n")

plot_cm(m_lr, "Logistic Regression",
        file.path(PATHS$out_brca_figures,"fig08_lr_confusion.png"))
cat("Saved: fig08_lr_confusion.png\n")

saveRDS(lr_model,
        file.path(PATHS$out_brca_models,"ml1_logistic_regression.rds"))
all_metrics[["LR"]] <- m_lr
roc_objects[["Logistic Regression"]] <- m_lr$roc_obj


# ═══════════════════════════════════════════════════════════════
# ML-2: ELASTIC NET (glmnet alpha=0.5)
# ═══════════════════════════════════════════════════════════════

cat("\n──────────────────────────────────────────────\n")
cat(" ML-2: Elastic Net (alpha=0.5, L1+L2)\n")
cat("──────────────────────────────────────────────\n")

set.seed(42)
en_cv <- cv.glmnet(X_train, y_train, alpha=0.5,
                   family="binomial", nfolds=5,
                   type.measure="auc")

cat(sprintf("Best lambda (min): %.6f\n", en_cv$lambda.min))
cat(sprintf("Best CV AUC      : %.4f\n", max(en_cv$cvm)))

en_prob <- predict(en_cv, X_test, s="lambda.min",
                   type="response")[,1]

m_en <- compute_metrics(y_test, en_prob, 0.5, "Elastic Net")
cat("\n--- Elastic Net Test Set Performance ---\n")
print_cm(m_en)

# EN coefficients
en_coef <- coef(en_cv, s="lambda.min")
en_coef_df <- tibble(
  feature   = rownames(en_coef)[-1],
  coef      = as.numeric(en_coef[-1])
) |> filter(coef != 0) |>
  mutate(direction=if_else(coef>0,"Reduces survival",
                           "Improves survival")) |>
  arrange(desc(abs(coef)))

cat(sprintf("Non-zero coefficients: %d\n", nrow(en_coef_df)))
cat("--- Elastic Net non-zero features ---\n")
print(en_coef_df)

p_en_imp <- ggplot(en_coef_df |>
                     mutate(feature=fct_reorder(feature,coef)),
                   aes(x=coef,y=feature,fill=direction)) +
  geom_col(width=0.7,colour="white") +
  geom_vline(xintercept=0,colour="grey40",linewidth=0.8) +
  scale_fill_manual(values=c("Improves survival"="#0F6E56",
                              "Reduces survival"="#EF4444")) +
  labs(title="Elastic Net — Non-Zero Coefficients",
       subtitle=sprintf("alpha=0.5 | lambda.min | %d features selected",
                        nrow(en_coef_df)),
       x="Coefficient", y=NULL, fill=NULL) +
  theme_minimal(base_size=11) +
  theme(plot.background=element_rect(fill="white",NA),
        legend.position="bottom")

print(p_en_imp)
ggsave(file.path(PATHS$out_brca_figures,"fig08_en_importance.png"),
       p_en_imp, width=10, height=7, dpi=300)
cat("Saved: fig08_en_importance.png\n")

plot_roc(m_en$roc_obj,"Elastic Net",m_en$AUC,
         file.path(PATHS$out_brca_figures,"fig08_en_roc.png"))
cat("Saved: fig08_en_roc.png\n")

plot_cm(m_en,"Elastic Net",
        file.path(PATHS$out_brca_figures,"fig08_en_confusion.png"))
cat("Saved: fig08_en_confusion.png\n")

saveRDS(en_cv,
        file.path(PATHS$out_brca_models,"ml2_elastic_net.rds"))
all_metrics[["EN"]] <- m_en
roc_objects[["Elastic Net"]] <- m_en$roc_obj


# ═══════════════════════════════════════════════════════════════
# ML-3: RANDOM FOREST
# ═══════════════════════════════════════════════════════════════

cat("\n──────────────────────────────────────────────\n")
cat(" ML-3: Random Forest\n")
cat("──────────────────────────────────────────────\n")

set.seed(42)
rf_model <- randomForest(
  x         = X_train,
  y         = factor(y_train, levels=c(0,1)),
  ntree     = 500,
  mtry      = floor(sqrt(ncol(X_train))),
  classwt   = c("0"=sum(y_train==1)/sum(y_train==0), "1"=1),
  importance= TRUE
)

cat("--- Random Forest summary ---\n")
print(rf_model)

cat(sprintf("\nOOB error rate : %.4f\n",
            rf_model$err.rate[500,"OOB"]))
cat(sprintf("Class 0 OOB    : %.4f\n",
            rf_model$err.rate[500,"0"]))
cat(sprintf("Class 1 OOB    : %.4f\n",
            rf_model$err.rate[500,"1"]))

rf_prob <- predict(rf_model, X_test, type="prob")[,2]
m_rf <- compute_metrics(y_test, rf_prob, 0.5, "Random Forest")
cat("\n--- RF Test Set Performance ---\n")
print_cm(m_rf)

# RF importance
rf_imp <- importance(rf_model) |>
  as.data.frame() |>
  rownames_to_column("feature") |>
  arrange(desc(MeanDecreaseAccuracy)) |>
  head(20)

cat("--- RF Top 20 features (MDA) ---\n")
print(rf_imp)

p_rf_imp <- ggplot(rf_imp |>
                     mutate(feature=fct_reorder(
                       feature,MeanDecreaseAccuracy)),
                   aes(x=MeanDecreaseAccuracy,y=feature)) +
  geom_col(fill="#185FA5",width=0.7,colour="white") +
  geom_text(aes(label=round(MeanDecreaseAccuracy,2)),
            hjust=-0.1,size=3.2) +
  scale_x_continuous(expand=expansion(mult=c(0,0.2))) +
  labs(title="Random Forest — Feature Importance (MDA)",
       subtitle="Mean Decrease in Accuracy — higher = more important",
       x="Mean Decrease Accuracy",y=NULL,
       caption="Top 20 features | ntree=500") +
  theme_minimal(base_size=11) +
  theme(plot.background=element_rect(fill="white",NA))

print(p_rf_imp)
ggsave(file.path(PATHS$out_brca_figures,"fig08_rf_importance.png"),
       p_rf_imp, width=10, height=8, dpi=300)
cat("Saved: fig08_rf_importance.png\n")

# OOB error convergence plot
oob_df <- tibble(trees=1:500,
                 OOB=rf_model$err.rate[,1],
                 Died=rf_model$err.rate[,2],
                 Survived=rf_model$err.rate[,3])

p_rf_oob <- oob_df |>
  pivot_longer(-trees,names_to="Class",values_to="Error") |>
  ggplot(aes(x=trees,y=Error,colour=Class,
             linewidth=Class)) +
  geom_line() +
  scale_colour_manual(values=c("OOB"="#1F2937",
                                "Died"="#EF4444",
                                "Survived"="#185FA5")) +
  scale_linewidth_manual(values=c("OOB"=1.2,
                                   "Died"=0.8,"Survived"=0.8)) +
  labs(title="Random Forest — OOB Error Convergence",
       subtitle="500 trees is sufficient when error stabilises",
       x="Number of Trees",y="Error Rate",colour=NULL,
       linewidth=NULL) +
  theme_minimal(base_size=12) +
  theme(plot.background=element_rect(fill="white",NA),
        legend.position="bottom")

print(p_rf_oob)
ggsave(file.path(PATHS$out_brca_figures,"fig08_rf_oob.png"),
       p_rf_oob, width=9, height=6, dpi=300)
cat("Saved: fig08_rf_oob.png\n")

plot_roc(m_rf$roc_obj,"Random Forest",m_rf$AUC,
         file.path(PATHS$out_brca_figures,"fig08_rf_roc.png"))
cat("Saved: fig08_rf_roc.png\n")

plot_cm(m_rf,"Random Forest",
        file.path(PATHS$out_brca_figures,"fig08_rf_confusion.png"))
cat("Saved: fig08_rf_confusion.png\n")

saveRDS(rf_model,
        file.path(PATHS$out_brca_models,"ml3_random_forest.rds"))
all_metrics[["RF"]] <- m_rf
roc_objects[["Random Forest"]] <- m_rf$roc_obj


# ═══════════════════════════════════════════════════════════════
# ML-4: XGBOOST
# ═══════════════════════════════════════════════════════════════

cat("\n──────────────────────────────────────────────\n")
cat(" ML-4: XGBoost\n")
cat("──────────────────────────────────────────────\n")

spw <- sum(y_train==1)/sum(y_train==0)
cat(sprintf("scale_pos_weight: %.2f (handles class imbalance)\n",spw))

dtrain <- xgb.DMatrix(X_train, label=y_train)
dtest  <- xgb.DMatrix(X_test,  label=y_test)

set.seed(42)
xgb_model <- xgb.train(
  params = list(
    objective        = "binary:logistic",
    eval_metric      = "auc",
    max_depth        = 4,
    eta              = 0.05,
    subsample        = 0.8,
    colsample_bytree = 0.8,
    min_child_weight = 5,
    scale_pos_weight = spw,
    seed             = 42
  ),
  data                  = dtrain,
  nrounds               = 300,
  evals                 = list(train=dtrain, test=dtest),
  early_stopping_rounds = 30,
  verbose               = 1,
  print_every_n         = 50
)

cat(sprintf("\nBest iteration : %d\n",xgb_model$best_iteration))
cat(sprintf("Best test AUC  : %.4f\n",xgb_model$best_score))

xgb_prob <- predict(xgb_model, dtest)
m_xgb <- compute_metrics(y_test, xgb_prob, 0.5, "XGBoost")
cat("\n--- XGBoost Test Set Performance ---\n")
print_cm(m_xgb)

# XGBoost importance
xgb_imp <- xgb.importance(feature_names=FEATURE_COLS,
                           model=xgb_model) |>
  as.data.frame() |>
  head(20)

cat("--- XGBoost Top 20 features (Gain) ---\n")
print(xgb_imp)

p_xgb_imp <- ggplot(xgb_imp |>
                      mutate(Feature=fct_reorder(Feature,Gain)),
                    aes(x=Gain,y=Feature)) +
  geom_col(fill="#0F6E56",width=0.7,colour="white") +
  geom_text(aes(label=round(Gain,3)),hjust=-0.1,size=3.2) +
  scale_x_continuous(expand=expansion(mult=c(0,0.2))) +
  labs(title="XGBoost — Feature Importance (Gain)",
       subtitle="Gain = relative contribution to split decisions",
       x="Gain",y=NULL,
       caption="Top 20 features") +
  theme_minimal(base_size=11) +
  theme(plot.background=element_rect(fill="white",NA))

print(p_xgb_imp)
ggsave(file.path(PATHS$out_brca_figures,"fig08_xgb_importance.png"),
       p_xgb_imp, width=10, height=8, dpi=300)
cat("Saved: fig08_xgb_importance.png\n")

plot_roc(m_xgb$roc_obj,"XGBoost",m_xgb$AUC,
         file.path(PATHS$out_brca_figures,"fig08_xgb_roc.png"))
cat("Saved: fig08_xgb_roc.png\n")

plot_cm(m_xgb,"XGBoost",
        file.path(PATHS$out_brca_figures,"fig08_xgb_confusion.png"))
cat("Saved: fig08_xgb_confusion.png\n")

xgb.save(xgb_model,
         file.path(PATHS$out_brca_models,"ml4_xgboost.model"))
all_metrics[["XGBoost"]] <- m_xgb
roc_objects[["XGBoost"]] <- m_xgb$roc_obj


# ═══════════════════════════════════════════════════════════════
# ML-5: SUPPORT VECTOR MACHINE
# ═══════════════════════════════════════════════════════════════

cat("\n──────────────────────────────────────────────\n")
cat(" ML-5: Support Vector Machine (RBF kernel)\n")
cat("──────────────────────────────────────────────\n")

set.seed(42)
svm_model <- svm(
  x          = X_train,
  y          = factor(y_train),
  kernel     = "radial",
  probability= TRUE,
  class.weights=c("0"=sum(y_train==1)/sum(y_train==0),
                  "1"=1),
  cost       = 1,
  gamma      = 1/ncol(X_train)
)

cat("--- SVM summary ---\n")
print(svm_model)

svm_pred <- predict(svm_model, X_test, probability=TRUE)
svm_prob <- attr(svm_pred, "probabilities")[,"1"]
m_svm <- compute_metrics(y_test, svm_prob, 0.5,
                          "Support Vector Machine")
cat("\n--- SVM Test Set Performance ---\n")
print_cm(m_svm)

plot_roc(m_svm$roc_obj,"Support Vector Machine",m_svm$AUC,
         file.path(PATHS$out_brca_figures,"fig08_svm_roc.png"))
cat("Saved: fig08_svm_roc.png\n")

plot_cm(m_svm,"Support Vector Machine",
        file.path(PATHS$out_brca_figures,"fig08_svm_confusion.png"))
cat("Saved: fig08_svm_confusion.png\n")

saveRDS(svm_model,
        file.path(PATHS$out_brca_models,"ml5_svm.rds"))
all_metrics[["SVM"]] <- m_svm
roc_objects[["SVM"]] <- m_svm$roc_obj


# ═══════════════════════════════════════════════════════════════
# COMPARISON: ROC OVERLAY + METRICS TABLE + BAR CHART
# ═══════════════════════════════════════════════════════════════

cat("\n──────────────────────────────────────────────\n")
cat(" SECTION A FINAL: Model Comparison\n")
cat("──────────────────────────────────────────────\n")

# Performance table
perf_tbl <- lapply(names(all_metrics), function(nm) {
  m <- all_metrics[[nm]]
  tibble(Model=nm, AUC=m$AUC, Accuracy=m$Accuracy,
         Recall=m$Recall, Precision=m$Precision,
         Specificity=m$Specificity, F1=m$F1, MCC=m$MCC,
         TP=m$TP, TN=m$TN, FP=m$FP, FN=m$FN)
}) |> bind_rows() |> arrange(desc(AUC))

cat("\n===== ML MODEL COMPARISON — TEST SET =====\n")
print(perf_tbl)

write_csv(perf_tbl,
          file.path(PATHS$out_brca_tables,
                    "08_model_performance_all.csv"))
cat("Saved: 08_model_performance_all.csv\n")

# ROC overlay
model_cols <- c(
  "Logistic Regression"="navy",
  "Elastic Net"        ="#5B21B6",
  "Random Forest"      ="#185FA5",
  "XGBoost"            ="#0F6E56",
  "SVM"                ="#F59E0B"
)

roc_all_df <- imap(roc_objects, function(roc_obj, nm) {
  auc_v <- as.numeric(pROC::auc(roc_obj))
  tibble(
    fpr  = 1 - roc_obj$specificities,
    tpr  = roc_obj$sensitivities,
    model= sprintf("%s (AUC=%.4f)", nm, auc_v)
  )
}) |> bind_rows()

p_roc_overlay <- ggplot(roc_all_df,
                         aes(x=fpr,y=tpr,colour=model)) +
  geom_line(linewidth=1.2) +
  geom_abline(linetype="dashed",colour="grey50") +
  scale_colour_manual(
    values=setNames(model_cols,
                    names(roc_objects) |>
                      map_chr(~sprintf("%s (AUC=%.4f)",
                                       .x,
                                       as.numeric(pROC::auc(
                                         roc_objects[[.x]])))))) +
  scale_x_continuous(labels=percent_format(),limits=c(0,1)) +
  scale_y_continuous(labels=percent_format(),limits=c(0,1)) +
  labs(title="ROC Curves — All ML Models",
       subtitle=sprintf(
         "Test set: n=%d | TCGA-BRCA | 5-year survival",
         length(y_test)),
       x="False Positive Rate (1-Specificity)",
       y="True Positive Rate (Sensitivity)",
       colour=NULL,
       caption="All models trained on same 80% split | Evaluated on same 20% test set") +
  theme_minimal(base_size=12) +
  theme(plot.background=element_rect(fill="white",NA),
        legend.position="right")

print(p_roc_overlay)
ggsave(file.path(PATHS$out_brca_figures,"fig08_roc_overlay.png"),
       p_roc_overlay, width=11, height=8, dpi=300)
cat("Printed and saved: fig08_roc_overlay.png\n")

# Metrics bar chart
p_compare <- perf_tbl |>
  mutate(Model=fct_reorder(Model,AUC)) |>
  select(Model,AUC,Accuracy,Recall,Specificity,F1,MCC) |>
  pivot_longer(-Model,names_to="Metric",values_to="Value") |>
  mutate(Metric=factor(Metric,
    levels=c("AUC","Accuracy","Recall",
             "Specificity","F1","MCC"))) |>
  ggplot(aes(x=Model,y=Value,fill=Metric)) +
  geom_col(position="dodge",width=0.75,colour="white") +
  geom_text(aes(label=round(Value,3)),
            position=position_dodge(0.75),
            vjust=-0.4,size=2.3,fontface="bold") +
  geom_hline(yintercept=0.5,linetype="dotted",
             colour="grey50",linewidth=0.6) +
  facet_wrap(~Metric,ncol=3,scales="free_y") +
  coord_flip() +
  scale_fill_brewer(palette="Blues") +
  scale_y_continuous(labels=number_format(accuracy=0.01),
                     expand=expansion(mult=c(0,0.25))) +
  labs(title="Model Performance Comparison — 5-Year Survival",
       subtitle=sprintf("Test set: n=%d patients",length(y_test)),
       x=NULL,y="Score",fill="Metric",
       caption="TCGA-BRCA | Confirmed-outcome cohort (n=305)") +
  theme_minimal(base_size=9) +
  theme(plot.background=element_rect(fill="white",NA),
        legend.position="none",
        strip.text=element_text(face="bold",size=9))

print(p_compare)
ggsave(file.path(PATHS$out_brca_figures,"fig08_model_comparison.png"),
       p_compare, width=14, height=10, dpi=300)
cat("Printed and saved: fig08_model_comparison.png\n")


# ══════════════════════════════════════════════════════════════
# ████████████████  SECTION B — COX REGRESSION  ██████████████
# ══════════════════════════════════════════════════════════════

cat("\n")
cat("╔══════════════════════════════════════════════╗\n")
cat("║        SECTION B — COX REGRESSION           ║\n")
cat("║  All 964 patients including censored         ║\n")
cat("╚══════════════════════════════════════════════╝\n\n")

cat("Cox handles censored patients mathematically.\n")
cat("No patient excluded — censoring accounted for\n")
cat("in the partial likelihood estimation.\n\n")

# Prepare Cox features (same features, no binary outcome)
cox_feature_cols <- setdiff(colnames(feat_cox),
                             c("patient_id","os_months",
                               "os_event"))

cox_df <- feat_cox |>
  select(patient_id, os_months, os_event,
         all_of(cox_feature_cols)) |>
  filter(!is.na(os_months), !is.na(os_event))

cat(sprintf("Cox dataset: %d patients, %d features\n",
            nrow(cox_df), length(cox_feature_cols)))

# Impute NA with median
for (col in cox_feature_cols) {
  na_n <- sum(is.na(cox_df[[col]]))
  if (na_n > 0) {
    cox_df[[col]][is.na(cox_df[[col]])] <-
      median(cox_df[[col]], na.rm=TRUE)
    cat(sprintf("Imputed %d NA in %s\n", na_n, col))
  }
}

# Univariate Cox — test each feature individually
cat("\n--- Univariate Cox p-values (top 20) ---\n")
cat("(Identifies which features associate with survival)\n\n")

uni_cox_results <- lapply(cox_feature_cols, function(f) {
  if (var(cox_df[[f]], na.rm=TRUE) == 0) return(NULL)
  form <- as.formula(sprintf(
    "Surv(os_months, os_event) ~ %s", f))
  tryCatch({
    fit <- coxph(form, data=cox_df)
    s   <- summary(fit)
    tibble(
      feature = f,
      HR      = round(exp(coef(fit)),3),
      CI_lower= round(exp(confint(fit)[1]),3),
      CI_upper= round(exp(confint(fit)[2]),3),
      p_value = round(s$coefficients[,"Pr(>|z|)"],4)
    )
  }, error=function(e) NULL)
}) |> bind_rows() |>
  arrange(p_value) |>
  head(20)

print(uni_cox_results)

# Multivariate Cox — significant features (p < 0.1)
sig_features <- uni_cox_results |>
  filter(p_value < 0.1) |>
  pull(feature)

cat(sprintf("\nSignificant features (p<0.1): %d\n",
            length(sig_features)))
print(sig_features)

if (length(sig_features) >= 2) {
  cox_formula <- as.formula(
    sprintf("Surv(os_months, os_event) ~ %s",
            paste(sig_features, collapse=" + ")))

  multi_cox <- coxph(cox_formula, data=cox_df)

  cat("\n--- Multivariate Cox summary ---\n")
  print(summary(multi_cox))

  # C-index
  cindex <- concordance(multi_cox)$concordance
  cat(sprintf("\nC-index (concordance): %.4f\n", cindex))
  cat("(0.5=random, 0.7+=good, 1.0=perfect)\n")

  # Hazard ratios table
  cox_hr <- tidy(multi_cox, exponentiate=TRUE,
                  conf.int=TRUE) |>
    mutate(
      HR       = round(estimate,3),
      CI_lower = round(conf.low,3),
      CI_upper = round(conf.high,3),
      p_value  = round(p.value,4),
      significance = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        p.value < 0.1   ~ ".",
        TRUE            ~ ""
      )
    ) |>
    select(term,HR,CI_lower,CI_upper,p_value,significance)

  cat("\n--- Multivariate Cox Hazard Ratios ---\n")
  print(cox_hr)

  write_csv(cox_hr,
            file.path(PATHS$out_brca_tables,
                      "08_cox_hazard_ratios.csv"))
  cat("Saved: 08_cox_hazard_ratios.csv\n")

  # Forest plot of hazard ratios
  p_cox_forest <- ggplot(cox_hr,
                          aes(x=HR,y=fct_reorder(term,HR))) +
    geom_vline(xintercept=1,linetype="dashed",
               colour="grey50",linewidth=0.9) +
    geom_point(size=3.5,colour="#185FA5") +
    geom_errorbarh(aes(xmin=CI_lower,xmax=CI_upper),
                   height=0.3,colour="#185FA5",
                   linewidth=0.9) +
    geom_text(aes(label=sprintf("HR=%.2f\np=%.3f",
                                HR,p_value)),
              hjust=-0.15,size=3,colour="#4B5563") +
    scale_x_continuous(
      limits=c(min(cox_hr$CI_lower)*0.8,
               max(cox_hr$CI_upper)*1.5)) +
    labs(title="Cox Regression — Hazard Ratios (Forest Plot)",
         subtitle=sprintf(
           "n=%d patients (incl censored) | C-index=%.3f",
           nrow(cox_df), cindex),
         x="Hazard Ratio (HR) with 95% CI",
         y=NULL,
         caption="HR>1=higher death risk | HR<1=protective | 964 patients") +
    theme_minimal(base_size=12) +
    theme(plot.background=element_rect(fill="white",NA))

  print(p_cox_forest)
  ggsave(file.path(PATHS$out_brca_figures,
                   "fig08_cox_forest_plot.png"),
         p_cox_forest, width=11, height=7, dpi=300)
  cat("Printed and saved: fig08_cox_forest_plot.png\n")

  # ── FIXED: Kaplan-Meier chunk ─────────────────────────────────
  
  cat("\n--- Kaplan-Meier curves by stage (fixed) ---\n")
  
  # Step 1: Load stage_clean from clinical_cleaned
  clinical_stage <- read_csv(
    file.path(PATHS$proc_brca, "clinical_cleaned.csv"),
    show_col_types=FALSE) |>
    select(patient_id, stage_clean)
  
  # Step 2: Join stage_clean onto cox_df cleanly
  km_data <- cox_df |>
    left_join(clinical_stage, by="patient_id") |>
    filter(
      stage_clean %in% c("Stage I","Stage II",
                         "Stage III","Stage IV"),
      !is.na(os_months),
      !is.na(os_event)
    )
  
  cat(sprintf("KM dataset: %d patients\n", nrow(km_data)))
  cat("--- Stage breakdown ---\n")
  print(table(km_data$stage_clean))
  
  # Step 3: Fit KM
  km_fit <- survfit(
    Surv(os_months, os_event) ~ stage_clean,
    data = km_data
  )
  
  print(km_fit)
  
  # Step 4: Plot
  p_km <- ggsurvplot(
    km_fit,
    data         = km_data,   # ← km_data not cox_df
    pval         = TRUE,
    conf.int     = FALSE,
    risk.table   = TRUE,
    title        = "Kaplan-Meier Survival Curves by Cancer Stage",
    subtitle     = "TCGA-BRCA | n=964 patients (all, incl censored)",
    xlab         = "Time (months)",
    ylab         = "Survival Probability",
    legend.title = "Stage",
    palette      = c("#5DCAA5","#185FA5","#F59E0B","#EF4444"),
    ggtheme      = theme_minimal(base_size=12)
  )
  
  print(p_km)
  
  png(file.path(PATHS$out_brca_figures, "fig08_km_curves.png"),
      width=3200, height=2400, res=300)
  print(p_km)
  dev.off()
  cat("Printed and saved: fig08_km_curves.png\n")
  
  png(file.path(PATHS$out_brca_figures,"fig08_km_curves.png"),
      width=3200, height=2400, res=300)
  print(p_km)
  dev.off()
  cat("Printed and saved: fig08_km_curves.png\n")

  saveRDS(multi_cox,
          file.path(PATHS$out_brca_models,"cox_model.rds"))
  cat("Saved: cox_model.rds\n")
}


# ═══════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════

cat("\n====================================================\n")
cat(" SCRIPT 08 COMPLETE\n")
cat("====================================================\n")
cat("\n  SECTION A — ML Models (305 patients):\n")
cat(sprintf("  Train: %d | Test: %d\n",
            nrow(X_train), nrow(X_test)))
cat(sprintf("  Features: %d\n", length(FEATURE_COLS)))
cat("\n  ── Performance (test set) ──\n")
cat(sprintf("  %-22s  AUC    Recall  Spec   F1\n","Model"))
cat(sprintf("  %s\n", strrep("-",58)))
for (nm in perf_tbl$Model) {
  r <- perf_tbl[perf_tbl$Model==nm,]
  cat(sprintf("  %-22s  %.4f  %.4f  %.4f  %.4f\n",
              nm, r$AUC, r$Recall, r$Specificity, r$F1))
}
cat("\n  SECTION B — Cox Regression (964 patients):\n")
if (exists("cindex"))
  cat(sprintf("  C-index: %.4f\n", cindex))
cat(sprintf("  Significant features: %d\n", length(sig_features)))
cat("\n  Plots printed + saved:\n")
cat("    Individual: ROC, CM, importance per model\n")
cat("    fig08_roc_overlay.png\n")
cat("    fig08_model_comparison.png\n")
cat("    fig08_cox_forest_plot.png\n")
cat("    fig08_km_curves.png\n")
cat("\nNext: source('scripts/09_test_second_dataset.R')\n")
cat("====================================================\n")

