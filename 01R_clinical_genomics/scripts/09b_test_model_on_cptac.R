# ============================================================
# Script 09b: Apply Trained Models to CPTAC-BRCA 2020
# Author  : Muteeba Azhar | MOI Internship | July 2026
#
# CPTAC available columns: PATIENT_ID, AGE only
# MISSING: STAGE, OS_STATUS, OS_MONTHS, PRIMARY_DIAGNOSIS
# Mode: PREDICTION ONLY (no AUC computation)
# ============================================================

library(tidyverse); library(pROC); library(randomForest)
library(xgboost);   library(glmnet); library(e1071)
library(scales);    library(pheatmap)

setwd("E:/R_Studio/Projects/MOI/01R_clinical_genomics")

PATHS <- list(
  raw_cptac        = "data/raw/TEST_brca_cptac_2020",
  proc_brca        = "data/processed/TCGA-BRCA",
  out_brca_tables  = "outputs/TCGA-BRCA/tables",
  out_brca_figures = "outputs/TCGA-BRCA/figures",
  out_brca_models  = "outputs/TCGA-BRCA/models"
)

NA_STRINGS <- c("","NA","N/A","n/a","[Not Available]",
                "[Not Applicable]","[Unknown]","[Discrepancy]")

TOP20 <- c("TP53","PIK3CA","TTN","CDH1","GATA3","MUC16",
           "KMT2C","MAP3K1","HMCN1","FLG","RYR2","SYNE1",
           "PTEN","USH2A","NCOR1","ZFHX4","SPTA1","DMD","MUC17","NEB")

CODING <- c("Missense_Mutation","Nonsense_Mutation",
            "Frame_Shift_Del","Frame_Shift_Ins","Splice_Site",
            "In_Frame_Del","In_Frame_Ins","Nonstop_Mutation",
            "Translation_Start_Site")

read_cbio <- function(p) read_tsv(p,comment="#",na=NA_STRINGS,
                                   show_col_types=FALSE,progress=FALSE)

cat("====================================================\n")
cat(" SCRIPT 09b: Trained Models → CPTAC-BRCA 2020\n")
cat(" Mode: PREDICTION ONLY\n")
cat(" Missing in CPTAC: STAGE, OS_STATUS, OS_MONTHS, DIAGNOSIS\n")
cat("====================================================\n\n")

# ── SECTION 1: Load models ────────────────────────────────────
cat("SECTION 1: Loading trained models\n")
lr_model  <- readRDS(file.path(PATHS$out_brca_models,"ml1_logistic_regression.rds"))
en_model  <- readRDS(file.path(PATHS$out_brca_models,"ml2_elastic_net.rds"))
rf_model  <- readRDS(file.path(PATHS$out_brca_models,"ml3_random_forest.rds"))
xgb_model <- xgb.load(file.path(PATHS$out_brca_models,"ml4_xgboost.model"))
svm_model <- readRDS(file.path(PATHS$out_brca_models,"ml5_svm.rds"))
scaler    <- readRDS(file.path(PATHS$out_brca_models,"feature_scaler.rds"))

FCOLS <- scaler$feature_cols
CCOLS <- scaler$cont_cols
cat(sprintf("Models loaded | %d features expected\n", length(FCOLS)))
cat("Feature columns:\n"); print(FCOLS)

# ── SECTION 2: Load CPTAC raw files ──────────────────────────
cat("\nSECTION 2: Loading CPTAC raw files\n")
clin <- read_cbio(file.path(PATHS$raw_cptac,"data_clinical_patient.txt"))
samp <- read_cbio(file.path(PATHS$raw_cptac,"data_clinical_sample.txt"))
maf  <- read_cbio(file.path(PATHS$raw_cptac,"data_mutations.txt"))

cat(sprintf("Clinical: %d x %d\n",nrow(clin),ncol(clin)))
cat(sprintf("Sample  : %d x %d\n",nrow(samp),ncol(samp)))
cat(sprintf("MAF     : %d x %d\n",nrow(maf),ncol(maf)))

cat("\n--- CPTAC clinical colnames ---\n"); print(colnames(clin))
cat("\n--- CPTAC clinical head(5) ---\n"); print(head(clin,5))
cat("\n--- CPTAC sample colnames ---\n");  print(colnames(samp))

# ── SECTION 3: Extract available clinical features ────────────
cat("\nSECTION 3: Extracting available features\n")
cat("Available: PATIENT_ID, AGE\n")
cat("Missing  : STAGE, OS_STATUS, OS_MONTHS, PRIMARY_DIAGNOSIS\n\n")

# Training medians for imputation (from scaler)
train_stage_median <- 2L   # Stage II most common in training
train_age_median   <- scaler$means["age"]  # training mean age

cat(sprintf("Impute stage with: %d (Stage II — training median)\n",
            train_stage_median))
cat(sprintf("Impute age with  : %.1f (training mean, if missing)\n",
            train_age_median))

cptac_clean <- tibble(
  patient_id    = str_trim(str_to_upper(clin$PATIENT_ID)),
  age           = as.numeric(clin$AGE),
  stage_numeric = train_stage_median,   # imputed — NOT in CPTAC
  stage_clean   = "Stage II (imputed)", # for labelling only
  primary_diagnosis = NA_character_     # not available
)

# Impute missing age with training mean
n_age_na <- sum(is.na(cptac_clean$age))
cptac_clean$age[is.na(cptac_clean$age)] <- train_age_median
cat(sprintf("\nAge: %d imputed with training mean (%.1f)\n",
            n_age_na, train_age_median))

cat("\n--- cptac_clean str() ---\n"); str(cptac_clean)
cat(sprintf("Patients: %d\n", nrow(cptac_clean)))
cat(sprintf("Age summary:\n")); print(summary(cptac_clean$age))

# Check sample file for TMB or stage
tmb_col   <- colnames(samp)[str_detect(str_to_upper(colnames(samp)),"TMB")]
stage_col <- colnames(samp)[str_detect(str_to_upper(colnames(samp)),"STAGE|TUMOR_STAGE")]

if(length(tmb_col)>0) {
  cat(sprintf("\nTMB found in sample file: %s\n",tmb_col[1]))
} else {
  cat("\nTMB not in sample file — will compute from MAF\n")
}

if(length(stage_col)>0) {
  cat(sprintf("Stage found in sample file: %s\n",stage_col[1]))
  cat("--- Stage values ---\n")
  print(table(samp[[stage_col[1]]],useNA="always"))
  # Update stage_numeric if found
  stage_map <- case_when(
    str_detect(coalesce(as.character(samp[[stage_col[1]]]),""),
               regex("stage iv",TRUE))  ~ 4L,
    str_detect(coalesce(as.character(samp[[stage_col[1]]]),""),
               regex("stage iii",TRUE)) ~ 3L,
    str_detect(coalesce(as.character(samp[[stage_col[1]]]),""),
               regex("stage ii",TRUE))  ~ 2L,
    str_detect(coalesce(as.character(samp[[stage_col[1]]]),""),
               regex("stage i",TRUE))   ~ 1L,
    TRUE ~ NA_integer_)

  stage_from_sample <- tibble(
    patient_id    = str_trim(str_to_upper(samp$PATIENT_ID)),
    stage_numeric_samp = stage_map
  ) |> group_by(patient_id) |> slice(1) |> ungroup()

  cptac_clean <- cptac_clean |>
    left_join(stage_from_sample, by="patient_id") |>
    mutate(stage_numeric = coalesce(stage_numeric_samp,
                                     as.integer(train_stage_median))) |>
    select(-stage_numeric_samp)

  cat(sprintf("\nStage updated from sample file\n"))
  print(table(cptac_clean$stage_numeric, useNA="always"))
}

# ── SECTION 4: Mutation counts + TMB ─────────────────────────
cat("\nSECTION 4: Mutation counts and TMB\n")
coding <- maf |>
  filter(Variant_Classification %in% CODING) |>
  mutate(patient_id=str_to_upper(str_trim(
    str_sub(Tumor_Sample_Barcode,1,12))))

cat(sprintf("Coding mutations: %d / %d\n",nrow(coding),nrow(maf)))
cat(sprintf("Unique patients in MAF: %d\n",n_distinct(coding$patient_id)))

mut_counts <- coding |> group_by(patient_id) |>
  summarise(
    n_total_coding   = n(),
    n_nonsense       = sum(Variant_Classification=="Nonsense_Mutation"),
    n_frameshift_del = sum(Variant_Classification=="Frame_Shift_Del"),
    n_frameshift_ins = sum(Variant_Classification=="Frame_Shift_Ins"),
    n_splice_site    = sum(Variant_Classification=="Splice_Site"),
    n_truncating     = sum(Variant_Classification %in%
      c("Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins")),
    tmb_nonsynonymous= round(n()/38,4),
    .groups="drop")

cat("--- TMB summary ---\n"); print(summary(mut_counts$tmb_nonsynonymous))
cat(sprintf("Mean TMB: %.3f | Median: %.3f\n",
    mean(mut_counts$tmb_nonsynonymous),
    median(mut_counts$tmb_nonsynonymous)))

# ── SECTION 5: Gene flags ─────────────────────────────────────
cat("\nSECTION 5: Gene flags (top 20 training genes)\n")
gf <- coding |> filter(Hugo_Symbol %in% TOP20) |>
  group_by(patient_id,Hugo_Symbol) |>
  summarise(m=1L,.groups="drop") |>
  pivot_wider(id_cols=patient_id,names_from=Hugo_Symbol,
              values_from=m,values_fill=0L)

for(g in TOP20) if(!g %in% colnames(gf)) gf[[g]] <- 0L

gf_full <- tibble(patient_id=cptac_clean$patient_id) |>
  left_join(gf,by="patient_id") |>
  mutate(across(all_of(TOP20),~coalesce(.,0L))) |>
  select(patient_id,all_of(TOP20))

cat(sprintf("Gene flags: %d x %d\n",nrow(gf_full),length(TOP20)))
cat("\n--- Gene frequencies in CPTAC ---\n")
print(gf_full |>
  summarise(across(all_of(TOP20),~round(mean(.)*100,1))) |>
  pivot_longer(everything(),names_to="Gene",values_to="pct") |>
  arrange(desc(pct)))

# ── SECTION 6: Diagnosis columns — all 0 (not in CPTAC) ───────
cat("\nSECTION 6: Diagnosis columns\n")

# Load training diagnosis column names
diag_ref <- tryCatch(
  read_csv(file.path(PATHS$out_brca_tables,"07_diagnosis_columns.csv"),
           show_col_types=FALSE)$column,
  error=function(e) {
    cat("07_diagnosis_columns.csv not found — detecting from FCOLS\n")
    FCOLS[str_starts(FCOLS,"diag_")]
  })

cat(sprintf("Training diag columns: %s\n", paste(diag_ref,collapse=", ")))
cat("All set to 0 for CPTAC (PRIMARY_DIAGNOSIS not available)\n")

diag_zero <- tibble(patient_id=cptac_clean$patient_id) |>
  bind_cols(lapply(diag_ref,function(col) tibble(!!col := 0L)) |>
              bind_cols())

cat("\n--- Diagnosis encoding head(3) ---\n"); print(head(diag_zero,3))

# ── SECTION 7: Assemble feature table ─────────────────────────
cat("\nSECTION 7: Assembling CPTAC feature table\n")

gcols <- c("tmb_nonsynonymous","n_truncating","n_nonsense",
           "n_frameshift_del","n_frameshift_ins","n_splice_site")

cptac_feat <- cptac_clean |>
  select(patient_id,age,stage_numeric) |>
  left_join(mut_counts|>select(patient_id,all_of(gcols)),
            by="patient_id") |>
  left_join(gf_full,by="patient_id") |>
  left_join(diag_zero,by="patient_id") |>
  mutate(across(all_of(gcols),~coalesce(.,0)))

cat(sprintf("cptac_feat: %d x %d\n",nrow(cptac_feat),ncol(cptac_feat)))
cat("\n--- str() ---\n"); str(cptac_feat)
cat("\n--- head(3) ---\n"); print(head(cptac_feat,3))
cat("\n--- Missing values ---\n")
for(f in c("age","stage_numeric","tmb_nonsynonymous"))
  cat(sprintf("  %-25s: %d missing\n",f,sum(is.na(cptac_feat[[f]]))))

write_csv(cptac_feat,
  file.path(PATHS$out_brca_tables,"09b_cptac_feature_table.csv"))
cat("Saved: 09b_cptac_feature_table.csv\n")

# ── SECTION 8: Build model matrix ────────────────────────────
cat("\nSECTION 8: Building model input matrix\n")

# Add any missing FCOLS as 0
for(col in FCOLS) {
  if(!col %in% colnames(cptac_feat)) {
    cptac_feat[[col]] <- 0
    cat(sprintf("  Added missing: %s = 0\n",col))
  }
}

X_raw <- cptac_feat[,FCOLS]

# Impute with TRAINING means (not CPTAC means)
cat("Imputing NA with training means:\n")
for(col in CCOLS) {
  na_n <- sum(is.na(X_raw[[col]]))
  if(na_n>0) {
    X_raw[[col]][is.na(X_raw[[col]])] <- scaler$means[col]
    cat(sprintf("  %s: %d imputed → %.3f\n",col,na_n,scaler$means[col]))
  }
}

# Apply training scaler
cat("Applying training scaler:\n")
X_sc <- as.data.frame(X_raw)
for(col in CCOLS)
  X_sc[[col]] <- (X_raw[[col]]-scaler$means[col])/scaler$sds[col]

X_cptac <- as.matrix(X_sc)
cat(sprintf("X_cptac: %d x %d\n",nrow(X_cptac),ncol(X_cptac)))

# ── SECTION 9: Apply all 5 models ────────────────────────────
cat("\nSECTION 9: Applying all 5 trained models\n")

# LR — add NA. if model expects it
lr_terms <- names(coef(lr_model))
X_lr     <- as.data.frame(X_cptac)
if("NA." %in% lr_terms && !"NA." %in% colnames(X_lr)) {
  X_lr[["NA."]] <- 0
  cat("  LR: added NA. column = 0\n")
}

# EN — add NA. if model expects it
en_nvar <- length(en_model$glmnet.fit$beta[,1])
X_en    <- X_cptac
if(en_nvar == ncol(X_en)+1) {
  X_en_df        <- as.data.frame(X_en)
  X_en_df[["NA."]] <- 0
  X_en           <- as.matrix(X_en_df)
  cat("  EN: added NA. column = 0\n")
}

# RF — reorder to match training column names exactly
rf_exp <- rownames(rf_model$importance)
X_rf   <- as.data.frame(X_cptac)
for(col in setdiff(rf_exp,colnames(X_rf)))
  X_rf[[col]] <- 0
X_rf <- as.matrix(X_rf[,rf_exp,drop=FALSE])
cat(sprintf("  RF: %d features matched\n",ncol(X_rf)))

# SVM — reorder to match training column names
svm_exp <- colnames(svm_model$SV)
X_svm   <- as.data.frame(X_cptac)
for(col in setdiff(svm_exp,colnames(X_svm)))
  X_svm[[col]] <- 0
X_svm <- as.matrix(X_svm[,svm_exp,drop=FALSE])
cat(sprintf("  SVM: %d features matched\n",ncol(X_svm)))

# XGBoost
X_xgb <- X_cptac

# Predict
p_lr  <- predict(lr_model,newdata=X_lr,type="response")
cat(sprintf("LR  done: %.3f – %.3f\n",min(p_lr),max(p_lr)))

p_en  <- predict(en_model,X_en,s="lambda.min",type="response")[,1]
cat(sprintf("EN  done: %.3f – %.3f\n",min(p_en),max(p_en)))

p_rf  <- predict(rf_model,X_rf,type="prob")[,2]
cat(sprintf("RF  done: %.3f – %.3f\n",min(p_rf),max(p_rf)))

p_xgb <- predict(xgb_model,xgb.DMatrix(X_xgb))
cat(sprintf("XGB done: %.3f – %.3f\n",min(p_xgb),max(p_xgb)))

p_svm <- attr(predict(svm_model,X_svm,probability=TRUE),
              "probabilities")[,"1"]
cat(sprintf("SVM done: %.3f – %.3f\n",min(p_svm),max(p_svm)))

p_ens <- rowMeans(cbind(p_lr,p_en,p_rf,p_xgb,p_svm))
cat(sprintf("ENS done: %.3f – %.3f\n",min(p_ens),max(p_ens)))

# ── SECTION 10: Compile predictions ──────────────────────────
cat("\nSECTION 10: Compiling predictions\n")

preds <- cptac_clean |>
  select(patient_id,age,stage_clean,stage_numeric) |>
  left_join(mut_counts|>select(patient_id,tmb_nonsynonymous),
            by="patient_id") |>
  mutate(tmb_nonsynonymous=coalesce(tmb_nonsynonymous,0)) |>
  mutate(
    prob_LR       = round(p_lr,4),
    prob_EN       = round(p_en,4),
    prob_RF       = round(p_rf,4),
    prob_XGB      = round(p_xgb,4),
    prob_SVM      = round(p_svm,4),
    prob_Ensemble = round(p_ens,4),
    pred_LR  = as.integer(p_lr >=.5),
    pred_EN  = as.integer(p_en >=.5),
    pred_RF  = as.integer(p_rf >=.5),
    pred_XGB = as.integer(p_xgb>=.5),
    pred_SVM = as.integer(p_svm>=.5),
    pred_Ens = as.integer(p_ens>=.5),
    n_models_agree = pred_LR+pred_EN+pred_RF+pred_XGB+pred_SVM,
    risk_category  = case_when(
      p_ens>=.75 ~ "High survival predicted",
      p_ens>=.50 ~ "Moderate survival predicted",
      TRUE       ~ "Low survival predicted")
  )

cat(sprintf("preds: %d patients\n",nrow(preds)))
cat("\n--- head(10) ---\n")
print(head(preds|>select(patient_id,age,stage_clean,
                          prob_LR,prob_RF,prob_XGB,
                          prob_Ensemble,risk_category),10))
cat("\n--- Risk category distribution ---\n")
print(table(preds$risk_category))
cat("\n--- Model agreement (0=all predict death, 5=all predict survival) ---\n")
print(table(preds$n_models_agree))

write_csv(preds,
  file.path(PATHS$out_brca_tables,"09b_cptac_predictions.csv"))
cat("Saved: 09b_cptac_predictions.csv\n")

# ── SECTION 11: PLOTS ─────────────────────────────────────────
cat("\nSECTION 11: Plots (printed then saved)\n")

# PLOT 1: Risk distribution — 6 panels
p1 <- preds |>
  select(patient_id,prob_LR,prob_EN,prob_RF,
         prob_XGB,prob_SVM,prob_Ensemble) |>
  pivot_longer(-patient_id,names_to="Model",values_to="Prob") |>
  mutate(Model=recode(Model,"prob_LR"="Logistic Regression",
    "prob_EN"="Elastic Net","prob_RF"="Random Forest",
    "prob_XGB"="XGBoost","prob_SVM"="SVM",
    "prob_Ensemble"="Ensemble")) |>
  ggplot(aes(x=Prob,fill=Model))+
  geom_histogram(bins=25,colour="white",linewidth=.1,alpha=.85)+
  geom_vline(xintercept=.5,linetype="dashed",
             colour="#1F2937",linewidth=.9)+
  facet_wrap(~Model,ncol=3)+
  scale_fill_manual(values=c("Logistic Regression"="navy",
    "Elastic Net"="#5B21B6","Random Forest"="#185FA5",
    "XGBoost"="#0F6E56","SVM"="#F59E0B","Ensemble"="#EF4444"))+
  scale_x_continuous(labels=percent_format())+
  labs(title="Predicted 5-Year Survival Probability — CPTAC-BRCA 2020",
       subtitle=sprintf("n=%d | 5 models + ensemble | Dashed = 0.5 threshold",nrow(preds)),
       x="P(Survived 5 years)",y="Patients",
       caption=paste0("Trained on TCGA-BRCA (n=305 confirmed outcomes)\n",
         "Applied to CPTAC-BRCA 2020 — prediction only (no OS data available)"))+
  theme_minimal(base_size=11)+
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="none",
        strip.text=element_text(face="bold",size=9))

print(p1)
ggsave(file.path(PATHS$out_brca_figures,"fig09b_risk_distribution.png"),
       p1,width=13,height=9,dpi=300)
cat("Saved: fig09b_risk_distribution.png\n")

# PLOT 2: Model agreement scatter
p2 <- ggplot(preds,aes(x=prob_XGB,y=prob_Ensemble,
                        colour=risk_category))+
  geom_point(size=3,alpha=.8)+
  geom_abline(linetype="dashed",colour="grey50")+
  geom_hline(yintercept=.5,linetype="dotted",colour="#9CA3AF")+
  geom_vline(xintercept=.5,linetype="dotted",colour="#9CA3AF")+
  scale_colour_manual(values=c(
    "High survival predicted"="#0F6E56",
    "Moderate survival predicted"="#F59E0B",
    "Low survival predicted"="#EF4444"))+
  scale_x_continuous(labels=percent_format(),limits=c(0,1))+
  scale_y_continuous(labels=percent_format(),limits=c(0,1))+
  labs(title="Model Agreement — XGBoost vs Ensemble",
       subtitle="CPTAC-BRCA 2020 | Diagonal = perfect agreement",
       x="XGBoost Probability",y="Ensemble Probability",
       colour="Predicted Survival")+
  theme_minimal(base_size=12)+
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="bottom")

print(p2)
ggsave(file.path(PATHS$out_brca_figures,"fig09b_model_agreement.png"),
       p2,width=8,height=7,dpi=300)
cat("Saved: fig09b_model_agreement.png\n")

# PLOT 3: Age vs Predicted risk
p3 <- ggplot(preds,aes(x=age,y=prob_Ensemble,colour=risk_category))+
  geom_point(size=2.5,alpha=.75)+
  geom_smooth(method="loess",se=TRUE,colour="#1F2937",
              fill="grey80",linewidth=1)+
  geom_hline(yintercept=.5,linetype="dashed",
             colour="#EF4444",linewidth=.9)+
  scale_colour_manual(values=c(
    "High survival predicted"="#0F6E56",
    "Moderate survival predicted"="#F59E0B",
    "Low survival predicted"="#EF4444"))+
  scale_y_continuous(labels=percent_format(),limits=c(0,1))+
  labs(title="Age vs Predicted 5-Year Survival — CPTAC-BRCA 2020",
       subtitle=sprintf("n=%d | Ensemble model | Trend line = LOESS",
                        nrow(preds)),
       x="Age at Diagnosis",y="P(Survived 5 years)",
       colour="Predicted Survival",
       caption="Red dashed = 0.5 decision threshold")+
  theme_minimal(base_size=12)+
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="bottom")

print(p3)
ggsave(file.path(PATHS$out_brca_figures,"fig09b_age_vs_risk.png"),
       p3,width=9,height=6,dpi=300)
cat("Saved: fig09b_age_vs_risk.png\n")

# PLOT 4: TMB vs Predicted risk
p4 <- ggplot(preds,aes(x=tmb_nonsynonymous,y=prob_Ensemble,
                        colour=risk_category))+
  geom_point(size=2.5,alpha=.75)+
  geom_smooth(method="loess",se=TRUE,colour="#1F2937",
              fill="grey80",linewidth=1)+
  geom_hline(yintercept=.5,linetype="dashed",
             colour="#EF4444",linewidth=.9)+
  scale_colour_manual(values=c(
    "High survival predicted"="#0F6E56",
    "Moderate survival predicted"="#F59E0B",
    "Low survival predicted"="#EF4444"))+
  scale_x_continuous(trans="log1p",
    breaks=c(0,1,2,5,10,20,50))+
  scale_y_continuous(labels=percent_format(),limits=c(0,1))+
  labs(title="TMB vs Predicted Risk — CPTAC-BRCA 2020",
       subtitle="log1p scale | Higher TMB patients may have different risk",
       x="TMB (mut/Mb, log1p)",y="P(Survived 5 years)",
       colour="Predicted Survival")+
  theme_minimal(base_size=12)+
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="bottom")

print(p4)
ggsave(file.path(PATHS$out_brca_figures,"fig09b_tmb_vs_risk.png"),
       p4,width=9,height=6,dpi=300)
cat("Saved: fig09b_tmb_vs_risk.png\n")

# PLOT 5: TCGA test set performance (reference)
tcga_perf <- read_csv(
  file.path(PATHS$out_brca_tables,"08_model_performance_all.csv"),
  show_col_types=FALSE)

p5 <- tcga_perf |>
  select(Model,AUC,Recall,Specificity,F1) |>
  pivot_longer(-Model,names_to="Metric",values_to="Value") |>
  mutate(Metric=factor(Metric,levels=c("AUC","Recall","Specificity","F1")),
         Model=fct_reorder(Model,Value)) |>
  ggplot(aes(x=Model,y=Value,fill=Model))+
  geom_col(width=.7,colour="white")+
  geom_text(aes(label=round(Value,3)),vjust=-0.4,size=3,fontface="bold")+
  facet_wrap(~Metric,ncol=2)+
  coord_flip()+
  scale_y_continuous(expand=expansion(mult=c(0,.2)))+
  scale_fill_manual(values=c(
    "LR"="navy","EN"="#5B21B6","RF"="#185FA5",
    "XGBoost"="#0F6E56","SVM"="#F59E0B"))+
  labs(title="TCGA-BRCA Test Set Performance (Reference)",
       subtitle=paste0("Internal validation only — CPTAC used for prediction\n",
         "CPTAC AUC not computed: OS_STATUS and OS_MONTHS not available"),
       x=NULL,y="Score",
       caption="Note: CPTAC-BRCA 2020 lacks survival data — external AUC validation not possible")+
  theme_minimal(base_size=9)+
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="none",
        strip.text=element_text(face="bold",size=10))

print(p5)
ggsave(file.path(PATHS$out_brca_figures,"fig09b_tcga_reference_performance.png"),
       p5,width=11,height=8,dpi=300)
cat("Saved: fig09b_tcga_reference_performance.png\n")

# PLOT 6: Gene flags heatmap
fm  <- gf_full|>column_to_rownames("patient_id")|>as.matrix()
ann <- preds|>select(patient_id,risk_category)|>
  column_to_rownames("patient_id")
ann_c <- list(risk_category=c(
  "High survival predicted"="#0F6E56",
  "Moderate survival predicted"="#F59E0B",
  "Low survival predicted"="#EF4444"))

pheatmap(t(fm),color=c("white","#0F6E56"),breaks=c(-.5,.5,1.5),
         annotation_col=ann,annotation_colors=ann_c,
         show_colnames=FALSE,cluster_rows=TRUE,cluster_cols=TRUE,
         main="Gene Mutation Flags — CPTAC-BRCA 2020",
         fontsize=8,border_color="grey90")

png(file.path(PATHS$out_brca_figures,"fig09b_gene_flags_heatmap.png"),
    width=3000,height=2400,res=300)
pheatmap(t(fm),color=c("white","#0F6E56"),breaks=c(-.5,.5,1.5),
         annotation_col=ann,annotation_colors=ann_c,
         show_colnames=FALSE,cluster_rows=TRUE,cluster_cols=TRUE,
         main="Gene Mutation Flags — CPTAC-BRCA 2020",
         fontsize=8,border_color="grey90")
dev.off()
cat("Saved: fig09b_gene_flags_heatmap.png\n")

# ── FINAL SUMMARY ─────────────────────────────────────────────
cat("\n====================================================\n")
cat(" SCRIPT 09b COMPLETE\n")
cat("====================================================\n")
cat(sprintf("  CPTAC patients    : %d\n",nrow(preds)))
cat(sprintf("  Models applied    : 5 + ensemble\n"))
cat(sprintf("  Mode              : PREDICTION ONLY\n"))
cat(sprintf("  OS data available : NO\n"))
cat(sprintf("  AUC computed      : NO (no outcome to compare)\n"))
cat("\n  Risk distribution:\n")
rt <- table(preds$risk_category)
for(n in names(rt))
  cat(sprintf("    %-35s: %d (%.1f%%)\n",
              n,rt[n],rt[n]/sum(rt)*100))
cat("\n  Model agreement (5=all agree survival):\n")
print(table(preds$n_models_agree))
cat("\n  Figures saved:\n")
for(f in c("fig09b_risk_distribution.png",
           "fig09b_model_agreement.png",
           "fig09b_age_vs_risk.png",
           "fig09b_tmb_vs_risk.png",
           "fig09b_tcga_reference_performance.png",
           "fig09b_gene_flags_heatmap.png"))
  cat(sprintf("    %s\n",f))
cat("\n  Limitation statement for report:\n")
cat("  'The CPTAC-BRCA 2020 dataset does not contain\n")
cat("   OS_STATUS or OS_MONTHS. Therefore, formal AUC\n")
cat("   evaluation was not possible. The pipeline was\n")
cat("   applied to predict 5-year survival probability\n")
cat("   for each patient. Stage was imputed (Stage II)\n")
cat("   and diagnosis features set to 0, as these were\n")
cat("   not available in the CPTAC clinical file.'\n")
cat("====================================================\n")
