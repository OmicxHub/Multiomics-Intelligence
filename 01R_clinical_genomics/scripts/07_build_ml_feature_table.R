# ============================================================
# Script 07: ML-Ready Feature Table Construction
# Project : Clinical-Genomics Harmonisation Pipeline
# Author  : Muteeba Azhar | MOI Internship | July 2026
# FIX     : Remove "NA" column from diagnosis one-hot
# ============================================================

library(tidyverse)
library(scales)

setwd("E:/R_Studio/Projects/MOI/01R_clinical_genomics")

PATHS <- list(
  proc_brca        = "data/processed/TCGA-BRCA",
  out_brca_tables  = "outputs/TCGA-BRCA/tables",
  out_brca_figures = "outputs/TCGA-BRCA/figures"
)

TOP20_GENES <- c("TP53","PIK3CA","TTN","CDH1","GATA3","MUC16",
                 "KMT2C","MAP3K1","HMCN1","FLG","RYR2","SYNE1",
                 "PTEN","USH2A","NCOR1","ZFHX4","SPTA1","DMD",
                 "MUC17","NEB")

cat("====================================================\n")
cat(" SCRIPT 07: ML Feature Table Construction\n")
cat("====================================================\n\n")

# ── Load data ─────────────────────────────────────────────────
cat("SECTION 1: Loading data\n")
merged <- read_csv(file.path(PATHS$proc_brca,
  "clinical_genomics_merged.csv"), show_col_types=FALSE)
genomics_cleaned <- read_csv(file.path(PATHS$proc_brca,
  "genomics_cleaned.csv"), show_col_types=FALSE)

cat(sprintf("merged          : %d x %d\n",nrow(merged),ncol(merged)))
cat(sprintf("genomics_cleaned: %d mutations\n",nrow(genomics_cleaned)))
cat("\n--- merged str() ---\n"); str(merged)
cat("\n--- merged head(3) ---\n"); print(head(merged,3))

# ── five_year_survival — CORRECT definition ───────────────────
cat("\nSECTION 2: Defining five_year_survival\n")
cat("60 months = 5 years × 12. OS_MONTHS from diagnosis date.\n\n")

merged <- merged |>
  mutate(
    five_year_survival = case_when(
      os_event==0 & os_months >= 60 ~ 1L,
      os_event==1 & os_months >= 60 ~ 1L,
      os_event==1 & os_months <  60 ~ 0L,
      os_event==0 & os_months <  60 ~ NA_integer_
    ),
    outcome_label = case_when(
      os_event==0 & os_months >= 60 ~ "Confirmed survived (alive >=60m)",
      os_event==1 & os_months >= 60 ~ "Confirmed survived (died after 5yr)",
      os_event==1 & os_months <  60 ~ "Confirmed died within 5yr",
      os_event==0 & os_months <  60 ~ "CENSORED — excluded from ML"
    )
  )

outcome_summary <- merged |>
  count(outcome_label, five_year_survival) |>
  mutate(pct=round(n/nrow(merged)*100,1),
         in_ML=!is.na(five_year_survival)) |>
  arrange(desc(in_ML), five_year_survival)

cat("--- five_year_survival categories ---\n")
print(outcome_summary)

n_ml       <- sum(!is.na(merged$five_year_survival))
n_survived <- sum(merged$five_year_survival==1, na.rm=TRUE)
n_died     <- sum(merged$five_year_survival==0, na.rm=TRUE)
n_censored <- sum(is.na(merged$five_year_survival))

cat(sprintf("\n  Total matched  : %d\n", nrow(merged)))
cat(sprintf("  ML cohort      : %d (%.1f%%)\n", n_ml, n_ml/nrow(merged)*100))
cat(sprintf("  Survived (=1)  : %d (%.1f%%)\n", n_survived, n_survived/n_ml*100))
cat(sprintf("  Died (=0)      : %d (%.1f%%)\n", n_died, n_died/n_ml*100))
cat(sprintf("  Class balance  : %.1f:1\n", n_survived/n_died))
cat(sprintf("  Censored       : %d (excluded)\n", n_censored))

write_csv(outcome_summary,
  file.path(PATHS$out_brca_tables,"07_five_year_outcome_summary.csv"))

# ── Gene flags ────────────────────────────────────────────────
cat("\nSECTION 3: Gene flags (top 20 genes)\n")

gene_flags <- genomics_cleaned |>
  filter(Hugo_Symbol %in% TOP20_GENES) |>
  group_by(patient_id, Hugo_Symbol) |>
  summarise(mutated=1L, .groups="drop") |>
  pivot_wider(id_cols=patient_id, names_from=Hugo_Symbol,
              values_from=mutated, values_fill=0L)

# Add missing genes
for(g in TOP20_GENES)
  if(!g %in% colnames(gene_flags)) gene_flags[[g]] <- 0L

# Ensure all 964 patients
gene_flags_complete <- tibble(patient_id=merged$patient_id) |>
  left_join(gene_flags, by="patient_id") |>
  mutate(across(all_of(TOP20_GENES), ~coalesce(.,0L))) |>
  select(patient_id, all_of(TOP20_GENES))

cat(sprintf("gene_flags: %d x %d\n",
    nrow(gene_flags_complete), ncol(gene_flags_complete)-1))

gene_freq <- gene_flags_complete |>
  summarise(across(all_of(TOP20_GENES),~round(mean(.)*100,1))) |>
  pivot_longer(everything(),names_to="gene",values_to="pct_mutated") |>
  arrange(desc(pct_mutated))
cat("\n--- Gene frequencies ---\n"); print(gene_freq)

write_csv(gene_freq,
  file.path(PATHS$out_brca_tables,"07_gene_flag_frequency.csv"))

# ── Diagnosis one-hot ─────────────────────────────────────────
cat("\nSECTION 4: Encoding primary_diagnosis\n")
cat("--- primary_diagnosis values ---\n")
print(table(merged$primary_diagnosis, useNA="always"))

diag_encoded <- merged |>
  select(patient_id, primary_diagnosis) |>
  mutate(
    diag_col = primary_diagnosis |>
      str_replace_all("[^a-zA-Z0-9]","_") |>
      str_replace_all("_+","_") |>
      str_to_lower() |>
      str_c("diag_", x=_)
  ) |>
  mutate(value=1L) |>
  pivot_wider(id_cols=patient_id,
              names_from=diag_col,
              values_from=value,
              values_fill=0L)

# ══ CRITICAL FIX: Remove "NA" column from invalid diagnoses ════
before_cols <- ncol(diag_encoded) - 1
diag_encoded <- diag_encoded |>
  select(-any_of(c("NA","na","diag_na","diag_NA","")))
after_cols <- ncol(diag_encoded) - 1
cat(sprintf("\nDiagnosis columns: %d (removed %d invalid)\n",
            after_cols, before_cols - after_cols))
# ══════════════════════════════════════════════════════════════

diag_cols <- setdiff(colnames(diag_encoded), "patient_id")
cat(sprintf("Diagnosis columns kept: %s\n",
            paste(diag_cols, collapse=", ")))

# Save diagnosis column reference (needed by Script 09b)
write_csv(tibble(column=diag_cols),
  file.path(PATHS$out_brca_tables,"07_diagnosis_columns.csv"))
cat("Saved: 07_diagnosis_columns.csv\n")

# ── Assemble feature table ────────────────────────────────────
cat("\nSECTION 5: Assembling feature table\n")

base_features <- merged |>
  select(patient_id, age, stage_numeric,
         tmb_nonsynonymous, n_truncating, n_nonsense,
         n_frameshift_del, n_frameshift_ins, n_splice_site,
         os_months, os_event, five_year_survival)

feature_table <- base_features |>
  left_join(gene_flags_complete, by="patient_id") |>
  left_join(diag_encoded, by="patient_id")

cat(sprintf("feature_table: %d x %d\n",
    nrow(feature_table), ncol(feature_table)))
cat("\n--- column names ---\n"); print(colnames(feature_table))
cat("\n--- str() ---\n"); str(feature_table)
cat("\n--- head(3) ---\n"); print(head(feature_table,3))

# Missing values check
cat("\n--- Missing values in continuous features ---\n")
for(f in c("age","stage_numeric","tmb_nonsynonymous",
           "n_truncating","n_nonsense")) {
  cat(sprintf("  %-25s: %d missing\n",f,sum(is.na(feature_table[[f]]))))
}

# ── ML dataset (305 patients — confirmed outcomes only) ───────
cat("\nSECTION 6: ML dataset (confirmed outcomes only)\n")

feature_table_ml <- feature_table |>
  filter(!is.na(five_year_survival)) |>
  select(-os_months, -os_event)

cat(sprintf("feature_table_ml: %d x %d\n",
    nrow(feature_table_ml), ncol(feature_table_ml)))
cat(sprintf("  Survived (=1): %d\n",
    sum(feature_table_ml$five_year_survival==1)))
cat(sprintf("  Died (=0)    : %d\n",
    sum(feature_table_ml$five_year_survival==0)))

cat("\n--- ML column names ---\n")
print(colnames(feature_table_ml))

write_csv(feature_table_ml,
  file.path(PATHS$proc_brca,"feature_table_ml.csv"))
cat(sprintf("Saved: feature_table_ml.csv\n"))

# ── Cox dataset (all 964) ─────────────────────────────────────
cat("\nSECTION 7: Cox dataset (all 964 patients)\n")

feature_table_cox <- feature_table |>
  select(-five_year_survival)

cat(sprintf("feature_table_cox: %d x %d\n",
    nrow(feature_table_cox), ncol(feature_table_cox)))

write_csv(feature_table_cox,
  file.path(PATHS$proc_brca,"feature_table_cox.csv"))
cat(sprintf("Saved: feature_table_cox.csv\n"))

# ── Feature summary ───────────────────────────────────────────
cat("\nSECTION 8: Feature summary\n")

ml_feats <- feature_table_ml |>
  select(-patient_id,-five_year_survival)

feat_summary <- tibble(
  Feature  = colnames(ml_feats),
  Type = case_when(
    Feature %in% c("age","tmb_nonsynonymous","n_truncating",
                   "n_nonsense","n_frameshift_del",
                   "n_frameshift_ins","n_splice_site") ~ "Continuous",
    Feature=="stage_numeric" ~ "Ordinal (1-4)",
    Feature %in% TOP20_GENES ~ "Binary (gene flag)",
    str_starts(Feature,"diag_") ~ "Binary (diagnosis)",
    TRUE ~ "Other"),
  N_missing = sapply(ml_feats, function(x) sum(is.na(x))),
  Mean      = sapply(ml_feats, function(x) round(mean(x,na.rm=TRUE),3))
)
cat("--- Feature summary ---\n"); print(feat_summary, n=Inf)

write_csv(feat_summary,
  file.path(PATHS$out_brca_tables,"07_feature_summary.csv"))

# ── Plots ─────────────────────────────────────────────────────
cat("\nSECTION 9: Plots\n")

# Plot 1: Class balance
p1 <- feature_table_ml |>
  mutate(Outcome=if_else(five_year_survival==1,
    "Survived (=1)","Did not survive (=0)")) |>
  count(Outcome) |>
  mutate(pct=round(n/sum(n)*100,1)) |>
  ggplot(aes(x=Outcome,y=n,fill=Outcome)) +
  geom_col(width=.55,colour="white") +
  geom_text(aes(label=sprintf("%d\n(%.1f%%)",n,pct)),
            vjust=-0.3,size=5,fontface="bold") +
  scale_fill_manual(
    values=c("Survived (=1)"="#0F6E56",
             "Did not survive (=0)"="#EF4444")) +
  scale_y_continuous(expand=expansion(mult=c(0,.2))) +
  labs(title="Target Variable — Five-Year Survival",
       subtitle=sprintf("ML cohort: n=%d confirmed outcomes | Balance: %.1f:1",
                        n_ml, n_survived/n_died),
       x=NULL, y="Patients",
       caption=sprintf(
         "0=deceased <60m | 1=alive/survived >=60m\n%d censored excluded",
         n_censored)) +
  theme_minimal(base_size=13) +
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="none")

print(p1)
ggsave(file.path(PATHS$out_brca_figures,"fig07_class_balance.png"),
       p1,width=8,height=6,dpi=300)
cat("Saved: fig07_class_balance.png\n")

# Plot 2: Survival by stage
p2_df <- feature_table_ml |>
  left_join(merged|>select(patient_id,stage_clean),by="patient_id") |>
  filter(stage_clean %in% c("Stage I","Stage II","Stage III","Stage IV")) |>
  mutate(stage_clean=factor(stage_clean,
    levels=c("Stage I","Stage II","Stage III","Stage IV"))) |>
  group_by(stage_clean) |>
  summarise(N=n(), pct_survived=round(mean(five_year_survival)*100,1),
            .groups="drop")

cat("\n--- 5-year survival by stage (ML cohort) ---\n")
print(p2_df)

p2 <- ggplot(p2_df,aes(x=stage_clean,y=pct_survived,fill=stage_clean)) +
  geom_col(width=.65,colour="white") +
  geom_text(aes(label=sprintf("%.1f%%\n(n=%d)",pct_survived,N)),
            vjust=-0.3,size=4.5,fontface="bold") +
  scale_fill_manual(values=c("Stage I"="#5DCAA5","Stage II"="#185FA5",
                              "Stage III"="#F59E0B","Stage IV"="#EF4444")) +
  scale_y_continuous(limits=c(0,115),labels=function(x) paste0(x,"%")) +
  labs(title="Five-Year Survival Rate by Cancer Stage",
       subtitle=sprintf("ML cohort: n=%d confirmed-outcome patients",n_ml),
       x="Cancer Stage",y="% Survived 5 Years",
       caption="Stage differences justify stage_numeric as key ML feature") +
  theme_minimal(base_size=13) +
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="none")

print(p2)
ggsave(file.path(PATHS$out_brca_figures,"fig07_survival_by_stage.png"),
       p2,width=9,height=6,dpi=300)
cat("Saved: fig07_survival_by_stage.png\n")

# Plot 3: Feature distributions
cont_feats <- c("age","stage_numeric","tmb_nonsynonymous",
                "n_nonsense","n_truncating","n_frameshift_del",
                "n_frameshift_ins","n_splice_site")

p3 <- feature_table_ml |>
  select(all_of(cont_feats),five_year_survival) |>
  mutate(Outcome=factor(five_year_survival,
    levels=c(0,1),labels=c("Did NOT survive","Survived 5yr"))) |>
  pivot_longer(-c(five_year_survival,Outcome),
               names_to="Feature",values_to="Value") |>
  ggplot(aes(x=Value,fill=Outcome)) +
  geom_histogram(bins=30,position="identity",alpha=.7,
                 colour="white",linewidth=.1) +
  facet_wrap(~Feature,scales="free",ncol=3) +
  scale_fill_manual(values=c("Did NOT survive"="#EF4444",
                              "Survived 5yr"="#185FA5")) +
  labs(title="Feature Distributions by 5-Year Survival",
       subtitle=sprintf("n=%d | Red=died | Blue=survived",n_ml),
       x="Value",y="Count",fill="5-Year Survival") +
  theme_minimal(base_size=10) +
  theme(plot.background=element_rect(fill="white",colour=NA),
        legend.position="bottom")

print(p3)
ggsave(file.path(PATHS$out_brca_figures,"fig07_feature_distributions.png"),
       p3,width=12,height=9,dpi=300)
cat("Saved: fig07_feature_distributions.png\n")

# Plot 4: Gene frequency
p4 <- gene_freq |>
  mutate(gene=fct_reorder(gene,pct_mutated)) |>
  ggplot(aes(x=gene,y=pct_mutated)) +
  geom_col(fill="#185FA5",width=.75,colour="white") +
  geom_text(aes(label=sprintf("%.1f%%",pct_mutated)),
            hjust=-0.1,size=3.5,fontface="bold") +
  coord_flip() +
  scale_y_continuous(limits=c(0,max(gene_freq$pct_mutated)*1.3),
                     labels=function(x)paste0(x,"%")) +
  labs(title="Top 20 Gene Mutation Frequency — ML Cohort",
       subtitle=sprintf("n=%d confirmed-outcome patients",n_ml),
       x=NULL,y="% Patients Mutated",
       caption="All 20 encoded as binary features (1=mutated, 0=not)") +
  theme_minimal(base_size=12) +
  theme(plot.background=element_rect(fill="white",colour=NA))

print(p4)
ggsave(file.path(PATHS$out_brca_figures,"fig07_gene_frequency.png"),
       p4,width=10,height=8,dpi=300)
cat("Saved: fig07_gene_frequency.png\n")

# ── FINAL SUMMARY ─────────────────────────────────────────────
cat("\n====================================================\n")
cat(" SCRIPT 07 COMPLETE\n")
cat("====================================================\n")
cat(sprintf("  ML dataset  : %d x %d (feature_table_ml.csv)\n",
            nrow(feature_table_ml),ncol(feature_table_ml)))
cat(sprintf("  Cox dataset : %d x %d (feature_table_cox.csv)\n",
            nrow(feature_table_cox),ncol(feature_table_cox)))
cat(sprintf("  Survived    : %d | Died: %d | Balance: %.1f:1\n",
            n_survived, n_died, n_survived/n_died))
cat(sprintf("  Diag cols   : %d (NA column removed)\n",length(diag_cols)))
cat("\nNext: source('scripts/08_train_ml_and_cox_models.R')\n")
cat("====================================================\n")
