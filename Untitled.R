# ============================================================
#  THESIS TABLES: Impact of Conflict/Peace on Education in Colombia
#  DID / TWFE Analysis
#  Tables: 1) Summary Stats | 2) Main No Controls | 3) Main Controls
#          4) Ref Year Sensitivity | 5) Cutoff Sensitivity
# ============================================================

# ---- 0. Packages -------------------------------------------
library(tidyverse)
library(fixest)       # TWFE models
library(modelsummary) # Table output
library(kableExtra)   # Table formatting
library(haven)        # If loading .dta
library(gt)           # Optional: pretty tables

# ---- 0.1 Load Data -----------------------------------------
##load the dataset
edu <- read_csv("edu_combined2.csv")
conflict <- read_csv("conflict_combined.csv")
##removing NA codmpios and municipalities that don't have any data
conflict <- conflict %>%
  filter(!codmpio %in% 97777)
edu <- edu %>%
  filter(!codmpio %in% 97777,
         !codmpio %in% NA)

##merging the datasets
merged <- edu %>%
  left_join(conflict, by = c("dept_code","codmpio", "ano"))
colnames(merged)
##using data from 2010-2020 and selecting only the variables i need
merged <- merged %>% 
  filter(ano >= 2010 & ano <= 2020,
         !dept_code %in% 91)

merged <- merged %>% 
  select(
    dept_code,
    codmpio,
    ano,
    students_total,
    teachers_total,
    s11_score_total,
    s11_score_math,
    s11_score_english,
    net_cov_edu,
    gross_cov_edu,
    stud_teach_ratio,
    AMPI,
    rural_i, ##Rurality index (Rural population/Total population)
    rur_student, ##Total school students in rural area
    effect_i, ##Effectiveness index municipality government
    efficien_i ##Efficiency index municipality government
  )

# ---- 0.2 Shared Settings -----------------------------------
outcomes     <- c("s11_score_total", "s11_score_math", "s11_score_english",
                  "net_cov_edu", "gross_cov_edu", "stud_teach_ratio")

outcome_labels <- c("Total Score", "Math Score", "English Score",
                    "Net Coverage", "Gross Coverage", "Student-Teacher Ratio")

controls <- c("rural_i", "effect_i", "efficien_i")
ctrl_formula <- paste(controls, collapse = " + ")

# Clean coefficient/row names for all tables
coef_map <- c(
  "treated_dept:post"     = "DID (Treated × Post)",
  "treated_dept_q75:post" = "DID (Treated × Post)",
  "treated_dept_q90:post" = "DID (Treated × Post)",
  "rural_i"               = "Rural",
  "effect_i"              = "Effectiveness",
  "efficien_i"            = "Efficiency"
)

gof_keep <- c("nobs", "r.squared", "adj.r.squared")

gof_map <- list(
  list(raw = "nobs",         clean = "Observations",  fmt = 0),
  list(raw = "r.squared",    clean = "R²",             fmt = 3),
  list(raw = "adj.r.squared",clean = "Adj. R²",        fmt = 3)
)

# ============================================================
# TABLE 1 — SUMMARY STATISTICS
# ============================================================

vars_summary <- c("s11_score_total", "s11_score_math", "s11_score_english",
                  "net_cov_edu", "gross_cov_edu", "stud_teach_ratio",
                  "AMPI", "rural_i", "effect_i", "efficien_i",
                  "treated_dept", "post")

var_labels_summary <- c(
  "s11_score_total"   = "Total Score (Saber 11)",
  "s11_score_math"    = "Math Score (Saber 11)",
  "s11_score_english" = "English Score (Saber 11)",
  "net_cov_edu"       = "Net Education Coverage",
  "gross_cov_edu"     = "Gross Education Coverage",
  "stud_teach_ratio"  = "Student-Teacher Ratio",
  "AMPI"              = "AMPI (Conflict Index)",
  "rural_i"           = "Rural Index",
  "effect_i"          = "Effectiveness Index",
  "efficien_i"        = "Efficiency Index",
  "treated_dept"      = "Treated (Median Cutoff)",
  "post"              = "Post-Treatment"
)

##Creating an avg AMPI for each dept, so as to use it for the defining the treatment and control groups
##dept level
merged <- merged %>%
  group_by(dept_code) %>% 
  mutate(
    avg_AMPI_dept = mean(AMPI[ano >= 2010 & ano <= 2015], na.rm = TRUE)
  ) %>% 
  ungroup()

cutoff_dept <- merged %>%
  distinct(dept_code, avg_AMPI_dept) %>%
  summarise(med = median(avg_AMPI_dept, na.rm = TRUE)) %>%
  pull(med)

cutoff_dept_q75 <- merged %>%
  distinct(dept_code, avg_AMPI_dept) %>%   # one value per municipality
  summarise(q75 = quantile(avg_AMPI_dept, 0.75, na.rm = TRUE)) %>%
  pull(q75)

cutoff_dept_q90 <- merged %>%
  distinct(dept_code, avg_AMPI_dept) %>%   # one value per municipality
  summarise(q90 = quantile(avg_AMPI_dept, 0.90, na.rm = TRUE)) %>%
  pull(q90)

##Assigning treatment and control depts for each cutoff level 
merged <- merged %>%
  mutate(treated_dept = ifelse(avg_AMPI_dept > cutoff_dept, 1, 0))

merged <- merged %>%
  mutate(treated_dept_q75 = ifelse(avg_AMPI_dept > cutoff_dept_q75, 1, 0))

merged <- merged %>%
  mutate(treated_dept_q90 = ifelse(avg_AMPI_dept > cutoff_dept_q90, 1, 0))

merged %>%
  distinct(dept_code, avg_AMPI_dept) %>%
  arrange(avg_AMPI_dept) %>% 
  print(n=32)

merged %>%
  distinct(dept_code, treated_dept) %>%
  count(treated_dept)

merged %>%
  distinct(dept_code, treated_dept_q75) %>%
  count(treated_dept_q75)

merged %>%
  distinct(dept_code, treated_dept_q90) %>%
  count(treated_dept_q90)
##16 control, 16 treatment 
merged %>%
  distinct(dept_code, treated_dept) %>%
  filter(treated_dept == 1)
##24 control 8 treatment 
merged %>%
  distinct(dept_code, treated_dept_q75) %>%
  filter(treated_dept_q75 == 1)
##28 control 4 treatment 
merged %>%
  distinct(dept_code, treated_dept_q90) %>%
  filter(treated_dept_q90 == 1)

##Adding a post-treatment dummy 
merged <- merged %>% 
  mutate(post = ifelse(ano > 2016, 1, 0))

tbl1 <- datasummary(
  All(merged[, vars_summary]) ~
    N + Mean + SD + Min + Median + Max,
  data    = merged[, vars_summary],
  title   = "Table 1. Summary Statistics",
  output  = "kableExtra"
) %>%
  kable_styling(latex_options = c("hold_position", "scale_down"),
                font_size = 10) %>%
  pack_rows("Outcome Variables", 1, 6) %>%
  pack_rows("Conflict Index",    7, 7) %>%
  pack_rows("Controls",          8, 10) %>%
  pack_rows("DID Variables",     11, 12)

# Save Table 1
tbl1 %>% save_kable("table1_summary_stats.pdf")
cat("Table 1 saved.\n")


# ============================================================
# TABLE 2 — MAIN DID RESULTS: TWFE, Median Cutoff, NO CONTROLS
# ============================================================

# Run TWFE for each outcome — median cutoff, no controls
models_main_nc <- lapply(outcomes, function(y) {
  fml <- as.formula(paste0(y, " ~ treated_dept:post | dept_code + ano"))
  feols(fml, data = merged, cluster = ~dept_code)
})
names(models_main_nc) <- outcome_labels

tbl2 <- modelsummary(
  models_main_nc,
  coef_map    = coef_map,
  gof_map     = gof_map,
  stars       = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  title       = "Table 2. Main DID Results — TWFE, Median Cutoff (No Controls)",
  notes       = "Standard errors clustered at department level. Department and year fixed effects included in all specifications.",
  output      = "kableExtra"
) %>%
  kable_styling(latex_options = c("hold_position", "scale_down"),
                font_size = 10) %>%
  add_header_above(c(" " = 1,
                     "Test Scores" = 3,
                     "Coverage" = 2,
                     "Teacher Supply" = 1))

tbl2 %>% save_kable("table2_main_no_controls.pdf")
cat("Table 2 saved.\n")


# ============================================================
# TABLE 3 — MAIN DID RESULTS: TWFE, Median Cutoff, WITH CONTROLS
# ============================================================

models_main_wc <- lapply(outcomes, function(y) {
  fml <- as.formula(paste0(y, " ~ treated_dept:post + ", ctrl_formula,
                           " | dept_code + ano"))
  feols(fml, data = merged, cluster = ~dept_code)
})
names(models_main_wc) <- outcome_labels

tbl3 <- modelsummary(
  models_main_wc,
  coef_map    = coef_map,
  gof_map     = gof_map,
  stars       = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  title       = "Table 3. Main DID Results — TWFE, Median Cutoff (With Controls)",
  notes       = "Standard errors clustered at department level. Department and year fixed effects included. Controls: rural index, effectiveness, and efficiency.",
  output      = "kableExtra"
) %>%
  kable_styling(latex_options = c("hold_position", "scale_down"),
                font_size = 10) %>%
  add_header_above(c(" " = 1,
                     "Test Scores" = 3,
                     "Coverage" = 2,
                     "Teacher Supply" = 1))

tbl3 %>% save_kable("table3_main_with_controls.pdf")
cat("Table 3 saved.\n")


# ============================================================
# TABLE 4 — REFERENCE YEAR SENSITIVITY (2014, 2015 vs 2016 baseline)
# ============================================================

# Helper: run TWFE for a given reference year and extract DID coef + SE
run_ref_year <- function(ref_year, outcome, with_controls = TRUE) {
  df_sub <- merged %>% filter(ano != ref_year | ano == ref_year)
  # Re-define post relative to reference year
  # Assumes peace agreement year is known — adjust cutoff_year as needed
  cutoff_year <- 2016  # year the peace agreement took effect
  df_sub <- df_sub %>%
    mutate(post_alt = ifelse(ano >= cutoff_year, 1, 0))
  
  ctrl_part <- if (with_controls) paste0(" + ", ctrl_formula) else ""
  fml <- as.formula(paste0(outcome, " ~ treated_dept:post_alt",
                           ctrl_part, " | dept_code + ano"))
  
  # Drop the reference year from estimation
  df_est <- df_sub %>% filter(ano != ref_year)
  
  feols(fml, data = df_est, cluster = ~dept_code)
}

ref_years  <- c(2014, 2015)
year_labels <- c("Ref. Year: 2014", "Ref. Year: 2015")

# For each outcome, run both reference years (with controls)
models_ref <- lapply(ref_years, function(ry) {
  lapply(outcomes, function(y) {
    run_ref_year(ref_year = ry, outcome = y, with_controls = TRUE)
  })
})

# Flatten into named list: alternating by ref year within each outcome
# Structure: for each outcome, two columns (2014, 2015)
models_ref_flat <- list()
for (i in seq_along(outcomes)) {
  models_ref_flat[[paste0(outcome_labels[i], " (2014)")]] <- models_ref[[1]][[i]]
  models_ref_flat[[paste0(outcome_labels[i], " (2015)")]] <- models_ref[[2]][[i]]
}

tbl4 <- modelsummary(
  models_ref_flat,
  coef_map    = coef_map,
  gof_map     = gof_map,
  stars       = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  title       = "Table 4. Sensitivity Analysis — Reference Year (2014 and 2015)",
  notes       = "Baseline reference year is 2016. All specifications include department and year FEs, all controls, clustered SEs at department level.",
  output      = "kableExtra"
) %>%
  kable_styling(latex_options = c("hold_position", "scale_down"),
                font_size = 9) %>%
  add_header_above(c(" "  = 1,
                     "Total Score"   = 2, "Math Score"  = 2,
                     "English Score" = 2, "Net Coverage" = 2,
                     "Gross Coverage"= 2, "Stud-Teacher" = 2))

tbl4 %>% save_kable("table4_ref_year_sensitivity.pdf")
cat("Table 4 saved.\n")


# ============================================================
# TABLE 5 — CUTOFF SENSITIVITY (Median vs Q75 vs Q90)
# ============================================================

cutoff_vars   <- c("treated_dept", "treated_dept_q75", "treated_dept_q90")
cutoff_labels <- c("Median", "Q75", "Q90")

# For each cutoff × outcome, run TWFE with controls
models_cutoff <- lapply(cutoff_vars, function(tv) {
  lapply(outcomes, function(y) {
    fml <- as.formula(paste0(y, " ~ ", tv, ":post + ", ctrl_formula,
                             " | dept_code + ano"))
    feols(fml, data = merged, cluster = ~dept_code)
  })
})

# Flatten: for each outcome, three columns (median, q75, q90)
coef_map_cutoff <- c(
  "treated_dept:post"     = "DID Coefficient",
  "treated_dept_q75:post" = "DID Coefficient",
  "treated_dept_q90:post" = "DID Coefficient"
)

models_cutoff_flat <- list()
for (i in seq_along(outcomes)) {
  for (j in seq_along(cutoff_vars)) {
    models_cutoff_flat[[paste0(outcome_labels[i], " (", cutoff_labels[j], ")")]] <-
      models_cutoff[[j]][[i]]
  }
}

tbl5 <- modelsummary(
  models_cutoff_flat,
  coef_map    = coef_map_cutoff,
  gof_map     = gof_map,
  stars       = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  title       = "Table 5. Sensitivity Analysis — Treatment Cutoff Definition",
  notes       = "Treatment defined by Adjusted Mazziota-Pareto Index (AMPI) at median, 75th, and 90th percentiles. All specifications include department and year FEs and all controls. SEs clustered at department level.",
  output      = "kableExtra"
) %>%
  kable_styling(latex_options = c("hold_position", "scale_down"),
                font_size = 9) %>%
  add_header_above(c(" "  = 1,
                     "Total Score"   = 3, "Math Score"    = 3,
                     "English Score" = 3, "Net Coverage"  = 3,
                     "Gross Coverage"= 3, "Stud-Teacher"  = 3))

tbl5 %>% save_kable("table5_cutoff_sensitivity.pdf")
cat("Table 5 saved.\n")
