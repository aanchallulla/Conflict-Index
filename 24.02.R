library(dplyr)
library(haven)
library(labelled)
library(stringr)  
library(ggplot2)  
library(tidyr)  
library(readxl)
library(skimr)
library(tidyverse)
library(sf)
library(ggthemes)
library(rnaturalearth)
library(dlookr)
library(scales)
library(patchwork)
library(plm)
library(fixest)
library(corrplot)
library(gtsummary)
library(gt)
library(tidyr)
library(knitr)
library(kableExtra)

##Load the datasets and combine them 
conflict <- read_csv("conflict_combined.csv")
##Threats, Forced displacement and Homicides - key conflict indicators 
##their pct and AMPI 
edu <- read_csv("edu_combined2.csv")
poverty <- read_csv("poverty.csv")

##filter the codmpios with no data and NAs
conflict <- conflict %>%
  filter(!codmpio %in% 97777)
edu <- edu %>%
  filter(!codmpio %in% 97777,
         !codmpio %in% NA)

##merging 
merged <- edu %>%
  left_join(conflict, by = c("dept_code","codmpio", "ano"))

merged <- merged %>% 
  left_join(poverty,by = c("codmpio", "ano"))

##data that will be used for analysis
merged <- merged %>% 
  select(
    dept_code,
    codmpio,
    ano,
    AMPI,
    s11_score_total,
    s11_score_math,
    s11_score_english,
    gross_cov_edu,
    net_cov_edu,
    stud_teach_ratio,
    rural_i, ##Rurality index (Rural population/Total population)
    rur_student, ##Total school students in rural area
    effect_i, ##Effectiveness index municipality government
    efficien_i ##Efficiency index municipality government
  )

#We only need data from 2010 to 2020
merged <- merged %>% 
  filter(ano >= 2010 & ano <= 2020,
         !dept_code %in% 91)

##Creating an avg AMPI for each dept, so as to use it for the defining the treatment and control groups
##dept level
merged <- merged %>%
  group_by(dept_code) %>% 
  mutate(
    avg_AMPI_dept = mean(AMPI[ano >= 2010 & ano <= 2015], na.rm = TRUE)
  ) %>% 
  ungroup()

##Cutoffs for regressions - dept level

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

##Checking the counts of treated and control depts
merged %>%
  distinct(dept_code, avg_AMPI_dept) %>%
  arrange(avg_AMPI_dept) %>% 
  print(n=32)

##cutoff med - 16 control, 16 treat
merged %>%
  distinct(dept_code, treated_dept) %>%
  count(treated_dept)
##cutoff 75 - 24 control 8 treatment 
merged %>%
  distinct(dept_code, treated_dept_q75) %>%
  count(treated_dept_q75)
##cutoff 90 -28 control 4 treatment 
merged %>%
  distinct(dept_code, treated_dept_q90) %>%
  count(treated_dept_q90)

##Adding a post-treatment dummy 
merged <- merged %>% 
  mutate(post = ifelse(ano > 2016, 1, 0))

##Checking the parallel trend assumptions
###CUT OFF - MEDIAN
##Checking the parallel trends assumptions - dept level 
outcomes <- c("s11_score_total", "s11_score_math", "s11_score_english", "net_cov_edu", "gross_cov_edu", "stud_teach_ratio") 
trend_data <- merged %>% 
  group_by(ano, treated_dept) %>% 
  summarise( across(all_of(outcomes), ~ mean(.x, na.rm = TRUE)), .groups = "drop" ) %>% 
  pivot_longer(cols = all_of(outcomes), names_to = "outcome", values_to = "mean_value") 

ggplot(trend_data, 
       aes(ano, mean_value, color = factor(treated_dept))) + 
  geom_line(size = 1.1) + 
  geom_vline(xintercept = 2016, linetype = "dashed") + 
  facet_wrap(~ outcome, scales = "free_y") + 
  labs(color = "Treated") + 
  theme_minimal()


##Summary statistics
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

tab1 <- merged %>%
  select(all_of(vars_summary)) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "x") %>%
  group_by(Variable) %>%
  summarise(
    Obs.     = sum(!is.na(x)),
    Mean     = mean(x, na.rm = TRUE),
    `Est.Dev.` = sd(x, na.rm = TRUE),
    Min      = min(x, na.rm = TRUE),
    Median = median(x, na.rm = TRUE),
    Max      = max(x, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Variable = recode(Variable, !!!var_labels_summary),
    across(c(Mean, `Est.Dev.`, Min, Median, Max), ~ round(.x, 3))
  )

# 3) Render (LaTeX/PDF look)
kable(tab1, format = "latex", booktabs = TRUE,
      caption = "Descriptive statistics.",
      align = c("l","r","r","r","r","r")) %>%
  kable_styling(latex_options = c("hold_position")) %>%
  row_spec(0, bold = TRUE)


##Running DID model - no controls - median cutoff
models_dept <- feols(
  c(s11_score_total, s11_score_math, s11_score_english, net_cov_edu, gross_cov_edu, stud_teach_ratio) ~ treated_dept * post |
    dept_code + ano,
  cluster = ~dept_code,
  data = merged
)

etable(models_dept)

##Running DID model - all controls - median cutoff
models_dept_controls <- feols(
  c(s11_score_total, s11_score_math, s11_score_english,
    net_cov_edu, gross_cov_edu, stud_teach_ratio) ~
    treated_dept * post +
    rural_i + effect_i + efficien_i |
    dept_code + ano,
  cluster = ~dept_code,
  data = merged
)

etable(models_dept_controls)

##Cutoff 75
models_dept_controls_q75 <- feols(
  c(s11_score_total, s11_score_math, s11_score_english,
    net_cov_edu, gross_cov_edu, stud_teach_ratio) ~
    treated_dept_q75 * post +
    rural_i + effect_i + efficien_i |
    dept_code + ano,
  cluster = ~dept_code,
  data = merged
)

etable(models_dept_controls_q75)

##Cutoff 90
models_dept_controls_q90 <- feols(
  c(s11_score_total, s11_score_math, s11_score_english,
    net_cov_edu, gross_cov_edu, stud_teach_ratio) ~
    treated_dept_q90 * post +
    rural_i + effect_i + efficien_i |
    dept_code + ano,
  cluster = ~dept_code,
  data = merged
)

etable(models_dept_controls_q90)

models_dept_controls_ampi <- feols(
  c(s11_score_total, s11_score_math, s11_score_english,
    net_cov_edu, gross_cov_edu, stud_teach_ratio) ~
    AMPI * post +
    rural_i + effect_i + efficien_i |
    dept_code + ano,
  cluster = ~dept_code,
  data = merged
)

etable(models_dept_controls_ampi)

##Sensitivity Analysis, ref year - 2014
merged <- merged %>% 
  mutate(post_14 = ifelse(ano > 2014, 1, 0))

merged <- merged %>%
  group_by(dept_code) %>% 
  mutate(
    avg_AMPI_dept14 = mean(AMPI[ano >= 2010 & ano <= 2013], na.rm = TRUE)
  ) %>% 
  ungroup()


cutoff_dept14 <- merged %>%
  distinct(dept_code, avg_AMPI_dept14) %>%
  summarise(med13 = median(avg_AMPI_dept14, na.rm = TRUE)) %>%
  pull(med13)

merged <- merged %>%
  mutate(treated_dept14 = ifelse(avg_AMPI_dept14 > cutoff_dept14, 1, 0))

models_dept_controls_14 <- feols(
  c(s11_score_total, s11_score_math, s11_score_english,
    net_cov_edu, gross_cov_edu, stud_teach_ratio) ~
    treated_dept14 * post_14 +
    rural_i + effect_i + efficien_i |
    dept_code + ano,
  cluster = ~dept_code,
  data = merged
)

etable(models_dept_controls_14)


##Sensitivity Analysis, ref year 2015
merged <- merged %>% 
  mutate(post_15 = ifelse(ano > 2015, 1, 0))

merged <- merged %>%
  group_by(dept_code) %>% 
  mutate(avg_AMPI_dept15 = mean(AMPI[ano >= 2010 & ano <= 2014], na.rm = TRUE)
  ) %>% 
  ungroup()

cutoff_dept15 <- merged %>%
  distinct(dept_code, avg_AMPI_dept15) %>%
  summarise(med15 = median(avg_AMPI_dept15, na.rm = TRUE)) %>%
  pull(med15)

merged <- merged %>%
  mutate(treated_dept15 = ifelse(avg_AMPI_dept15 > cutoff_dept15, 1, 0))

models_dept_controls_15 <- feols(
  c(s11_score_total, s11_score_math, s11_score_english,
    net_cov_edu, gross_cov_edu, stud_teach_ratio) ~
    treated_dept15 * post_15 +
    rural_i + effect_i + efficien_i |
    dept_code + ano,
  cluster = ~dept_code,
  data = merged
)

etable(models_dept_controls_15)

##Continous DID - as robustness?
models_dept_conti <- feols(
  c(s11_score_total, s11_score_math, s11_score_english,
    net_cov_edu, gross_cov_edu, stud_teach_ratio) ~
    AMPI * post +
    rural_i + effect_i + efficien_i |
    dept_code + ano,
  cluster = ~dept_code,
  data = merged
)

etable(models_dept_conti)


models_dynamic <- feols(
  c(s11_score_total, s11_score_math, s11_score_english,
    net_cov_edu, gross_cov_edu, stud_teach_ratio) ~ 
    i(ano, treated_dept, ref = 2015) + 
    rural_i + effect_i + efficien_i | 
    dept_code + ano,
  cluster = ~dept_code,
  data = merged
)

etable(models_dynamic)

