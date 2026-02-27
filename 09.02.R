library(dplyr)
library(haven)
library(labelled)
library(stringr)  
library(ggplot2)  
library(tidyr)  
library(readxl)
library(skimr)
library(ggplot2)
library(tidyverse)
library(sf)
library(ggthemes)
library(rnaturalearth)
library(dlookr)
library(scales)
library(patchwork)
library(viridis)
library(tidyverse)
library(gifski)
library(corrplot)
library(fixest)

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

##REGRESSION ANALYIS : Fixed Effects OLS and TWFE
##Start with some basic regressions
## Municipality level variations, Dept level variations first w municipality and dept fixed effects and then with these plus year fixed effects
dept_model1 <- feols(s11_score_total ~ AMPI | dept_code, 
                     data = merged, 
                     cluster = ~dept_code)
summary(dept_model1)

dept_model1.1 <- feols(s11_score_math ~ AMPI | dept_code, 
                     data = merged, 
                     cluster = ~dept_code)
summary(dept_model1.1)



dept_model1.2 <- feols(s11_score_english ~ AMPI | dept_code, 
                     data = merged, 
                     cluster = ~dept_code)
summary(dept_model1)




dept_model2 <- feols(net_cov_edu ~ AMPI | dept_code, 
                     data = merged, 
                     cluster = ~dept_code)
summary(dept_model2)

dept_model3 <- feols(gross_cov_edu ~ AMPI | dept_code, 
                     data = merged, 
                     cluster = ~dept_code)
summary(dept_model3)


dept_model4 <- feols(students_total ~ AMPI | dept_code, 
                     data = merged, 
                     cluster = ~dept_code)
summary(dept_model4)

dept_model5 <- feols(teachers_total ~ AMPI | dept_code, 
                     data = merged, 
                     cluster = ~dept_code)
summary(dept_model5)


dept_model6 <- feols(stud_teach_ratio ~ AMPI | dept_code, 
                     data = merged, 
                     cluster = ~dept_code)
summary(dept_model6)

# Compare all department models
etable(dept_model1,dept_model1.1,dept_model1.2, dept_model2, dept_model3, dept_model4, dept_model5,dept_model6,
       headers = c("Test Scores","Test Scores English", "Test Scores Math" , "Net Coverage", "Gross Coverage", 
                   "Students", "Teachers", "STR"),
       title = "Department Fixed Effects Models")

##Models with both year and location FE
# Model 1: Test Scores
dept_model7 <- feols(s11_score_total ~ AMPI | dept_code + ano, 
                     data = merged, 
                     cluster = ~dept_code)
summary(dept_model7)

dept_model7.1 <- feols(s11_score_english ~ AMPI | dept_code + ano, 
                     data = merged, 
                     cluster = ~dept_code)
summary(dept_model7.1)

dept_model7.2 <- feols(s11_score_math ~ AMPI | dept_code + ano, 
                     data = merged, 
                     cluster = ~dept_code)
summary(dept_model7.2)

# Model 2: Net Coverage
dept_model8 <- feols(net_cov_edu ~ AMPI | dept_code + ano, 
                     data = merged, 
                     cluster = ~dept_code)
summary(dept_model8)
# Model 3: Gross Coverage
dept_model9 <- feols(gross_cov_edu ~ AMPI | dept_code + ano, 
                     data = merged, 
                     cluster = ~dept_code)
summary(dept_model9)

# Model 4: Student-Teacher Ratio
dept_model10 <- feols(stud_teach_ratio ~ AMPI | dept_code + ano, 
                      data = merged, 
                      cluster = ~dept_code)
summary(dept_model10)

# Compare all department models
etable(dept_model7,dept_model7.1,dept_model7.2, dept_model8, dept_model9, dept_model10,
       headers = c("Test Scores","Test Scores English", "Test Scores Math", "Net Coverage", "Gross Coverage", 
                   "STR"),
       title = "Department Fixed Effects Models (with Year FE)")

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

##Running DID model 
models_dept <- feols(
  c(s11_score_total, s11_score_math, s11_score_english, net_cov_edu, gross_cov_edu, stud_teach_ratio) ~ treated_dept * post |
    dept_code + ano,
  cluster = ~dept_code,
  data = merged
)

etable(models_dept)

###CUTOFF - Q75 
##Checking the parallel trends assumptions - dept level - q75 treatment 
outcomes <- c("s11_score_total", "s11_score_math", "s11_score_english", "net_cov_edu", "gross_cov_edu", "stud_teach_ratio") 
trend_data <- merged %>% 
  group_by(ano, treated_dept_q75) %>% 
  summarise( across(all_of(outcomes), ~ mean(.x, na.rm = TRUE)), .groups = "drop" ) %>% 
  pivot_longer(cols = all_of(outcomes), names_to = "outcome", values_to = "mean_value") 

ggplot(trend_data, 
       aes(ano, mean_value, color = factor(treated_dept_q75))) + 
  geom_line(size = 1.1) + 
  geom_vline(xintercept = 2016, linetype = "dashed") + 
  facet_wrap(~ outcome, scales = "free_y") + 
  labs(color = "Treated") + 
  theme_minimal()

models_dept_q75 <- feols(
  c(s11_score_total, s11_score_math, s11_score_english, net_cov_edu, gross_cov_edu, stud_teach_ratio) ~ treated_dept_q75:post | dept_code + ano,
      cluster = ~dept_code,
      data = merged)
etable(models_dept_q75)

##Cutoff Q90
models_dept_q90 <- feols(
  c(s11_score_total, s11_score_math, s11_score_english, net_cov_edu, gross_cov_edu, stud_teach_ratio) ~ treated_dept_q90:post | dept_code + ano,
  cluster = ~dept_code,
  data = merged)
etable(models_dept_q90)

skim(merged) %>% 
  arrange(desc(complete_rate))

# Department level models
dept_models <- feols(
  .[outcomes] ~ AMPI + csw0(rural_i, effect_i, efficien_i, rural_i + effect_i + efficien_i) | dept_code + ano,
  data = merged,
  cluster = ~dept_code
)
etable(dept_models)

# Method 1: Create separate model sets for each control specification
dept_models_no <- feols(
  .[outcomes] ~ AMPI | dept_code + ano,
  data = merged,
  cluster = ~dept_code
) 

dept_models_rural <- feols(
  .[outcomes] ~ AMPI + rural_i | dept_code + ano,
  data = merged,
  cluster = ~dept_code
)

dept_models_effect <- feols(
  .[outcomes] ~ AMPI + effect_i | dept_code + ano,
  data = merged,
  cluster = ~dept_code
)

dept_models_efficien <- feols(
  .[outcomes] ~ AMPI + efficien_i | dept_code + ano,
  data = merged,
  cluster = ~dept_code
)

dept_models_all <- feols(
  .[outcomes] ~ AMPI + rural_i + effect_i + efficien_i | dept_code + ano,
  data = merged,
  cluster = ~dept_code
)

# Now create your 4 tables
table0 <-  etable(dept_models_no, title = "No controls")
table1 <- etable(dept_models_rural, title = "Rural Control Only")
table2 <- etable(dept_models_effect, title = "Effect Control Only")
table3 <- etable(dept_models_efficien, title = "Efficiency Control Only")
table4 <- etable(dept_models_all, title = "All Controls")
# Print tables
table0
table1
table2
table3
table4

colnames(merged)

summary(merged)

#  Tables: 1) Summary Stats | 2) Main No Controls | 3) Main Controls
#          4) Ref Year Sensitivity | 5) Cutoff Sensitivity
vars_summary <- c("s11_score_total", "s11_score_math", "s11_score_english",
                  "net_cov_edu", "gross_cov_edu", "stud_teach_ratio",
                  "AMPI", "rural_i", "effect_i", "efficien_i",
                  "treated_dept", "post")

tbl1 <- datasummary(
  All(merged[, vars_summary]) ~
    N + Mean + SD + Min + Median + Max,
  data    = merged[, vars_summary],
  title   = "Table 1. Summary Statistics")
tbl1

models_dept <- feols(
  c(s11_score_total, s11_score_math, s11_score_english, net_cov_edu, gross_cov_edu, stud_teach_ratio) ~ treated_dept * post |
    dept_code + ano,
  cluster = ~dept_code,
  data = merged
)
etable(models_dept)

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

# ---- 0) Basic setup ----
# df: your panel municipality-year dataset
# Required columns: codmpio, ano, ampi, dept_code, outcomes + controls

# Choose your policy/pivot year (peace agreement implementation)
post_year <- 2016

# Define pre-period used to compute baseline intensity
pre_start <- 2010
pre_end   <- 2015   # last pre year (change if you treat 2015 as "transition")

df2 <- merged %>%
  mutate(
    post = as.integer(ano >= post_year)
  )

# ---- 1) Construct baseline (pre-period) AMPI by municipality ----
ampi_pre <- df2 %>%
  filter(ano >= pre_start, ano <= pre_end) %>%
  group_by(dept_code) %>%
  summarise(ampi_pre = mean(AMPI, na.rm = TRUE), .groups = "drop")

df2 <- df2 %>%
  left_join(ampi_pre, by = "dept_code")

# Optional: standardize intensity so coefficients are "per 1 SD"
df2 <- df2 %>%
  mutate(ampi_pre_z = as.numeric(scale(ampi_pre)))

# ---- 2) Continuous DID regression ----
# Controls: adjust to your exact column names
controls <- c("rural_i", "efficien_i", "effect_i")

# Build formula programmatically
rhs <- paste0("post:ampi_pre_z + ", paste(controls, collapse = " + "))
fml <- as.formula(paste("s11_score_total ~", rhs, "| dept_code + ano"))

m1 <- feols(
  fml,
  data = df2,
  cluster = ~ dept_code  # or ~codmpio; or multiway: cluster = ~codmpio + ano
)

summary(m1)

# Interpretation:
# coef on post:ampi_pre_z = change in outcome in post-period per 1 SD higher baseline conflict
# Choose a reference year in the pre-period (commonly last pre year)
ref_year <- 2015

fml_es <- as.formula(paste(
  "s11_score_total ~ i(ano, ampi_pre_z, ref = ", ref_year, ") + ",
  paste(controls, collapse = " + "),
  " | dept_code + ano"
))

m_es <- feols(fml_es, data = df2, cluster = ~ dept_code)

summary(m_es)

# Plot dynamic effects (per 1 SD higher baseline AMPI, relative to ref_year)
iplot(m_es, ref.line = 0,
      main = "Event-study: Differential change by baseline conflict (1 SD)",
      xlab = "Year", ylab = "Coefficient (vs ref year)")

outcomes <- c("s11_score_total", "s11_score_math", "s11_score_english",
              "net_cov_edu", "gross_cov_edu", "stud_teach_ratio")

models <- lapply(outcomes, function(y){
  fml_y <- as.formula(paste(y, "~", rhs, "| dept_code + ano"))
  feols(fml_y, data = df2, cluster = ~ dept_code)
})

etable(models, se = "cluster")


# Choose a reference year in the pre-period (commonly last pre year)
ref_year <- 2015

fml_es <- as.formula(paste(
  "s11_score_total ~ i(ano, ampi_pre_z, ref = ", ref_year, ") + ",
  paste(controls, collapse = " + "),
  " | dept_code + ano"
))

m_es <- feols(fml_es, data = df2, cluster = ~ dept_code)

summary(m_es)

# Plot dynamic effects (per 1 SD higher baseline AMPI, relative to ref_year)
iplot(m_es, ref.line = 0,
      main = "Event-study: Differential change by baseline conflict (1 SD)",
      xlab = "Year", ylab = "Coefficient (vs ref year)")


m_tv <- feols(
  s11_score_total ~ post:AMPI + rural_i + efficien_i + effect_i |
    dept_code + ano,
  data = df2,
  cluster = ~ dept_code
)
summary(m_tv)



