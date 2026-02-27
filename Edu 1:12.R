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


### LOADING THE DATASET
### Variables available students_total,teachers_total,students_official,teachers_official,students_preschool,students_middle,students_secondary,
### students_rural,teachers_rural,students_urban,teachers_urban,schools_total,s11_score_total,s11_score_total_female, s11_score_total_male
los_andes <-  read_csv("edu_clean.csv")
str(los_andes)

### Variables available Net coverage in education - Total , Gross coverage in education - Total
terridata <- read_csv("terridata_clean.csv")
str(terridata)

### Variables available     entity_code, education_income, population, education_index, rural_index
efd <- read_csv("efd_clean.csv")
str(efd)

##Making it ready to be merged
los_andes <- los_andes %>% 
  rename(
    codmpio = municipality_code,
    ano = year
  )

terridata <- terridata %>% 
  rename(
    codmpio = 'Entity Code',
    ano = Year,
    net_cov_edu = `Net coverage in education - Total`,
    gross_cov_edu = `Gross coverage in education - Total` 
  )

terridata <- terridata %>% 
  mutate(
    codmpio = as.numeric(codmpio)
  )

efd <- efd %>% 
  rename(
    codmpio = entity_code,
    ano = year,
  )

efd <- efd %>% 
  mutate(
    codmpio = as.numeric(codmpio)
  )

##Merging all datasets
edu <- los_andes %>%
  left_join(efd, by = c("codmpio", "ano")) %>%
  left_join(terridata, by = c("codmpio", "ano"))

str(edu)

##Filtering the data to have years from 1993 - 2020
edu <- edu %>% 
  filter(
    ano >= 1993 & ano <= 2020,
    !codmpio %in% 27086
  )

skim(edu)

###Starting with the exploratory analysis

n_distinct(edu$codmpio) ##1124 muncipalities
table(edu$ano)  # observations per year
edu %>% 
  count(codmpio) %>% 
  summary()  # obs per municipality

##Grouping s11 scores by year - we only have data from 2000 - 2020
edu_year <- edu %>%
  group_by(ano) %>%
  summarise(
    mean_s11 = mean(s11_score_total, na.rm = TRUE),
    median_s11 = median(s11_score_total, na.rm = TRUE),
    sd_s11 = sd(s11_score_total, na.rm = TRUE)
  )

ggplot(edu_year, aes(x = ano, y = mean_s11)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Average s11 score over time",
    x = "Year", y = "Mean s11 score"
  )



##Grouping s11 scores by year - we only have data from 2000 - 2020
s11_year <- edu %>%
  group_by(ano) %>%
  summarise(
    mean_s11 = mean(s11_score_total, na.rm = TRUE),
    median_s11 = median(s11_score_total, na.rm = TRUE),
    sd_s11 = sd(s11_score_total, na.rm = TRUE)
  )

ggplot(s11_year, aes(x = ano, y = mean_s11)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Average s11 score over time",
    x = "Year", y = "Mean s11 score"
  )


##Grouping student totals by year
student_year <- edu %>%
  group_by(ano) %>%
  summarise(
    mean_stud = mean(students_total, na.rm = TRUE),
    median_stud = median(students_total, na.rm = TRUE),
    sd_stud = sd(students_total, na.rm = TRUE)
  )

ggplot(student_year, aes(x = ano, y = mean_stud)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Average students over time",
    x = "Year", y = "Mean Students"
  )


##Grouping gross coverage totals by year
gross_cov_year <- edu %>%
  group_by(ano) %>%
  summarise(
    mean_gross_cov = mean(gross_cov_edu, na.rm = TRUE),
    median_gross_cov = median(gross_cov_edu, na.rm = TRUE),
    sd_gross_cov = sd(gross_cov_edu, na.rm = TRUE)
  )

ggplot(gross_cov_year, aes(x = ano, y = mean_gross_cov)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Average Gross Cov Edu over time",
    x = "Year", y = "Mean Gross Cov Edu"
  )

##Grouping gross coverage totals by year
net_cov_year <- edu %>%
  group_by(ano) %>%
  summarise(
    mean_net_cov = mean(net_cov_edu, na.rm = TRUE),
    median_net_cov = median(net_cov_edu, na.rm = TRUE),
    sd_net_cov = sd(net_cov_edu, na.rm = TRUE)
  )

ggplot(net_cov_year, aes(x = ano, y = mean_net_cov)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Average Net Cov Edu over time",
    x = "Year", y = "Mean Net Cov Edu"
  )





