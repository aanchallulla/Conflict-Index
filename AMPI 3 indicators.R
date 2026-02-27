library(tidyverse)
library(haven)
library(labelled)
library(stringr)  
library(tidyr)  
library(readxl)
library(skimr)
library(ggplot2)
library(sf)
library(ggthemes)
library(rnaturalearth)
library(viridis)
library(gganimate)
library(gifski)
library(corrplot)


##Loading the dataset and choosing event variables
conflict <- read_dta("PANEL_CONFLICTO_Y_VIOLENCIA(2024).dta")
conflict <-  conflict %>% 
  select(
    codmpio,
    ano,
    e_amenaza, ##threats events
    e_desplaza, ##forced displacement 
    e_homic ##homicide 
  )
str(conflict)

##MPI 
##Selecting the indicator variables for the conflict index
indicator_vars <- conflict %>% 
  select(-codmpio, -ano) %>% 
  names()

##Checking the correlation between the indicators
cor_matrix <- cor(conflict[, indicator_vars], use = "pairwise.complete.obs")
print(round(cor_matrix, 3))

##Visualising the indicators 
corrplot(cor_matrix, 
         method = "color",
         type = "upper",
         tl.col = "black",
         tl.cex = 0.7,
         tl.srt = 45,
         addCoef.col = "black",
         number.cex = 0.5,
         col = colorRampPalette(c("#6D9EC1", "white", "#E46726"))(200),
         title = "Correlation Matrix - Conflict Indicators",
         mar = c(0,0,2,0))

corrplot::corrplot(cor_matrix)

##Transforming the conflict indicators into their standardized z scores
conflict_standardized <- conflict %>%
  mutate(across(all_of(indicator_vars),
                ~(100 + 10 * ((. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE))),
                .names = "std_{.col}"))

##Creating a vector of the new standardized variables
std_vars <- paste0("std_", indicator_vars)

##Final step of creating the MAzziotta-Pareto Index
conflict_mpi <- conflict_standardized %>%
  mutate(
    #Mean of standardized indicators for all standardized indicators for each municipal-year unit 
    M = rowMeans(select(., all_of(std_vars)), na.rm = TRUE),    
    #Standard deviation for all standardized indicators for each municipal-year unit 
    S = apply(select(., all_of(std_vars)), 1,
              function(x) sd(x, na.rm = TRUE)),
    
    #coefficient of variation for all standardized indicators for each municipal-year unit  
    ##(sd/m)
    CV = (S / M) * 100     
  )

##Calculating the MPI using the formula
conflict_mpi <- conflict_mpi %>%
  mutate(
    MPI = M + (S * CV / 100), 
    MPI_normalized = ((MPI - min(MPI, na.rm = TRUE)) /
                        (max(MPI, na.rm = TRUE) - min(MPI, na.rm = TRUE))) * 100
  )
conflict_final <- conflict_mpi %>%
  select(codmpio, ano, M, CV, MPI, MPI_normalized)

summary(conflict_final %>% select(M, CV, MPI,MPI_normalized))

##Histogram of the MPI
hist(conflict_final$MPI, 
     main = "Distribution of Mazziotta-Pareto Index",
     xlab = "MPI Value",
     col = "lightblue")

conflict_final %>%
  arrange(desc(MPI)) %>%
  head(10)

##HEATMAP FOR TOP 30 CONFLICT AFFECTED MUNCIPALITIES
##creating a vector of municipalities with high levels of conflict 
top_municipalities <- conflict_final %>%
  group_by(codmpio) %>%
  summarise(avg_mpi = mean(MPI_normalized, na.rm = TRUE)) %>%
  arrange(desc(avg_mpi)) %>%
  head(30) %>%
  pull(codmpio)

##heatmap - top 30 muncipalities
p1 <- conflict_final %>%
  filter(codmpio %in% top_municipalities) %>%
  ggplot(aes(x = ano, y = as.factor(codmpio), fill = MPI_normalized)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_gradient2(low = "#2ca02c", mid = "#ff7f0e", high = "#d62728",
                       midpoint = 50, name = "MPI\n(0-100)") +
  labs(title = "Conflict Intensity Heatmap",
       subtitle = "Top 30 Most Affected Municipalities",
       x = "Year",
       y = "Municipality Code") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7))

print(p1)


##BOX PLOT - CONFLICT INTENSITY BY YEAR 
p2 <- conflict_final %>%
  ggplot(aes(x = as.factor(ano), y = MPI_normalized)) +
  geom_boxplot(fill = "#1f77b4", alpha = 0.7) +
  labs(title = "Distribution of Conflict Intensity by Year",
       x = "Year",
       y = "MPI (Normalized 0-100)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p2)


###MPI by year 
##standardizing the indicators of conflict after grouping them by year 
conflict_mpi_year <- conflict %>%
  group_by(ano) %>%
  mutate(
    across(all_of(indicator_vars),
           ~(100 + 10 * ((. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE))),
           .names = "std_{.col}")
  ) %>%
  ungroup()

std_vars <- paste0("std_", indicator_vars)


conflict_mpi_year <- conflict_mpi_year %>%
  mutate(
    M = rowMeans(select(., all_of(std_vars)), na.rm = TRUE),
    S = apply(select(., all_of(std_vars)), 1,
              function(x) sd(x, na.rm = TRUE)),
    CV = (S/M)*100   
  )


conflict_mpi_year <- conflict_mpi_year %>%
  mutate(
    MPI = M + (S * CV/100)
  )

conflict_mpi_year <- conflict_mpi_year %>%
  group_by(ano) %>%
  mutate(
    MPI_normalized = ((MPI - min(MPI, na.rm = TRUE)) /
                        (max(MPI, na.rm = TRUE) - min(MPI, na.rm = TRUE))) * 100
  ) %>%
  ungroup()

conflict_final_year <- conflict_mpi_year %>%
  select(codmpio, ano, M, CV, MPI, MPI_normalized)

temporal_summary <- conflict_final_year %>%
  group_by(ano) %>%
  summarise(
    mean_MPI = mean(MPI, na.rm = TRUE),
    median_MPI = median(MPI, na.rm = TRUE),
    sd_MPI = sd(MPI, na.rm = TRUE),
    mean_CV = mean(CV, na.rm = TRUE),
    n_municipalities = n(),
    high_conflict_munic = sum(MPI_normalized > 75, na.rm = TRUE)
  )
print(temporal_summary)

##How many municipalities have an MPI of over 75 from 1993 - 2023
p3 <- ggplot(temporal_summary, aes(x = ano, y = high_conflict_munic)) +
  geom_bar(stat = "identity", fill = "#d62728", alpha = 0.8) +
  geom_text(aes(label = high_conflict_munic), vjust = -0.5, size = 3) +
  labs(title = "High-Conflict Municipalities Over Time",
       subtitle = "Number of municipalities with MPI > 75",
       x = "Year",
       y = "Count of Municipalities") +
  theme_minimal()
print(p3)

p4 <- ggplot(temporal_summary, aes(x = ano, y = mean_CV)) +
  geom_line(color = "#ff7f0e", size = 1.2) +
  geom_point(color = "#ff7f0e", size = 3) +
  labs(title = "Average Imbalance in Conflict Types Over Time",
       subtitle = "Higher CV indicates more uneven distribution across conflict types",
       x = "Year",
       y = "Average Coefficient of Variation") +
  theme_minimal()
print(p4) 
##This graph tells us about the average CV per year : 
##it is the avg imbalance/disparity of conflict across municipalities in that year
##High CV : uneven conflict profile
##Low CV : balanced conflict profile



##Adjusted Mazziotta Pareto Index : AMPI
##The polarity is positive since, high homicides = more conflict
##The penalty represents the imbalance between indicators
##I will use the formula for positive polarity and thus subtract the penalty term

calculate_ampi <- function(conflict, indicator_vars, 
                           year_var = "ano",
                           use_reference = FALSE,
                           reference_year = NULL) {
  
  # Get unique years and process them chronologically
  years <- unique(conflict[[year_var]])
  years <- years[order(years)]
  
  # Initialize result columns
  ##Creating vectors for each obs
  n_obs <- nrow(conflict)
  AMPI <- numeric(n_obs)
  M_ri <- numeric(n_obs)
  S_ri <- numeric(n_obs)
  cv_ri <- numeric(n_obs)
  penalty <- numeric(n_obs)
  
  # Store all normalized matrices for output
  ##One R matrix per year
  ##Goalposts list will hold goalposts for each indicator
  R_list <- vector("list", length(years))
  goalposts_list <- vector("list", length(years))
  
  # Process each year separately
  for (year_idx in seq_along(years)) {
    year <- years[year_idx]
    
    # Get indices for this year
    year_indices <- which(conflict[[year_var]] == year)
    
    # Extract indicator matrix for this year
    X_year <- conflict %>% 
      filter(!!sym(year_var) == year) %>%
      select(all_of(indicator_vars))
    X_matrix <- as.matrix(X_year)
    
    # Step 1: Calculate goalposts for this year
    goalposts <- vector("list", ncol(X_matrix))
    
    for (j in seq_len(ncol(X_matrix))) {
      x_j <- X_matrix[, j]
      
      # Check if there are any non-NA values
      if (all(is.na(x_j))) {
        warning(paste("All values are NA for indicator", indicator_vars[j], 
                      "in year", year))
        goalposts[[j]] <- list(Min = NA, Max = NA)
        next
      }
      
      Inf_xj <- min(x_j, na.rm = TRUE)
      Sup_xj <- max(x_j, na.rm = TRUE)
      
      if (use_reference && !is.null(reference_year)) {
        # Reference based on this specific year's data
        Ref_xj <- mean(x_j, na.rm = TRUE)
        Delta <- (Sup_xj - Inf_xj) / 2
        Min_j <- Ref_xj - Delta
        Max_j <- Ref_xj + Delta
      } else {
        Min_j <- Inf_xj
        Max_j <- Sup_xj
      }
      
      goalposts[[j]] <- list(Min = Min_j, Max = Max_j)
    }
    
    # Step 2: Normalize to 70–130 scale for this year
    ##Create a matrix R filled with NAs
    R <- matrix(NA, nrow = nrow(X_matrix), ncol = ncol(X_matrix))
    
    for (j in seq_len(ncol(X_matrix))) {
      x_j <- X_matrix[, j]
      Min_j <- goalposts[[j]]$Min
      Max_j <- goalposts[[j]]$Max
      
      # Skip if goalposts are NA
      if (is.na(Min_j) || is.na(Max_j)) {
        R[, j] <- NA
        next
      }
      
      if (Max_j > Min_j) {
        R[, j] <- 60 * ((x_j - Min_j) / (Max_j - Min_j)) + 70
      } else {
        # If no variation, set to midpoint (but keep NAs as NA)
        R[, j] <- ifelse(is.na(x_j), NA, 100)
      }
    }
    
    # Step 3: Mean and standard deviation for this year
    # Only calculate if there are at least some non-NA values
    M_ri_year <- rowMeans(R, na.rm = TRUE)
    S_ri_year <- apply(R, 1, sd, na.rm = TRUE)
    
    # Handle cases where all values are NA
    M_ri_year[is.nan(M_ri_year)] <- NA
    S_ri_year[is.na(S_ri_year)] <- 0
    
    # Calculate CV
    cv_ri_year <- ifelse(M_ri_year > 0, S_ri_year / M_ri_year, 0)
    cv_ri_year[!is.finite(cv_ri_year)] <- 0
    cv_ri_year[is.na(M_ri_year)] <- NA
    
    # Step 4: Penalty and AMPI for this year
    penalty_year <- S_ri_year * cv_ri_year
    penalty_year[is.na(penalty_year)] <- 0
    
    AMPI_year <- M_ri_year - penalty_year
    
    # Store results in the corresponding positions
    AMPI[year_indices] <- AMPI_year
    M_ri[year_indices] <- M_ri_year
    S_ri[year_indices] <- S_ri_year
    cv_ri[year_indices] <- cv_ri_year
    penalty[year_indices] <- penalty_year
    
    # Store normalized matrix and goalposts
    R_list[[year_idx]] <- R
    goalposts_list[[year_idx]] <- goalposts
    names(R_list)[year_idx] <- as.character(year)
    names(goalposts_list)[year_idx] <- as.character(year)
  }
  
  # Output
  return(list(
    AMPI = AMPI,
    M_ri = M_ri,
    S_ri = S_ri,
    cv_ri = cv_ri,
    penalty = penalty,
    R_by_year = R_list,
    goalposts_by_year = goalposts_list
  ))
}

# Usage example:
ampi_results <- calculate_ampi(conflict, indicator_vars, 
                               year_var = "ano",
                               use_reference = FALSE,
                               reference_year = NULL)

# Add results to the dataframe
conflict_with_ampi <- conflict %>%
  mutate(
    AMPI = ampi_results$AMPI,
    AMPI_mean = ampi_results$M_ri,
    AMPI_penalty = ampi_results$penalty,
    AMPI_cv = ampi_results$cv_ri
  )

print(summary(conflict_with_ampi$AMPI))

#Check year-by-year statistics
ampi_summary_by_year <- conflict_with_ampi %>%
  group_by(ano) %>%
  summarise(
    n_obs = n(),
    n_complete = sum(!is.na(AMPI)),
    mean_AMPI = mean(AMPI, na.rm = TRUE),
    sd_AMPI = sd(AMPI, na.rm = TRUE),
    min_AMPI = min(AMPI, na.rm = TRUE),
    max_AMPI = max(AMPI, na.rm = TRUE)
  ) %>%
  ungroup()

# Diagnostic code to understand AMPI results with NA handling

diagnostic_check <- function(conflict_with_ampi, ampi_results, indicator_vars) {
  
  cat("=== DATA COMPLETENESS CHECK ===\n")
  cat("Total observations:", nrow(conflict_with_ampi), "\n")
  cat("Observations with AMPI:", sum(!is.na(conflict_with_ampi$AMPI)), "\n")
  cat("Observations with NA AMPI:", sum(is.na(conflict_with_ampi$AMPI)), "\n\n")
  
  # Check NAs in original indicators
  cat("=== NA COUNT BY INDICATOR ===\n")
  for (var in indicator_vars) {
    na_count <- sum(is.na(conflict_with_ampi[[var]]))
    na_pct <- round(100 * na_count / nrow(conflict_with_ampi), 1)
    cat(var, ":", na_count, "NAs (", na_pct, "%)\n")
  }
  cat("\n")
  
  cat("=== AMPI RANGE CHECK (excluding NAs) ===\n")
  cat("AMPI range:", range(conflict_with_ampi$AMPI, na.rm = TRUE), "\n")
  cat("Mean range:", range(conflict_with_ampi$AMPI_mean, na.rm = TRUE), "\n")
  cat("Penalty range:", range(conflict_with_ampi$AMPI_penalty, na.rm = TRUE), "\n\n")
  
  # Find the observation with highest AMPI (excluding NAs)
  valid_ampi <- !is.na(conflict_with_ampi$AMPI)
  if (sum(valid_ampi) == 0) {
    cat("ERROR: No valid AMPI values found!\n")
    return(invisible(NULL))
  }
  
  max_idx <- which(conflict_with_ampi$AMPI == max(conflict_with_ampi$AMPI, na.rm = TRUE))[1]
  
  cat("=== HIGHEST AMPI OBSERVATION ===\n")
  cat("Row index:", max_idx, "\n")
  if ("ano" %in% names(conflict_with_ampi)) {
    cat("Year:", conflict_with_ampi$ano[max_idx], "\n")
  }
  cat("AMPI:", conflict_with_ampi$AMPI[max_idx], "\n")
  cat("Mean (M_ri):", conflict_with_ampi$AMPI_mean[max_idx], "\n")
  cat("Penalty:", conflict_with_ampi$AMPI_penalty[max_idx], "\n")
  cat("SD (S_ri):", ampi_results$S_ri[max_idx], "\n")
  cat("CV:", conflict_with_ampi$AMPI_cv[max_idx], "\n\n")
  
  # Get the year of max AMPI observation
  max_year <- conflict_with_ampi$ano[max_idx]
  year_data <- conflict_with_ampi %>% filter(ano == max_year)
  max_year_idx <- which(year_data$AMPI == max(year_data$AMPI, na.rm = TRUE))[1]
  
  # Get normalized values for this observation
  year_position <- which(names(ampi_results$R_by_year) == as.character(max_year))
  if (length(year_position) == 0) {
    cat("ERROR: Could not find year in results\n")
    return(invisible(NULL))
  }
  
  R_matrix <- ampi_results$R_by_year[[year_position]]
  
  cat("=== NORMALIZED VALUES (R) FOR HIGHEST AMPI ===\n")
  cat("Indicator normalized values:\n")
  normalized_vals <- R_matrix[max_year_idx, ]
  for (i in seq_along(indicator_vars)) {
    cat("  ", indicator_vars[i], ":", normalized_vals[i], "\n")
  }
  cat("\nRange of normalized values for this observation:", 
      range(R_matrix[max_year_idx, ], na.rm = TRUE), "\n")
  cat("SD of normalized values:", sd(R_matrix[max_year_idx, ], na.rm = TRUE), "\n")
  cat("Number of non-NA indicators:", sum(!is.na(R_matrix[max_year_idx, ])), "\n\n")
  
  # Check if any observation has all values at 130
  cat("=== CHECKING FOR PERFECT SCORES ===\n")
  found_perfect <- FALSE
  for (y_idx in seq_along(ampi_results$R_by_year)) {
    R <- ampi_results$R_by_year[[y_idx]]
    
    # Check if anyone is at max (130) for all indicators (ignoring NAs)
    all_max <- apply(R, 1, function(row) {
      non_na <- row[!is.na(row)]
      if (length(non_na) == 0) return(FALSE)
      all(non_na >= 129.9)
    })
    
    if (any(all_max, na.rm = TRUE)) {
      cat("Year", names(ampi_results$R_by_year)[y_idx], 
          "- Found", sum(all_max, na.rm = TRUE), 
          "observation(s) with all indicators at max\n")
      found_perfect <- TRUE
    }
  }
  if (!found_perfect) {
    cat("No observations found with all indicators at maximum (130)\n")
  }
  cat("\n")
  
  # Check the distribution of maximum values per indicator per year
  cat("=== MAXIMUM NORMALIZED VALUES BY INDICATOR BY YEAR ===\n")
  for (y_idx in seq_along(ampi_results$R_by_year)) {
    R <- ampi_results$R_by_year[[y_idx]]
    year_name <- names(ampi_results$R_by_year)[y_idx]
    
    cat("Year", year_name, ":\n")
    for (j in 1:ncol(R)) {
      max_val <- max(R[, j], na.rm = TRUE)
      min_val <- min(R[, j], na.rm = TRUE)
      cat("  ", indicator_vars[j], "- Min:", round(min_val, 2), 
          "Max:", round(max_val, 2), "\n")
    }
    cat("\n")
  }
  
  # Theoretical maximum AMPI
  cat("=== WHY AMPI < 130 ===\n")
  cat("Theoretical maximum: 130 (requires ALL indicators = 130 AND zero variation)\n")
  cat("Your maximum AMPI:", max(conflict_with_ampi$AMPI, na.rm = TRUE), "\n\n")
  cat("Most likely reasons:\n")
  cat("1. No single observation has ALL indicators at maximum simultaneously\n")
  cat("2. When indicators vary (e.g., 130, 120, 110), penalty reduces AMPI\n")
  cat("3. Missing data (NAs) may reduce the effective maximum\n")
  cat("4. AMPI formula: Mean - (SD × CV) penalizes imbalance\n")
}

# Run the diagnostic
diagnostic_check(conflict_with_ampi, ampi_results, indicator_vars)



##Heatmap
top_municipalities_ampi <- conflict_with_ampi %>%
  group_by(codmpio) %>%
  summarise(avg_ampi = mean(AMPI, na.rm = TRUE)) %>%
  arrange(desc(avg_ampi)) %>%
  head(30) %>%
  pull(codmpio)

##Creating a range for a more meaningful heat map
rng <- range(conflict_with_ampi$AMPI, na.rm = TRUE)
# e.g., set visually useful limits
lims <- c(70, max(85, rng[2]))

q1 <- conflict_with_ampi %>%
  filter(codmpio %in% top_municipalities_ampi) %>%
  ggplot(aes(x = ano, y = as.factor(codmpio), fill = AMPI)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_gradientn(
    colours = c("#2ca02c", "#ff7f0e", "#d62728"),
    values = scales::rescale(c(lims[1], (lims[1]+lims[2])/2, lims[2])),
    limits = lims,
    name = "AMPI\n(70–130)"
  ) +
  labs(title = "Conflict Intensity Heatmap - AMPI",
       subtitle = "Top 30 Most Affected Municipalities",
       x = "Year", y = "Municipality Code") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7))

print(q1)

##Conflict index over time - with Inter quartile range
temporal_trends <- conflict_with_ampi %>%
  group_by(ano) %>%
  summarise(
    mean_ampi = mean(AMPI, na.rm = TRUE),
    median_ampi = median(AMPI, na.rm = TRUE),
    sd_ampi = sd(AMPI, na.rm = TRUE),
    q25 = quantile(AMPI, 0.25, na.rm = TRUE),
    q75 = quantile(AMPI, 0.75, na.rm = TRUE)
  )

q2 <- ggplot(temporal_trends, aes(x = ano)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2, fill = "darkred") +
  geom_line(aes(y = mean_ampi, color = "Mean"), linewidth = 1.2) +
  geom_line(aes(y = median_ampi, color = "Median"), linewidth = 1.2) +
  scale_color_manual(values = c("Mean" = "darkred", "Median" = "steelblue")) +
  labs(
    title = "Conflict Index Evolution Over Time",
    subtitle = "Average AMPI across all municipalities (shaded area = IQR)",
    x = "Year",
    y = "AMPI (Conflict Index)",
    color = "Measure"
  ) +
  theme(legend.position = "bottom")

print(q2)
##Shaded area - high variability across municipalities
##Conflict index is not high but the peaks align with the theory

##MAPPING
dane_shp <- st_read("Conteo_Mpio/MGN_MPIO.shp")
print(names(dane_shp))
dane_shp <- dane_shp %>% rename(codmpio_original = MPIO_CDPMP)

dane_shp$codmpio <- as.numeric(dane_shp$codmpio_original)

conflict_with_ampi$codmpio <- as.numeric(conflict_with_ampi$codmpio)

# Check how many municipalities will match
map_codes <- unique(dane_shp$codmpio)
data_codes <- unique(conflict_with_ampi$codmpio)
matched <- sum(data_codes %in% map_codes)
print(matched)
unmatched <- setdiff(data_codes, map_codes)
print(unmatched)


# Calculate average AMPI per municipality per year
avg_conflict <- conflict_with_ampi %>%
  group_by(codmpio) %>%
  summarise(
    avg_ampi = if (all(is.na(AMPI))) NA else mean(AMPI, na.rm = TRUE),
    max_ampi = if (all(is.na(AMPI))) NA else max(AMPI, na.rm = TRUE),
    min_ampi = if (all(is.na(AMPI))) NA else min(AMPI, na.rm = TRUE),
    years_count = sum(!is.na(AMPI))
  )

# Join with map
map_data_avg <- dane_shp %>%
  left_join(avg_conflict, by = "codmpio")

sum(is.na(map_data_avg$avg_ampi))
## 12 municipalities do not have data

##visualing the AMPI on the map
x1 <- ggplot(map_data_avg) +
  geom_sf(aes(fill = avg_ampi), color = NA) +  
  scale_fill_viridis(option = "plasma", na.value = "grey90",
                     name = "Average\nAMPI") +
  labs(
    title = "Average Conflict Index by Municipality",
    subtitle = "AMPI averaged across all years",
    caption = "Higher values indicate higher conflict levels"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

print(x1)


##AMPI by year across municipalities
map_data_2020 <- dane_shp %>%
  left_join(conflict_with_ampi %>% filter(ano == 2020), by = "codmpio")

ggplot(map_data_2020) +
  geom_sf(aes(fill = AMPI), color = NA) +
  scale_fill_viridis(option = "plasma", na.value = "grey90",
                     name = "AMPI") +
  labs(
    title = "Conflict Index by Municipality - 2020",
    caption = "Higher values indicate higher conflict levels"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

print(map_data_2020)

map_data_2016 <- dane_shp %>%
  left_join(conflict_with_ampi %>% filter(ano == 2016), by = "codmpio")

ggplot(map_data_2016) +
  geom_sf(aes(fill = AMPI), color = NA) +
  scale_fill_viridis(option = "plasma", na.value = "grey90",
                     name = "AMPI") +
  labs(
    title = "Conflict Index by Municipality - 2016",
    caption = "Higher values indicate higher conflict levels"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

print(map_data_2016)

map_data_2002 <- dane_shp %>%
  left_join(conflict_with_ampi %>% filter(ano == 2002), by = "codmpio")

ggplot(map_data_2002) +
  geom_sf(aes(fill = AMPI), color = NA) +
  scale_fill_viridis(option = "plasma", na.value = "grey90",
                     name = "AMPI") +
  labs(
    title = "Conflict Index by Municipality - 2002",
    caption = "Higher values indicate higher conflict levels"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

print(map_data_2002)

map_data_2019 <- dane_shp %>%
  left_join(conflict_with_ampi %>% filter(ano == 2019), by = "codmpio")

ggplot(map_data_2019) +
  geom_sf(aes(fill = AMPI), color = NA) +
  scale_fill_viridis(option = "plasma", na.value = "grey90",
                     name = "AMPI") +
  labs(
    title = "Conflict Index by Municipality - 2019",
    caption = "Higher values indicate higher conflict levels"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

print(map_data_2019)


##Change in conflict over time
##Change in the conflict level over time 
conflict_change <- conflict_with_ampi %>%
  group_by(codmpio) %>%
  arrange(ano) %>%
  summarise(
    first_year_ampi = first(AMPI[!is.na(AMPI)]),
    last_year_ampi = last(AMPI[!is.na(AMPI)]),
    change = last_year_ampi - first_year_ampi,
    percent_change = (change / first_year_ampi) * 100
  )

map_data_change <- dane_shp %>%
  left_join(conflict_change, by = "codmpio")

ggplot(map_data_change) +
  geom_sf(aes(fill = change), color = NA) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, na.value = "grey90",
                       name = "AMPI\nChange") +
  labs(title = "Change in Conflict Levels Over Time",
       subtitle = "Blue = Decreased, Red = Increased")
print(map_data_change)

quantile(conflict_final$MPI, na.rm = TRUE)
quantile(conflict_final$MPI_normalized, na.rm = TRUE)