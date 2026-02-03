# ============================================================================
# EMA Stress Item Multilevel Reliability Analysis
# Between-level and within-level reliability for single-item EMA measure
# Based on Cranford et al. (2006) and Bolger & Laurenceau (2013)
# ============================================================================

# Load required packages
library(psych)
library(tidyverse)
library(lme4)

# ----------------------------------------------------------------------------
# Load and prepare data
# ----------------------------------------------------------------------------

data <- read.csv("ndvi_stress_data.csv")

# Check data structure
cat("=== Data Overview ===\n")
cat("Total observations:", nrow(data), "\n")
cat("Number of participants:", length(unique(data$Participant_ID)), "\n")
cat("Stress variable summary:\n")
print(summary(data$stressed))

# Remove missing values for stress
stress_data <- data %>%
  filter(!is.na(stressed)) %>%
  select(Participant_ID, stressed)

cat("\nAfter removing missing values:\n")
cat("Total observations:", nrow(stress_data), "\n")
cat("Number of participants:", length(unique(stress_data$Participant_ID)), "\n")

# Check observations per participant
obs_per_person <- stress_data %>%
  group_by(Participant_ID) %>%
  summarise(n_obs = n())

cat("\nObservations per participant:\n")
cat("Mean:", round(mean(obs_per_person$n_obs), 2), "\n")
cat("SD:", round(sd(obs_per_person$n_obs), 2), "\n")
cat("Range:", min(obs_per_person$n_obs), "-", max(obs_per_person$n_obs), "\n")

# ----------------------------------------------------------------------------
# Multilevel Reliability Analysis
# Following Cranford et al. (2006) and Shrout & Lane (2012)
# ----------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("=== Multilevel Reliability Analysis for Single-Item EMA ===\n")
cat("============================================================\n\n")

# Fit null model (random intercept only) to decompose variance
null_model <- lmer(stressed ~ 1 + (1 | Participant_ID), data = stress_data)

cat("--- Null Model Summary ---\n")
print(summary(null_model))

# Extract variance components
var_components <- as.data.frame(VarCorr(null_model))
between_var <- var_components$vcov[1]  # Between-person variance (tau^2)
within_var <- var_components$vcov[2]   # Within-person variance (sigma^2)
total_var <- between_var + within_var

cat("\n=== Variance Components ===\n")
cat("Between-person variance (tau^2):", round(between_var, 4), "\n")
cat("Within-person variance (sigma^2):", round(within_var, 4), "\n")
cat("Total variance:", round(total_var, 4), "\n")

# ----------------------------------------------------------------------------
# ICC Calculation
# ----------------------------------------------------------------------------

cat("\n=== Intraclass Correlation Coefficient (ICC) ===\n")

# ICC = tau^2 / (tau^2 + sigma^2)
icc <- between_var / total_var
cat("ICC:", round(icc, 4), "\n")
cat("Interpretation: ", round(icc * 100, 1), "% of variance is between-person\n")
cat("               ", round((1 - icc) * 100, 1), "% of variance is within-person\n")

# ----------------------------------------------------------------------------
# Between-Person Reliability (RkRn or R1R)
# ----------------------------------------------------------------------------

cat("\n=== Between-Person Reliability ===\n")
cat("(Reliability of person means based on k observations)\n\n")

# Average number of observations per person
k_avg <- mean(obs_per_person$n_obs)

# RkRn = (k * ICC) / (1 + (k - 1) * ICC)
# This is the Spearman-Brown formula for the reliability of person means
RkRn <- (k_avg * icc) / (1 + (k_avg - 1) * icc)

cat("Average observations per person (k):", round(k_avg, 2), "\n")
cat("Between-person reliability (RkRn):", round(RkRn, 4), "\n")

# Also calculate for different k values
cat("\nBetween-person reliability at different k values:\n")
k_values <- c(5, 10, 15, 20, 25, 30, 40, 50)
for (k in k_values) {
  r_k <- (k * icc) / (1 + (k - 1) * icc)
  cat("  k =", sprintf("%2d", k), ": RkRn =", round(r_k, 4), "\n")
}

# ----------------------------------------------------------------------------
# Within-Person Reliability (Rc)
# ----------------------------------------------------------------------------

cat("\n=== Within-Person Reliability ===\n")
cat("(Reliability of detecting within-person change)\n\n")

# For single-item measures, within-person reliability can be estimated as:
# Rc = 1 - (sigma^2_error / sigma^2_within)
# However, for single items, we cannot separate true within-person variance
# from error variance without additional assumptions or data

# One approach: Use the proportion of within-person variance
# relative to total variance as a proxy
within_var_prop <- within_var / total_var
cat("Within-person variance proportion:", round(within_var_prop, 4), "\n")

# Alternative: Assume all within-person variance is reliable (upper bound)
cat("\nNote: For single-item measures, within-person reliability cannot be\n")
cat("directly estimated without additional information (e.g., multiple items\n")
cat("or test-retest data). The within-person variance proportion (",
    round(within_var_prop, 4), ")\n")
cat("represents the total within-person variability, which includes both\n")
cat("true fluctuations and measurement error.\n")

# ----------------------------------------------------------------------------
# ICC using psych package (for comparison)
# ----------------------------------------------------------------------------

cat("\n=== ICC Estimates using psych::ICC ===\n\n")

# Prepare data in wide format for psych::ICC
stress_wide <- stress_data %>%
  group_by(Participant_ID) %>%
  mutate(time = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = Participant_ID,
    names_from = time,
    values_from = stressed,
    names_prefix = "T"
  )

# Convert to matrix (remove ID column)
stress_matrix <- as.matrix(stress_wide[, -1])
rownames(stress_matrix) <- stress_wide$Participant_ID

# Run ICC analysis
icc_results <- ICC(stress_matrix, missing = TRUE, lmer = TRUE)
print(icc_results)

cat("\nKey ICC interpretations:\n")
cat("- ICC1: Single measurement reliability (", round(icc_results$results$ICC[1], 4), ")\n")
cat("- ICC(1,k): Average of k measurements reliability (", round(icc_results$results$ICC[2], 4), ")\n")
cat("  where k =", ncol(stress_matrix), "(max time points)\n")

# ----------------------------------------------------------------------------
# Summary Table
# ----------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("=== SUMMARY: Multilevel Reliability Estimates ===\n")
cat("============================================================\n\n")

summary_table <- data.frame(
  Metric = c(
    "Between-person variance (tau^2)",
    "Within-person variance (sigma^2)",
    "ICC",
    paste0("Between-person reliability (RkRn, k=", round(k_avg, 1), ")"),
    "Within-person variance proportion"
  ),
  Estimate = c(
    round(between_var, 4),
    round(within_var, 4),
    round(icc, 4),
    round(RkRn, 4),
    round(within_var_prop, 4)
  )
)

print(summary_table, row.names = FALSE)

# Export results to CSV
write.csv(summary_table, "revision_analysis/ema_stress_reliability.csv", row.names = FALSE)
cat("\nResults saved to: revision_analysis/ema_stress_reliability.csv\n")

cat("\n--- Interpretation Guide ---\n")
cat("ICC: Proportion of variance due to stable between-person differences.\n")
cat("     Higher ICC = more trait-like (stable individual differences).\n")
cat("     Lower ICC = more state-like (within-person fluctuations dominate).\n\n")
cat("RkRn (Between-person reliability): Reliability of distinguishing\n")
cat("     individuals based on their average stress across k observations.\n")
cat("     Values > 0.70 are generally considered acceptable.\n\n")
cat("Within-person variance proportion: The portion of variance reflecting\n")
cat("     momentary fluctuations in stress (both true change and error).\n")

# ----------------------------------------------------------------------------
# References
# ----------------------------------------------------------------------------

cat("\n--- References ---\n")
cat("Cranford, J. A., et al. (2006). A procedure for evaluating sensitivity to\n")
cat("  within-person change. Behavior Research Methods, 38(4), 585-594.\n")
cat("Shrout, P. E., & Lane, S. P. (2012). Psychometrics. In M. R. Mehl & T. S.\n")
cat("  Conner (Eds.), Handbook of research methods for studying daily life.\n")
cat("Bolger, N., & Laurenceau, J. P. (2013). Intensive longitudinal methods.\n")
