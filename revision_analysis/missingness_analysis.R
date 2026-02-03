# ============================================================================
# CORRECTED MISSINGNESS ANALYSIS - Using True 60-Day Enrollment
# ============================================================================
# Uses actual start dates and 60-day protocol for all participants
# This gives the TRUE adherence rates including dropout
# ============================================================================

library(tidyverse)
library(lubridate)
library(lme4)
library(broom.mixed)

# Create output directory
output_dir <- "revision_analysis/missingness_output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Start output capture
output_file <- file.path(output_dir, "missingness_analysis.txt")
sink(output_file, split = TRUE)

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("CORRECTED MISSINGNESS ANALYSIS - TRUE 60-DAY ENROLLMENT\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")
cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ============================================================================
# 1. LOAD DATA AND CREATE COMPLETE 60-DAY SEQUENCES
# ============================================================================

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("DATA PREPARATION\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

# Load actual data
data <- read_csv("ndvi_stress_data.csv") %>%
  mutate(Date = as.Date(Date))

cat("Actual data loaded:", nrow(data), "observations\n")

# Load start dates
start_dates <- read_csv("start_dates.csv") %>%
  mutate(start_date = as.Date(start_date))

cat("Start dates loaded:", nrow(start_dates), "participants\n\n")

# Create race variable
data <- data %>%
  mutate(
    race_white = case_when(
      str_detect(race, "White|European") ~ "White",
      !is.na(race) ~ "Non-White",
      TRUE ~ NA_character_
    ),
    race_white = factor(race_white, levels = c("White", "Non-White"))
  )

# Create COMPLETE 60-day sequence for ALL participants
complete_data <- start_dates %>%
  rowwise() %>%
  mutate(
    end_date = start_date + days(59),  # 60 days total (0-59)
    Date = list(seq(start_date, end_date, by = "day"))
  ) %>%
  unnest(Date) %>%
  group_by(Participant_ID) %>%
  mutate(study_day = as.numeric(Date - min(Date)) + 1) %>%
  ungroup() %>%
  select(Participant_ID, Date, study_day, start_date, end_date)

cat("COMPLETE DATA STRUCTURE CREATED:\n")
cat("Total possible observations (178 × 60):", nrow(complete_data), "\n")
cat("Participants:", n_distinct(complete_data$Participant_ID), "\n")
cat("Study days per participant: 60 (for all)\n\n")

# Join with actual observations
complete_data <- complete_data %>%
  left_join(data, by = c("Participant_ID", "Date"))

# Create missingness indicators
complete_data <- complete_data %>%
  mutate(
    stress_missing = is.na(stressed),
    ndvi_missing = is.na(Daily_TimeWeighted_NDVI),
    both_present = !stress_missing & !ndvi_missing,
    only_stress = !stress_missing & ndvi_missing,
    only_ndvi = stress_missing & !ndvi_missing,
    both_missing = stress_missing & ndvi_missing
  )

# ============================================================================
# 2. TRUE ADHERENCE RATES
# ============================================================================

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("TRUE ADHERENCE RATES (BASED ON 60-DAY PROTOCOL)\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

# Overall patterns
pattern_summary <- complete_data %>%
  summarise(
    total_possible = n(),
    n_both_present = sum(both_present),
    n_only_stress = sum(only_stress),
    n_only_ndvi = sum(only_ndvi),
    n_both_missing = sum(both_missing),
    n_at_least_one = sum(!both_missing),
    pct_both_present = mean(both_present) * 100,
    pct_only_stress = mean(only_stress) * 100,
    pct_only_ndvi = mean(only_ndvi) * 100,
    pct_both_missing = mean(both_missing) * 100,
    pct_at_least_one = mean(!both_missing) * 100
  )

cat("OVERALL OBSERVATION PATTERNS:\n")
cat("------------------------------\n")
cat(sprintf("Total possible observations (178 × 60): %d\n", pattern_summary$total_possible))
cat(sprintf("  Both stress & NDVI present: %d (%.1f%%)\n",
            pattern_summary$n_both_present, pattern_summary$pct_both_present))
cat(sprintf("  Only stress present:        %d (%.1f%%)\n",
            pattern_summary$n_only_stress, pattern_summary$pct_only_stress))
cat(sprintf("  Only NDVI present:          %d (%.1f%%)\n",
            pattern_summary$n_only_ndvi, pattern_summary$pct_only_ndvi))
cat(sprintf("  Neither present (missing):  %d (%.1f%%)\n",
            pattern_summary$n_both_missing, pattern_summary$pct_both_missing))
cat(sprintf("\nObservations with ≥1 variable: %d (%.1f%%)\n\n",
            pattern_summary$n_at_least_one, pattern_summary$pct_at_least_one))

# Separate adherence by measure
adherence_by_measure <- complete_data %>%
  summarise(
    stress_available = sum(!stress_missing),
    ndvi_available = sum(!ndvi_missing),
    stress_adherence = mean(!stress_missing) * 100,
    ndvi_adherence = mean(!ndvi_missing) * 100
  )

cat("ADHERENCE BY MEASURE:\n")
cat("---------------------\n")
cat(sprintf("Stress (EMA): %d observations (%.1f%% of 10,680 possible days)\n",
            adherence_by_measure$stress_available,
            adherence_by_measure$stress_adherence))
cat(sprintf("NDVI (GPS):   %d observations (%.1f%% of 10,680 possible days)\n\n",
            adherence_by_measure$ndvi_available,
            adherence_by_measure$ndvi_adherence))

# Per-participant adherence
participant_adherence <- complete_data %>%
  group_by(Participant_ID) %>%
  summarise(
    Group = first(na.omit(Group)),
    total_days = n(),  # Should be 60 for all
    stress_days = sum(!stress_missing),
    ndvi_days = sum(!ndvi_missing),
    both_days = sum(both_present),
    at_least_one_days = sum(!both_missing),
    stress_adherence = mean(!stress_missing) * 100,
    ndvi_adherence = mean(!ndvi_missing) * 100,
    both_adherence = mean(both_present) * 100,
    at_least_one_adherence = mean(!both_missing) * 100,
    .groups = "drop"
  )

cat("PER-PARTICIPANT ADHERENCE SUMMARY:\n")
cat("----------------------------------\n")
cat(sprintf("Number of participants: %d\n", nrow(participant_adherence)))
cat(sprintf("Days per participant: %d (all participants)\n\n", unique(participant_adherence$total_days)))

cat("Stress (EMA) Adherence:\n")
cat(sprintf("  Mean: %.1f%% (SD = %.1f%%)\n",
            mean(participant_adherence$stress_adherence),
            sd(participant_adherence$stress_adherence)))
cat(sprintf("  Median: %.1f%% (IQR: %.1f%% - %.1f%%)\n",
            median(participant_adherence$stress_adherence),
            quantile(participant_adherence$stress_adherence, 0.25),
            quantile(participant_adherence$stress_adherence, 0.75)))
cat(sprintf("  Range: %.1f%% - %.1f%%\n\n",
            min(participant_adherence$stress_adherence),
            max(participant_adherence$stress_adherence)))

cat("NDVI (GPS) Adherence:\n")
cat(sprintf("  Mean: %.1f%% (SD = %.1f%%)\n",
            mean(participant_adherence$ndvi_adherence),
            sd(participant_adherence$ndvi_adherence)))
cat(sprintf("  Median: %.1f%% (IQR: %.1f%% - %.1f%%)\n",
            median(participant_adherence$ndvi_adherence),
            quantile(participant_adherence$ndvi_adherence, 0.25),
            quantile(participant_adherence$ndvi_adherence, 0.75)))
cat(sprintf("  Range: %.1f%% - %.1f%%\n\n",
            min(participant_adherence$ndvi_adherence),
            max(participant_adherence$ndvi_adherence)))

# Adherence by group
adherence_by_group <- participant_adherence %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    stress_mean = mean(stress_adherence),
    stress_sd = sd(stress_adherence),
    ndvi_mean = mean(ndvi_adherence),
    ndvi_sd = sd(ndvi_adherence),
    .groups = "drop"
  )

cat("ADHERENCE BY GROUP:\n")
cat("-------------------\n")
print(adherence_by_group)
cat("\n")

# Identify dropout (participants who stopped providing data)
dropout_analysis <- complete_data %>%
  group_by(Participant_ID) %>%
  arrange(study_day) %>%
  mutate(
    has_data = !both_missing,
    last_data_day = max(study_day[has_data], na.rm = TRUE)
  ) %>%
  ungroup() %>%
  group_by(Participant_ID) %>%
  summarise(
    last_data_day = first(last_data_day),
    completed_60_days = last_data_day >= 60,
    .groups = "drop"
  )

cat("DROPOUT ANALYSIS:\n")
cat("-----------------\n")
cat(sprintf("Participants with data through day 60: %d (%.1f%%)\n",
            sum(dropout_analysis$completed_60_days),
            mean(dropout_analysis$completed_60_days) * 100))
cat(sprintf("Participants who stopped early (dropout): %d (%.1f%%)\n",
            sum(!dropout_analysis$completed_60_days),
            mean(!dropout_analysis$completed_60_days) * 100))
cat(sprintf("Mean last day of data: %.1f (SD = %.1f)\n\n",
            mean(dropout_analysis$last_data_day),
            sd(dropout_analysis$last_data_day)))

# Save adherence tables
write_csv(participant_adherence,
          file.path(output_dir, "participant_adherence.csv"))
write_csv(pattern_summary,
          file.path(output_dir, "observation_patterns.csv"))
write_csv(dropout_analysis,
          file.path(output_dir, "dropout_analysis.csv"))

# ============================================================================
# 3. VISUALIZATIONS
# ============================================================================

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("CREATING VISUALIZATIONS\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

# Plot 1: Observation patterns
pattern_data <- data.frame(
  Pattern = c("Both Present", "Only Stress", "Only NDVI", "Both Missing"),
  Count = c(pattern_summary$n_both_present, pattern_summary$n_only_stress,
            pattern_summary$n_only_ndvi, pattern_summary$n_both_missing),
  Percentage = c(pattern_summary$pct_both_present, pattern_summary$pct_only_stress,
                 pattern_summary$pct_only_ndvi, pattern_summary$pct_both_missing)
) %>%
  mutate(Pattern = factor(Pattern, levels = c("Both Present", "Only Stress",
                                               "Only NDVI", "Both Missing")))

p1 <- ggplot(pattern_data, aes(x = Pattern, y = Percentage, fill = Pattern)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", Percentage, Count)),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Both Present" = "steelblue",
                               "Only Stress" = "orange",
                               "Only NDVI" = "purple",
                               "Both Missing" = "gray70")) +
  labs(title = "Distribution of Daily Observation Patterns",
       subtitle = "Based on true 60-day enrollment (178 × 60 = 10,680 days)",
       x = "Observation Pattern",
       y = "Percentage of Possible Days (%)") +
  ylim(0, max(pattern_data$Percentage) * 1.15) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "observation_patterns.png"), p1,
       width = 10, height = 6)

# Plot 2: Adherence over time
adherence_by_day <- complete_data %>%
  group_by(study_day) %>%
  summarise(
    stress_adherence = mean(!stress_missing) * 100,
    ndvi_adherence = mean(!ndvi_missing) * 100,
    n_participants = n(),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(stress_adherence, ndvi_adherence),
               names_to = "Measure", values_to = "Adherence")

p2 <- ggplot(adherence_by_day, aes(x = study_day, y = Adherence,
                                    color = Measure, linetype = Measure)) +
  geom_line(linewidth = 1) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = c("stress_adherence" = "steelblue",
                                "ndvi_adherence" = "purple"),
                     labels = c("NDVI (GPS)", "Stress (EMA)")) +
  scale_linetype_manual(values = c("stress_adherence" = "solid",
                                   "ndvi_adherence" = "dashed"),
                        labels = c("NDVI (GPS)", "Stress (EMA)")) +
  labs(title = "Adherence Rates Over 60-Day Study Period",
       subtitle = "Shows dropout over time",
       x = "Study Day",
       y = "Adherence Rate (%)") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "adherence_over_time.png"), p2,
       width = 10, height = 6)

cat("Visualizations saved.\n\n")

# ============================================================================
# 4. MISSINGNESS PREDICTORS WITH RACE
# ============================================================================

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("MISSINGNESS PREDICTORS ANALYSIS\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

# Prepare lagged variables
complete_data <- complete_data %>%
  group_by(Participant_ID) %>%
  arrange(Date) %>%
  mutate(
    stress_lag1 = lag(stressed, 1),
    ndvi_lag1 = lag(Daily_TimeWeighted_NDVI, 1),
    stress_missing_lag1 = lag(stress_missing, 1),
    ndvi_missing_lag1 = lag(ndvi_missing, 1)
  ) %>%
  ungroup()

# Get participant baseline data
participant_baseline <- data %>%
  group_by(Participant_ID) %>%
  summarise(
    Group = first(na.omit(Group)),
    age = first(na.omit(age)),
    sex = first(na.omit(sex)),
    race_white = first(na.omit(race_white)),
    clgry_total_baseline = first(na.omit(clgry_total_baseline)),
    bprs_total_baseline = first(na.omit(bprs_total_baseline)),
    oasis_total_baseline = first(na.omit(oasis_total_baseline)),
    .groups = "drop"
  )

# Join and prepare model data
model_data <- complete_data %>%
  left_join(participant_baseline, by = "Participant_ID") %>%
  mutate(
    Group = coalesce(Group.y, Group.x),
    age = coalesce(age.y, age.x),
    sex = coalesce(sex.y, sex.x),
    race_white = coalesce(race_white.y, race_white.x),
    clgry_total_baseline = coalesce(clgry_total_baseline.y, clgry_total_baseline.x),
    bprs_total_baseline = coalesce(bprs_total_baseline.y, bprs_total_baseline.x),
    oasis_total_baseline = coalesce(oasis_total_baseline.y, oasis_total_baseline.x),
    Group = factor(Group),
    sex = factor(sex)
  ) %>%
  select(-ends_with(".x"), -ends_with(".y")) %>%
  filter(study_day > 1)

# MODEL 1A: Baseline predictors of STRESS missingness
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("MODEL 1A: BASELINE PREDICTORS OF STRESS MISSINGNESS\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

model1a <- glmer(
  stress_missing ~
    Group + age + sex + race_white +
    clgry_total_baseline +
    bprs_total_baseline +
    oasis_total_baseline +
    (1 | Participant_ID),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

model1a_summary <- tidy(model1a, effects = "fixed", conf.int = TRUE,
                       exponentiate = TRUE)

cat("Odds Ratios for Stress Missingness:\n\n")
print(model1a_summary %>% select(term, estimate, conf.low, conf.high, p.value), n = Inf)
cat("\n")

write_csv(model1a_summary,
          file.path(output_dir, "model1a_stress_baseline.csv"))

# MODEL 1B: Time-varying predictors of STRESS missingness
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("MODEL 1B: TIME-VARYING PREDICTORS OF STRESS MISSINGNESS\n")
cat("INCLUDING PREVIOUS NDVI (KEY TEST)\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

model_data_complete <- model_data %>%
  filter(!is.na(stress_lag1), !is.na(ndvi_lag1))

cat(sprintf("Sample size: %d observations from %d participants\n\n",
            nrow(model_data_complete),
            n_distinct(model_data_complete$Participant_ID)))

model1b <- glmer(
  stress_missing ~
    Group + age + sex + race_white +
    clgry_total_baseline +
    bprs_total_baseline +
    oasis_total_baseline +
    study_day +
    stress_lag1 +
    ndvi_lag1 +
    (1 | Participant_ID),
  data = model_data_complete,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

model1b_summary <- tidy(model1b, effects = "fixed", conf.int = TRUE,
                       exponentiate = TRUE)

cat("Odds Ratios for Stress Missingness:\n\n")
print(model1b_summary %>% select(term, estimate, conf.low, conf.high, p.value), n = Inf)
cat("\n")

write_csv(model1b_summary,
          file.path(output_dir, "model1b_stress_timevarying.csv"))

# KEY FINDINGS
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("*** KEY FINDINGS FOR REVIEWER ***\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

ndvi_result <- model1b_summary %>% filter(term == "ndvi_lag1")
stress_result <- model1b_summary %>% filter(term == "stress_lag1")

if (nrow(ndvi_result) > 0) {
  cat("1. DOES GREENSPACE (NDVI) PREDICT STRESS MISSINGNESS?\n\n")
  if (ndvi_result$p.value < 0.05) {
    cat(sprintf("   YES: OR = %.3f (95%% CI: %.3f-%.3f), p = %.4f\n",
                ndvi_result$estimate, ndvi_result$conf.low,
                ndvi_result$conf.high, ndvi_result$p.value))
    cat("   INTERPRETATION: Greenspace exposure predicts missingness (MAR).\n\n")
  } else {
    cat(sprintf("   NO: OR = %.3f (95%% CI: %.3f-%.3f), p = %.4f\n",
                ndvi_result$estimate, ndvi_result$conf.low,
                ndvi_result$conf.high, ndvi_result$p.value))
    cat("   INTERPRETATION: Greenspace NOT related to missingness.\n\n")
  }
}

if (nrow(stress_result) > 0) {
  cat("2. DOES PREVIOUS STRESS PREDICT STRESS MISSINGNESS?\n\n")
  if (stress_result$p.value < 0.05) {
    cat(sprintf("   YES: OR = %.3f (95%% CI: %.3f-%.3f), p = %.4f\n",
                stress_result$estimate, stress_result$conf.low,
                stress_result$conf.high, stress_result$p.value))
    cat("   INTERPRETATION: Concern about MNAR mechanism.\n\n")
  } else {
    cat(sprintf("   NO: OR = %.3f (95%% CI: %.3f-%.3f), p = %.4f\n",
                stress_result$estimate, stress_result$conf.low,
                stress_result$conf.high, stress_result$p.value))
    cat("   INTERPRETATION: Less concern about MNAR.\n\n")
  }
}

# MODEL 2A: Baseline predictors of NDVI missingness
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("MODEL 2A: BASELINE PREDICTORS OF NDVI MISSINGNESS\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

model2a <- glmer(
  ndvi_missing ~
    Group + age + sex + race_white +
    clgry_total_baseline +
    bprs_total_baseline +
    oasis_total_baseline +
    (1 | Participant_ID),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

model2a_summary <- tidy(model2a, effects = "fixed", conf.int = TRUE,
                       exponentiate = TRUE)

cat("Odds Ratios for NDVI Missingness:\n\n")
print(model2a_summary %>% select(term, estimate, conf.low, conf.high, p.value), n = Inf)
cat("\n")

write_csv(model2a_summary,
          file.path(output_dir, "model2a_ndvi_baseline.csv"))

# MODEL 2B: Time-varying predictors of NDVI missingness
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("MODEL 2B: TIME-VARYING PREDICTORS OF NDVI MISSINGNESS\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

model2b <- glmer(
  ndvi_missing ~
    Group + age + sex + race_white +
    clgry_total_baseline +
    bprs_total_baseline +
    oasis_total_baseline +
    study_day +
    stress_lag1 +
    ndvi_lag1 +
    (1 | Participant_ID),
  data = model_data_complete,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

model2b_summary <- tidy(model2b, effects = "fixed", conf.int = TRUE,
                       exponentiate = TRUE)

cat("Odds Ratios for NDVI Missingness:\n\n")
print(model2b_summary %>% select(term, estimate, conf.low, conf.high, p.value), n = Inf)
cat("\n")

write_csv(model2b_summary,
          file.path(output_dir, "model2b_ndvi_timevarying.csv"))

# ============================================================================
# 5. CREATE SUMMARY TABLE
# ============================================================================

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("CREATING SUMMARY TABLE\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

# Function to format results
format_result <- function(summary_df, term_name) {
  result <- summary_df %>% filter(term == term_name)
  if (nrow(result) == 0) return("--")

  or <- sprintf("%.2f", result$estimate)
  ci_low <- sprintf("%.2f", result$conf.low)
  ci_high <- sprintf("%.2f", result$conf.high)
  p_val <- result$p.value

  if (p_val < 0.001) {
    p_str <- "<.001***"
  } else if (p_val < 0.01) {
    p_str <- sprintf("%.3f**", p_val)
  } else if (p_val < 0.05) {
    p_str <- sprintf("%.3f*", p_val)
  } else if (p_val < 0.10) {
    p_str <- sprintf("%.3f†", p_val)
  } else {
    p_str <- sprintf("%.3f", p_val)
  }

  return(sprintf("OR=%s [%s,%s], p=%s", or, ci_low, ci_high, p_str))
}

# Create comprehensive summary table
summary_table <- data.frame(
  Predictor = c(
    "Baseline Demographics:",
    "  Group (HC vs CHR)",
    "  Age",
    "  Sex (Male vs Female)",
    "  Race (Non-White vs White)",
    "Baseline Clinical:",
    "  Calgary Depression",
    "  BPRS (Psychosis)",
    "  OASIS (Anxiety)",
    "Time-Varying:",
    "  Study Day",
    "  Previous Stress",
    "  Previous NDVI (Greenspace)",
    "  Previous Missingness"
  ),
  Stress_Missingness_Baseline = c(
    "",
    format_result(model1a_summary, "GroupHC"),
    format_result(model1a_summary, "age"),
    format_result(model1a_summary, "sexM"),
    format_result(model1a_summary, "race_whiteNon-White"),
    "",
    format_result(model1a_summary, "clgry_total_baseline"),
    format_result(model1a_summary, "bprs_total_baseline"),
    format_result(model1a_summary, "oasis_total_baseline"),
    "",
    "--",
    "--",
    "--",
    "--"
  ),
  Stress_Missingness_TimeVarying = c(
    "",
    format_result(model1b_summary, "GroupHC"),
    format_result(model1b_summary, "age"),
    format_result(model1b_summary, "sexM"),
    format_result(model1b_summary, "race_whiteNon-White"),
    "",
    "--",
    "--",
    "--",
    "",
    format_result(model1b_summary, "study_day"),
    format_result(model1b_summary, "stress_lag1"),
    format_result(model1b_summary, "ndvi_lag1"),
    format_result(model1b_summary, "stress_missing_lag1TRUE")
  ),
  NDVI_Missingness_Baseline = c(
    "",
    format_result(model2a_summary, "GroupHC"),
    format_result(model2a_summary, "age"),
    format_result(model2a_summary, "sexM"),
    format_result(model2a_summary, "race_whiteNon-White"),
    "",
    format_result(model2a_summary, "clgry_total_baseline"),
    format_result(model2a_summary, "bprs_total_baseline"),
    format_result(model2a_summary, "oasis_total_baseline"),
    "",
    "--",
    "--",
    "--",
    "--"
  ),
  NDVI_Missingness_TimeVarying = c(
    "",
    format_result(model2b_summary, "GroupHC"),
    format_result(model2b_summary, "age"),
    format_result(model2b_summary, "sexM"),
    format_result(model2b_summary, "race_whiteNon-White"),
    "",
    "--",
    "--",
    "--",
    "",
    format_result(model2b_summary, "study_day"),
    format_result(model2b_summary, "stress_lag1"),
    format_result(model2b_summary, "ndvi_lag1"),
    format_result(model2b_summary, "ndvi_missing_lag1TRUE")
  )
)

write_csv(summary_table,
          file.path(output_dir, "missingness_predictors_summary.csv"))

cat("Summary table saved.\n\n")

# ============================================================================
# 6. FINAL SUMMARY
# ============================================================================

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("ANALYSIS COMPLETE\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

cat("Files generated:\n")
cat("  - participant_adherence.csv\n")
cat("  - observation_patterns.csv\n")
cat("  - dropout_analysis.csv\n")
cat("  - model1a_stress_baseline.csv\n")
cat("  - model1b_stress_timevarying.csv\n")
cat("  - model2a_ndvi_baseline.csv\n")
cat("  - model2b_ndvi_timevarying.csv\n")
cat("  - missingness_predictors_summary.csv\n")
cat("  - observation_patterns.png\n")
cat("  - adherence_over_time.png\n\n")

cat("This analysis uses the TRUE 60-day enrollment period for all participants,\n")
cat("providing accurate adherence rates and missingness patterns including dropout.\n\n")

sink()

cat("Analysis complete! Output saved to:", output_file, "\n")
