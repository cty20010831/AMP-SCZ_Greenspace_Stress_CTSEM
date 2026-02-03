# ============================================================================
# EMA Stress Item Convergent Validity Analysis
# Correlating aggregated daily EMA stress with monthly PSS assessments
# Based on Þórarinsdóttir et al. (2019) and Siepe et al. (2025)
# ============================================================================

# Load required packages
library(tidyverse)
library(lubridate)
library(psych)

# ----------------------------------------------------------------------------
# Load and prepare data
# ----------------------------------------------------------------------------

# Load EMA data
ema_data <- read.csv("ndvi_stress_data.csv")

# Load clinical monthly data
clinical_data <- read.csv("clinical_monthly.csv")

# Load start dates for proper time window alignment
start_dates <- read.csv("start_dates.csv")
start_dates$start_date <- as.Date(start_dates$start_date)

cat("=== Data Overview ===\n")
cat("EMA data: ", nrow(ema_data), " observations from ",
    length(unique(ema_data$Participant_ID)), " participants\n")
cat("Clinical data: ", nrow(clinical_data), " participants\n")
cat("Start dates: ", nrow(start_dates), " participants\n\n")

# Convert Date column
ema_data$Date <- as.Date(ema_data$Date)

# ----------------------------------------------------------------------------
# Define time windows based on PSS reference period
# ----------------------------------------------------------------------------
# PSS asks about stress "in the past month"
# - Baseline PSS: assessed at study start, reflects month BEFORE study
# - Month 1 PSS: assessed ~30 days in, reflects days 0-30 of EMA
# - Month 2 PSS: assessed ~60 days in, reflects days 30-60 of EMA

cat("=== Aggregating EMA Stress by Assessment Period ===\n\n")

# Merge EMA data with start dates
ema_data <- ema_data %>%
  left_join(start_dates, by = "Participant_ID") %>%
  filter(!is.na(start_date)) %>%
  mutate(
    days_since_start = as.numeric(Date - start_date)
  )

# Define assessment periods that match PSS reference windows
# Month 1 PSS (assessed ~day 30) reflects stress from days 0-30
# Month 2 PSS (assessed ~day 60) reflects stress from days 30-60
ema_data <- ema_data %>%
  mutate(
    assessment_period = case_when(
      days_since_start >= 0 & days_since_start < 30 ~ "m1",   # First month of EMA -> Month 1 PSS
      days_since_start >= 30 & days_since_start < 60 ~ "m2",  # Second month of EMA -> Month 2 PSS
      TRUE ~ NA_character_
    )
  )

# Calculate participant-level aggregates by assessment period
ema_monthly <- ema_data %>%
  filter(!is.na(stressed) & !is.na(assessment_period)) %>%
  group_by(Participant_ID, assessment_period) %>%
  summarise(
    stress_mean = mean(stressed, na.rm = TRUE),
    stress_sd = sd(stressed, na.rm = TRUE),
    n_days = n(),
    first_day = min(days_since_start),
    last_day = max(days_since_start),
    .groups = "drop"
  )

cat("EMA aggregates by assessment period:\n")
print(ema_monthly %>%
        group_by(assessment_period) %>%
        summarise(
          n_participants = n(),
          mean_stress = round(mean(stress_mean), 2),
          sd_stress = round(sd(stress_mean), 2),
          mean_days = round(mean(n_days), 1),
          min_days = min(n_days),
          max_days = max(n_days)
        ))

# ----------------------------------------------------------------------------
# Prepare PSS data
# ----------------------------------------------------------------------------

# Reshape PSS data to long format
pss_long <- clinical_data %>%
  select(Participant_ID, Group, pss_total_baseline, pss_total_m1, pss_total_m2) %>%
  pivot_longer(
    cols = starts_with("pss_total"),
    names_to = "timepoint",
    names_prefix = "pss_total_",
    values_to = "pss"
  )

# ----------------------------------------------------------------------------
# Merge EMA and PSS data
# ----------------------------------------------------------------------------

cat("\n=== Merging EMA and PSS Data ===\n\n")

# Merge by participant and matching time period
# Note: We match m1 PSS with m1 EMA (first month) and m2 PSS with m2 EMA (second month)
merged_data <- ema_monthly %>%
  inner_join(
    pss_long %>% filter(timepoint %in% c("m1", "m2")),
    by = c("Participant_ID", "assessment_period" = "timepoint")
  ) %>%
  filter(!is.na(pss))

cat("Merged data: ", nrow(merged_data), " participant-period observations\n")
cat("Participants with merged data: ", length(unique(merged_data$Participant_ID)), "\n\n")

cat("Observations by period:\n")
print(merged_data %>%
        group_by(assessment_period) %>%
        summarise(
          n = n(),
          mean_ema_stress = round(mean(stress_mean), 2),
          mean_pss = round(mean(pss), 2)
        ))

# ============================================================================
# CONVERGENT VALIDITY ANALYSIS
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("=== CONVERGENT VALIDITY: EMA Stress vs PSS ===\n")
cat("============================================================\n\n")

cat("Following Siepe et al. (2025): 'Here we investigate how validated\n")
cat("retrospective reports over a longer period relate to aggregates\n")
cat("of daily reports.'\n\n")

cat("PSS Reference Period: 'In the last month...'\n")
cat("- Month 1 PSS (assessed ~day 30) <-> EMA aggregated from days 0-30\n")
cat("- Month 2 PSS (assessed ~day 60) <-> EMA aggregated from days 30-60\n\n")

# --- Correlation by Assessment Period ---
cat("--- Correlation by Assessment Period ---\n\n")

results_by_period <- data.frame()

for (period in c("m1", "m2")) {
  subset_data <- merged_data %>% filter(assessment_period == period)
  n_obs <- nrow(subset_data)

  if (n_obs >= 10) {
    cor_result <- cor.test(subset_data$stress_mean, subset_data$pss,
                           use = "pairwise.complete.obs")

    period_label <- ifelse(period == "m1", "Month 1", "Month 2")

    results_by_period <- rbind(results_by_period, data.frame(
      Period = period_label,
      n = n_obs,
      mean_EMA_days = round(mean(subset_data$n_days), 1),
      r = cor_result$estimate,
      CI_lower = cor_result$conf.int[1],
      CI_upper = cor_result$conf.int[2],
      t = cor_result$statistic,
      df = cor_result$parameter,
      p_value = cor_result$p.value
    ))

    cat(period_label, ":\n")
    cat("  n =", n_obs, "participants\n")
    cat("  Mean EMA days =", round(mean(subset_data$n_days), 1), "\n")
    cat("  Pearson r =", round(cor_result$estimate, 4), "\n")
    cat("  95% CI: [", round(cor_result$conf.int[1], 4), ", ",
        round(cor_result$conf.int[2], 4), "]\n")
    cat("  t(", cor_result$parameter, ") =", round(cor_result$statistic, 3), "\n")
    cat("  p-value =", format.pval(cor_result$p.value, digits = 4), "\n\n")
  }
}

# --- Overall Correlation (Pooled across periods) ---
cat("--- Overall Correlation (Pooled M1 + M2) ---\n")
cor_overall <- cor.test(merged_data$stress_mean, merged_data$pss,
                         use = "pairwise.complete.obs")

cat("n =", nrow(merged_data), "observations\n")
cat("Pearson r =", round(cor_overall$estimate, 4), "\n")
cat("95% CI: [", round(cor_overall$conf.int[1], 4), ", ",
    round(cor_overall$conf.int[2], 4), "]\n")
cat("t(", cor_overall$parameter, ") =", round(cor_overall$statistic, 3), "\n")
cat("p-value =", format.pval(cor_overall$p.value, digits = 4), "\n")

# Add overall to results
results_by_period <- rbind(results_by_period, data.frame(
  Period = "Overall (M1+M2)",
  n = nrow(merged_data),
  mean_EMA_days = round(mean(merged_data$n_days), 1),
  r = cor_overall$estimate,
  CI_lower = cor_overall$conf.int[1],
  CI_upper = cor_overall$conf.int[2],
  t = cor_overall$statistic,
  df = cor_overall$parameter,
  p_value = cor_overall$p.value
))

# --- Between-Person Correlation (Person Means) ---
cat("\n--- Between-Person Correlation (Person Means) ---\n")
cat("Aggregating across both assessment periods per participant\n\n")

person_means <- merged_data %>%
  group_by(Participant_ID) %>%
  summarise(
    stress_person_mean = mean(stress_mean, na.rm = TRUE),
    pss_person_mean = mean(pss, na.rm = TRUE),
    n_periods = n(),
    .groups = "drop"
  )

cor_between <- cor.test(person_means$stress_person_mean, person_means$pss_person_mean,
                         use = "pairwise.complete.obs")

cat("n =", nrow(person_means), "participants\n")
cat("Pearson r =", round(cor_between$estimate, 4), "\n")
cat("95% CI: [", round(cor_between$conf.int[1], 4), ", ",
    round(cor_between$conf.int[2], 4), "]\n")
cat("t(", cor_between$parameter, ") =", round(cor_between$statistic, 3), "\n")
cat("p-value =", format.pval(cor_between$p.value, digits = 4), "\n")

# Add between-person to results
results_by_period <- rbind(results_by_period, data.frame(
  Period = "Between-Person",
  n = nrow(person_means),
  mean_EMA_days = NA,
  r = cor_between$estimate,
  CI_lower = cor_between$conf.int[1],
  CI_upper = cor_between$conf.int[2],
  t = cor_between$statistic,
  df = cor_between$parameter,
  p_value = cor_between$p.value
))

# ============================================================================
# BASELINE PSS vs OVERALL EMA STRESS (Following Siepe et al., 2025)
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("=== BASELINE PSS vs OVERALL EMA STRESS ===\n")
cat("============================================================\n\n")

cat("Following Siepe et al. (2025, Figure 10): 'Here we investigate how\n")
cat("validated retrospective reports over a longer period relate to\n")
cat("aggregates of daily reports.'\n\n")

cat("This analysis correlates baseline PSS (reflecting stress in the month\n")
cat("BEFORE the study) with overall aggregated EMA stress (across the entire\n")
cat("study period). This tests whether pre-study trait-like stress levels\n")
cat("predict daily stress experiences during the study.\n\n")

# Calculate overall EMA stress mean per participant (entire study period)
overall_ema_stress <- ema_data %>%
  filter(!is.na(stressed)) %>%
  group_by(Participant_ID) %>%
  summarise(
    ema_stress_overall = mean(stressed, na.rm = TRUE),
    ema_stress_sd = sd(stressed, na.rm = TRUE),
    n_ema_days = n(),
    .groups = "drop"
  )

# Get baseline PSS
baseline_pss <- clinical_data %>%
  select(Participant_ID, Group, pss_total_baseline) %>%
  filter(!is.na(pss_total_baseline))

# Merge
baseline_validity <- overall_ema_stress %>%
  inner_join(baseline_pss, by = "Participant_ID")

cat("Sample size: n =", nrow(baseline_validity), "participants\n")
cat("Mean EMA observations per participant:", round(mean(baseline_validity$n_ema_days), 1), "\n\n")

# Correlation: Baseline PSS vs Overall EMA Stress
cor_baseline <- cor.test(baseline_validity$ema_stress_overall,
                          baseline_validity$pss_total_baseline,
                          use = "pairwise.complete.obs")

cat("--- Baseline PSS vs Overall EMA Stress ---\n")
cat("Pearson r =", round(cor_baseline$estimate, 4), "\n")
cat("95% CI: [", round(cor_baseline$conf.int[1], 4), ", ",
    round(cor_baseline$conf.int[2], 4), "]\n")
cat("t(", cor_baseline$parameter, ") =", round(cor_baseline$statistic, 3), "\n")
cat("p-value =", format.pval(cor_baseline$p.value, digits = 4), "\n\n")

# Add to results
results_by_period <- rbind(results_by_period, data.frame(
  Period = "Baseline PSS vs Overall EMA",
  n = nrow(baseline_validity),
  mean_EMA_days = round(mean(baseline_validity$n_ema_days), 1),
  r = cor_baseline$estimate,
  CI_lower = cor_baseline$conf.int[1],
  CI_upper = cor_baseline$conf.int[2],
  t = cor_baseline$statistic,
  df = cor_baseline$parameter,
  p_value = cor_baseline$p.value
))

cat("Note: Siepe et al. (2025) found r = 0.47 between their stressed EMA\n")
cat("item and baseline PSS-10 in the WARN-D study (student sample).\n")

# ============================================================================
# VISUALIZATION
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("=== Creating Visualizations ===\n")
cat("============================================================\n\n")

# 1. Scatter plot: EMA stress vs PSS by assessment period
png("revision_analysis/validity_ema_pss_by_period.png", width = 1000, height = 500, res = 120)

p1 <- ggplot(merged_data, aes(x = pss, y = stress_mean)) +
  geom_point(aes(color = assessment_period), alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  facet_wrap(~assessment_period, labeller = labeller(
    assessment_period = c("m1" = "Month 1", "m2" = "Month 2")
  )) +
  geom_smooth(method = "lm", se = TRUE, aes(color = assessment_period)) +
  scale_color_manual(values = c("m1" = "#2E86AB", "m2" = "#A23B72"),
                     labels = c("m1" = "Month 1", "m2" = "Month 2")) +
  labs(
    title = "Convergent Validity: Aggregated EMA Stress vs Monthly PSS",
    subtitle = "EMA daily stress aggregated to match PSS 'past month' reference period",
    x = "PSS Score (past month)",
    y = "Mean EMA Stress (daily reports)",
    color = "Period"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11)
  )

print(p1)
dev.off()
cat("Saved: revision_analysis/validity_ema_pss_by_period.png\n")

# 2. Combined scatter plot with overall correlation
png("revision_analysis/validity_ema_pss_overall.png", width = 800, height = 600, res = 120)

p2 <- ggplot(merged_data, aes(x = pss, y = stress_mean)) +
  geom_point(aes(color = assessment_period, shape = assessment_period),
             alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = c("m1" = "#2E86AB", "m2" = "#A23B72"),
                     labels = c("m1" = "Month 1", "m2" = "Month 2")) +
  scale_shape_manual(values = c("m1" = 16, "m2" = 17),
                     labels = c("m1" = "Month 1", "m2" = "Month 2")) +
  labs(
    title = "Convergent Validity: EMA Stress vs PSS",
    subtitle = paste0("Overall: r = ", round(cor_overall$estimate, 3),
                      ", 95% CI [", round(cor_overall$conf.int[1], 3), ", ",
                      round(cor_overall$conf.int[2], 3), "]",
                      ", p ", ifelse(cor_overall$p.value < 0.001, "< 0.001",
                                     paste0("= ", round(cor_overall$p.value, 3)))),
    x = "PSS Score (Perceived Stress Scale)",
    y = "Mean Daily EMA Stress",
    color = "Assessment Period",
    shape = "Assessment Period"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(p2)
dev.off()
cat("Saved: revision_analysis/validity_ema_pss_overall.png\n")

# 3. Scatter plot: Baseline PSS vs Overall EMA Stress
png("revision_analysis/validity_baseline_pss_overall_ema.png", width = 800, height = 600, res = 120)

p3_baseline <- ggplot(baseline_validity, aes(x = pss_total_baseline, y = ema_stress_overall)) +
  geom_point(aes(color = Group), alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = c("CHR" = "#E74C3C", "HC" = "#3498DB")) +
  labs(
    title = "Baseline PSS vs Overall Aggregated EMA Stress",
    subtitle = paste0("r = ", round(cor_baseline$estimate, 3),
                      ", 95% CI [", round(cor_baseline$conf.int[1], 3), ", ",
                      round(cor_baseline$conf.int[2], 3), "]",
                      ", p ", ifelse(cor_baseline$p.value < 0.001, "< 0.001",
                                     paste0("= ", round(cor_baseline$p.value, 3)))),
    x = "Baseline PSS (month before study)",
    y = "Overall Mean EMA Stress (entire study period)",
    color = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(p3_baseline)
dev.off()
cat("Saved: revision_analysis/validity_baseline_pss_overall_ema.png\n")

# 4. Forest plot of correlations by period
png("revision_analysis/validity_correlation_forest.png", width = 900, height = 550, res = 120)

forest_data <- results_by_period %>%
  mutate(Period = factor(Period, levels = rev(c("Month 1", "Month 2",
                                                  "Overall (M1+M2)", "Between-Person",
                                                  "Baseline PSS vs Overall EMA"))))

p3 <- ggplot(forest_data, aes(x = r, y = Period)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, linewidth = 0.8) +
  geom_point(size = 3, color = "#2E86AB") +
  geom_text(aes(label = paste0("r = ", round(r, 2), " [", round(CI_lower, 2),
                                ", ", round(CI_upper, 2), "]")),
            hjust = -0.1, vjust = -0.8, size = 3.5) +
  scale_x_continuous(limits = c(-0.1, max(forest_data$CI_upper) + 0.3)) +
  labs(
    title = "Convergent Validity: EMA Stress - PSS Correlations",
    subtitle = "Pearson correlations with 95% confidence intervals",
    x = "Correlation (r)",
    y = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 11)
  )

print(p3)
dev.off()
cat("Saved: revision_analysis/validity_correlation_forest.png\n")

# ============================================================================
# SUMMARY TABLE
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("=== SUMMARY: Convergent Validity Results ===\n")
cat("============================================================\n\n")

# Format results table
summary_table <- results_by_period %>%
  mutate(
    across(c(r, CI_lower, CI_upper), ~round(., 3)),
    t = round(t, 2),
    p_value = format.pval(p_value, digits = 3)
  ) %>%
  select(Period, n, mean_EMA_days, r, CI_lower, CI_upper, p_value)

colnames(summary_table) <- c("Period", "n", "Mean EMA Days", "r", "CI Lower",
                              "CI Upper", "p-value")

print(summary_table)

# Export results
write.csv(results_by_period, "revision_analysis/ema_stress_validity.csv", row.names = FALSE)
cat("\nResults saved to: revision_analysis/ema_stress_validity.csv\n")

# ============================================================================
# INTERPRETATION
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("=== INTERPRETATION ===\n")
cat("============================================================\n\n")

cat("CONVERGENT VALIDITY INTERPRETATION:\n\n")

cat("The single-item daily EMA stress measure was validated against the\n")
cat("Perceived Stress Scale (PSS), a well-established 10-item measure\n")
cat("that assesses stress 'in the past month.'\n\n")

cat("TWO APPROACHES WERE USED:\n\n")

cat("1. CONCURRENT VALIDITY (Month 1 & Month 2):\n")
cat("   Following Siepe et al. (2025), we aggregated daily EMA stress reports\n")
cat("   to match the PSS reference period:\n")
cat("   - Month 1 PSS was correlated with EMA from days 0-30\n")
cat("   - Month 2 PSS was correlated with EMA from days 30-60\n")
cat("   This tests whether aggregated daily reports correspond to\n")
cat("   retrospective monthly assessments of the SAME time period.\n\n")

cat("2. BASELINE PSS vs OVERALL EMA (Siepe et al., 2025 Figure 10):\n")
cat("   Baseline PSS reflects stress in the month BEFORE the study started.\n")
cat("   Overall EMA stress is aggregated across the ENTIRE study period.\n")
cat("   This tests whether pre-study trait-like stress levels predict\n")
cat("   daily stress experiences during the study.\n")
cat("   Reference: Siepe et al. found r = 0.47 between EMA stressed and\n")
cat("   baseline PSS-10 in their student sample.\n\n")

cat("EXPECTED CORRELATION MAGNITUDE:\n")
cat("- Moderate correlations (r ~ 0.30-0.50) are expected and appropriate\n")
cat("- Lower correlations than multi-item scales are normal because:\n")
cat("  1. Different measurement modalities (EMA vs retrospective)\n")
cat("  2. Different temporal granularity (daily vs monthly recall)\n")
cat("  3. Retrospective bias in PSS vs real-time EMA capture\n")
cat("- Þórarinsdóttir et al. (2019) found significant but modest\n")
cat("  associations between single-item smartphone stress and PSS\n\n")

cat("VALIDITY EVIDENCE:\n")
cat("- Significant positive correlations support convergent validity\n")
cat("- The EMA stress item captures construct-relevant variance that\n")
cat("  aligns with the established PSS measure\n")
cat("- Consistency across concurrent (M1, M2) and baseline-overall\n")
cat("  analyses strengthens validity evidence\n\n")

cat("--- References ---\n")
cat("Siepe, B., et al. (2025). Understanding EMA item performance.\n")
cat("  Advances in Methods and Practices in Psychological Science.\n")
cat("Þórarinsdóttir, H., et al. (2019). Smartphone-based stress validation.\n")
cat("Cohen, S., et al. (1983). Perceived Stress Scale.\n")