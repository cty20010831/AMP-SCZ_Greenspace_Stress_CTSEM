Supplementary Materials: Psychometric Properties of the Single-Item EMA Stress Measure

S1. Reliability and Validity of the Daily Stress Item

S1.1 Overview

The daily stress item asked participants to rate "How stressed are you right now?" on a 7-point Likert scale (1 = not at all to 7 = very much). Single-item measures offer practical advantages for EMA research, including reduced participant burden and lower data processing costs (Diamantopoulos et al., 2012). However, given that single-item EMA measures cannot be evaluated using traditional internal consistency approaches (e.g., Cronbach's alpha), we employed multilevel reliability methods following established guidelines for intensive longitudinal data (Cranford et al., 2006; Shrout & Lane, 2012; Trull & Ebner-Priemer, 2024; Castro-Alvarez et al., 2024). Convergent validity was assessed by comparing EMA stress with validated retrospective measures (Siepe et al., 2025; Faurholt-Jepsen et al., 2019).

S1.2 Methods

S1.2.1 Reliability Analysis

Reliability estimates are not commonly reported in applied EMA research despite being a core psychometric feature (Trull & Ebner-Priemer, 2024). We assessed the reliability of the single-item stress measure using multilevel variance decomposition, which distinguishes errors in person-level measures from errors in the measurement of different moments within a person (Cranford et al., 2006; Bolger & Laurenceau, 2013; Castro-Alvarez et al., 2024). A null multilevel model with random intercepts was fitted to partition total variance into between-person and within-person components:

$$\text{Stress}_{ti} = \gamma_{00} + u_{0i} + e_{ti}$$

where $\gamma_{00}$ is the grand mean, $u_{0i}$ is the random intercept for person $i$ (with variance $\tau^2$), and $e_{ti}$ is the residual for person $i$ at time $t$ (with variance $\sigma^2$).

From this model, we calculated:

Intraclass Correlation Coefficient (ICC): The proportion of total variance attributable to stable between-person differences:

$$\text{ICC} = \frac{\tau^2}{\tau^2 + \sigma^2}$$

The ICC indicates the degree to which the stress item captures trait-like individual differences versus state-like momentary fluctuations (Trull & Ebner-Priemer, 2024). Between-person variability is assessed via ICCs, which help distinguish moment-to-moment fluctuation from trait-like differences (Brose et al., 2020).

Between-Person Reliability ($R_{kRn}$): For individual differences in means, ICCs with adjustment using the Spearman-Brown prediction formula are an indicator of reliability (Trull & Ebner-Priemer, 2024). The reliability of person-mean aggregates based on $k$ observations was calculated as:

$$R_{kRn} = \frac{k \times \text{ICC}}{1 + (k - 1) \times \text{ICC}}$$

This metric indicates how reliably the measure distinguishes between individuals when averaging across multiple daily assessments.

S1.2.2 Convergent Validity Analysis

Convergent validity was assessed by correlating aggregated EMA stress with the Perceived Stress Scale (PSS; Cohen et al., 1983), a well-validated 10-item retrospective measure that assesses stress "in the past month." This approach follows recommendations for comparing dynamic momentary items with slower-moving or trait-like measures (Siepe et al., 2025; Faurholt-Jepsen et al., 2019). In traditional validation work, such analyses establish concurrent validity (how well a daily item overlaps with less frequent assessments of the same construct) and predictive validity (how well baseline measures predict daily experiences; Siepe et al., 2025; Diamantopoulos et al., 2012).

We employed two complementary approaches:

Approach 1: Concurrent Validity (Time-Matched Aggregation). Daily EMA stress reports were aggregated to match the PSS reference period. Using each participant's study start date, we defined assessment windows as follows: Month 1 EMA from days 0-30 was correlated with Month 1 PSS (assessed around day 30), and Month 2 EMA from days 30-60 was correlated with Month 2 PSS (assessed around day 60). This approach tests whether aggregated daily reports correspond to retrospective monthly assessments of the same time period.

Approach 2: Baseline PSS vs. Overall EMA Stress. Following Siepe et al. (2025, Figure 10), we correlated baseline PSS (reflecting stress in the month before the study) with overall EMA stress aggregated across the entire study period. This approach tests whether pre-study trait-like stress levels predict daily stress experiences during the study, similar to the validation approach used by Faurholt-Jepsen et al. (2019) comparing daily smartphone-based stress ratings with the PSS.

Pearson correlations with 95% confidence intervals were calculated for all analyses.

S1.3 Results

S1.3.1 Reliability

Variance decomposition revealed that 37.1% of the total variance in daily stress ratings was attributable to stable between-person differences (ICC = 0.37), while 62.9% reflected within-person fluctuations. This pattern indicates that the stress item captures both trait-like individual differences and state-like momentary variations, consistent with stress being conceptualized as having both stable and dynamic components.

The within-person variance proportion (63%) is comparable to findings from coordinated analyses of affect in daily life. Brose et al. (2020) reported within-person variance ranging from 45-66% for negative affect and 25-74% for positive affect across seven EMA and diary studies. The substantial within-person variance observed in our data is desirable for intensive longitudinal analyses, as it indicates sufficient variability to detect within-person associations.

The between-person reliability ($R_{kRn}$) based on the average number of observations per participant (k = 26) was 0.94, exceeding the conventional threshold of 0.70 for acceptable reliability. This indicates that person-mean stress scores reliably distinguish between individuals when aggregated across the study period. This high reliability is consistent with Brose et al. (2020), who reported between-person reliability ranging from 0.96-1.00 across multiple EMA studies, noting that "reliability was adequate to high at all levels of analysis despite different items and designs."

Table S1. Multilevel Reliability Estimates for the Daily Stress Item

Between-person variance (tau-squared): 1.11
Within-person variance (sigma-squared): 1.89
Intraclass Correlation (ICC): 0.37
Between-person reliability (RkRn, k = 26): 0.94
Within-person variance proportion: 0.63

S1.3.2 Convergent Validity

The single-item EMA stress measure demonstrated significant convergent validity with the PSS across all analyses (Table S2).

Time-Matched Concurrent Validity: Aggregated EMA stress was significantly correlated with PSS at both Month 1 (r = 0.43, 95% CI [0.28, 0.55], p < .001) and Month 2 (r = 0.50, 95% CI [0.37, 0.62], p < .001). The overall correlation pooling both assessment periods was r = 0.46 (95% CI [0.36, 0.55], p < .001).

Baseline PSS vs. Overall EMA Stress: Baseline PSS was significantly correlated with overall aggregated EMA stress (r = 0.47, 95% CI [0.34, 0.58], p < .001), indicating that pre-study stress levels predicted daily stress experiences during the EMA period.

Table S2. Convergent Validity: Correlations Between EMA Stress and PSS

Month 1 (concurrent): n = 136, Mean EMA Days = 15.6, r = 0.43, 95% CI [0.28, 0.55], p < .001
Month 2 (concurrent): n = 140, Mean EMA Days = 13.5, r = 0.50, 95% CI [0.37, 0.62], p < .001
Overall (M1 + M2 pooled): n = 276, Mean EMA Days = 14.6, r = 0.46, 95% CI [0.36, 0.55], p < .001
Between-person: n = 168, r = 0.51, 95% CI [0.39, 0.62], p < .001
Baseline PSS vs. Overall EMA: n = 177, Mean EMA Days = 26.0, r = 0.47, 95% CI [0.34, 0.58], p < .001

The magnitude of these correlations is consistent with prior EMA validation studies. Siepe et al. (2025) reported r = 0.47 between a single-item EMA stress measure and baseline PSS-10 in a student sample, remarkably similar to our baseline-overall correlation of r = 0.47. Similarly, Siepe et al. (2025) found that mean daily depression correlated r = 0.56 with the sum score of the weekly PHQ-9, concluding that "the weekly measurement appears to approximate the daily scores fairly well." Faurholt-Jepsen et al. (2019) also demonstrated the validity of daily smartphone-based stress self-assessment by comparing it with the PSS and related constructs in healthy individuals.

The moderate correlation magnitudes are expected given differences in measurement modality (momentary vs. retrospective), temporal granularity (daily vs. monthly recall), and the known retrospective bias in traditional questionnaires compared to real-time EMA capture (Shiffman et al., 2008).

S1.4 Conclusion

The single-item daily stress measure demonstrated adequate psychometric properties for use in intensive longitudinal research. The ICC of 0.37 indicated that the item captures meaningful between-person differences while also being sensitive to within-person fluctuations—a desirable property for studying stress dynamics. The between-person reliability of 0.94 exceeded conventional thresholds when aggregating across multiple daily assessments. Convergent validity was supported by significant correlations with the PSS ranging from r = 0.43 to r = 0.51 across multiple analytic approaches, indicating that the EMA stress item captures construct-relevant variance consistent with established stress measures. These findings align with growing evidence supporting the psychometric adequacy of single-item EMA measures when appropriately evaluated using multilevel methods (Trull & Ebner-Priemer, 2024; Brose et al., 2020; Castro-Alvarez et al., 2024).

S1.5 References

Bolger, N., & Laurenceau, J. P. (2013). Intensive longitudinal methods: An introduction to diary and experience sampling research. Guilford Press.

Brose, A., Voelkle, M. C., Lövdén, M., Lindenberger, U., & Schmiedek, F. (2020). A coordinated analysis of variance in affect in daily life. Assessment, 27(8), 1683-1698. https://doi.org/10.1177/1073191118799460

Castro-Alvarez, S., Bringmann, L. F., Zajner, C., & Jeronimus, B. F. (2024). Assessing the internal consistency reliability of ecological momentary assessment measures: Insights from the WARN-D study. Psychological Assessment, 36(11-12), 738-753. https://doi.org/10.1037/pas0001410

Cohen, S., Kamarck, T., & Mermelstein, R. (1983). A global measure of perceived stress. Journal of Health and Social Behavior, 24(4), 385-396.

Cranford, J. A., Shrout, P. E., Iida, M., Rafaeli, E., Yip, T., & Bolger, N. (2006). A procedure for evaluating sensitivity to within-person change: Can mood measures in diary studies detect change reliably? Personality and Social Psychology Bulletin, 32(7), 917-929.

Diamantopoulos, A., Sarstedt, M., Fuchs, C., Wilczynski, P., & Kaiser, S. (2012). Single item measures in psychological science. European Journal of Psychological Assessment, 28(1), 1-8. https://doi.org/10.1027/1015-5759/a000699

Faurholt-Jepsen, M., Frost, M., Busk, J., Christensen, E. M., Bardram, J. E., Vinberg, M., & Kessing, L. V. (2019). The validity of daily self-assessed perceived stress measured using smartphones in healthy individuals: Cohort study. JMIR mHealth and uHealth, 7(8), e13418. https://doi.org/10.2196/13418

Shiffman, S., Stone, A. A., & Hufford, M. R. (2008). Ecological momentary assessment. Annual Review of Clinical Psychology, 4, 1-32.

Shrout, P. E., & Lane, S. P. (2012). Psychometrics. In M. R. Mehl & T. S. Conner (Eds.), Handbook of research methods for studying daily life (pp. 302-320). Guilford Press.

Siepe, B. S., Rieble, C. L., Tutunji, R., Rimpler, A., Marz, J., Proppert, R. K. K., & Fried, E. I. (2025). Understanding ecological-momentary-assessment data: A tutorial on exploring item performance in ecological-momentary-assessment data. Advances in Methods and Practices in Psychological Science, 8(1), 1-20. https://doi.org/10.1177/25152459241286877

Trull, T. J., & Ebner-Priemer, U. W. (2024). Evaluation of pressing issues in ecological momentary assessment. Annual Review of Clinical Psychology, 20, 107-135. https://doi.org/10.1146/annurev-clinpsy-080921-083128