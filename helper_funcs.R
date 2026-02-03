###################
# This functions serves to transform the raw data into the format ct-sem model expects
# Args:
#   1) raw_data (dataframe): dataframe used for model fitting
#   2) scale_iv_dv (boolean): whether to scale NDVI and stress (default is TRUE)
###################
data_prep <- function(raw_data, scale_iv_dv = TRUE){
  prepared_data <- raw_data %>%
    rename(
      id = Participant_ID,
      iv = Daily_TimeWeighted_NDVI,
      dv = stressed
    ) %>%
    # Convert dates to numeric time (days since first observation per person)
    mutate(Date = as.Date(Date)) %>%
    group_by(id) %>%
    arrange(Date) %>%
    mutate(
      time = as.numeric(Date - min(Date))
    ) %>%
    ungroup() %>%
    # Rename levels in Group column
    mutate(Group = forcats::fct_recode(Group, "CC" = "HC")) %>%
    
    # Handle demographics
    # Given the data mainly consists of white partitipants, 
    # we chose to binarize it into white versus non-white for later analysis
    mutate(
      race_numeric = ifelse(race == "White/European/North American/Australian", 0, 1)
    ) %>% 
    # Create baseline demographics (one value per person)
    group_by(id) %>%
    mutate(
      # Take first non-missing value for each person
      group_baseline = first(na.omit(Group)),
      sex_baseline = first(na.omit(sex)),
      age_baseline = first(na.omit(age)),
      education_baseline = first(na.omit(chrdemo_edu_max)),
      race_baseline = first(na.omit(race)),
      clgry_baseline = first(na.omit(clgry_total_baseline)),
      bprs_baseline = first(na.omit(bprs_total_baseline)),
      oasis_baseline = first(na.omit(oasis_total_baseline))
    ) %>%
    ungroup() %>%
    # Scale and center continuous variables
    mutate(
      iv = if (scale_iv_dv) as.numeric(scale(iv)) else as.numeric(iv),
      dv = if (scale_iv_dv) as.numeric(scale(dv)) else as.numeric(dv),
      
      # Continuous demographics
      age = as.numeric(scale(age_baseline)),
      education = as.numeric(scale(education_baseline)),
      clgry = as.numeric(scale(clgry_baseline)),
      bprs = as.numeric(scale(bprs_baseline)),
      oasis = as.numeric(scale(oasis_baseline))
    ) %>%
    mutate(
      # Code binary group and sex variables as numeric (0/1)
      group_numeric = ifelse(group_baseline == "CC", 0, 1),
      sex_numeric = ifelse(sex_baseline == "M", 0, 1)
    ) %>% 
    # (Only) mean-center binary variables 
    # Note: you could also not center the binary variables 
    # (in this case, the estimates would be the reference group)
    mutate(
      group = group_numeric - 0.5,
      sex = sex_numeric - 0.5,
      race = race_numeric - 0.5
    ) %>%
    # TD predictors
    mutate(
      day = scale(as.numeric(format(Date, "%w"))),  # Sunday as reference (0)
      hometime = scale(as.numeric(hometime)),
      entropy = scale(as.numeric(entropy)),
      locations = scale(as.numeric(locations)),
    ) %>%
    # Keep needed columns
    select(id, time, iv, dv,
           group, sex, age,
           race, education,
           clgry, bprs, oasis,
           day, hometime, entropy, locations) %>%

    # Remove rows with missing key variables
    # filter(!is.na(time), !is.na(id), !is.na(iv), !is.na(dv)) %>%
    arrange(id, time)

  return(prepared_data)
}

###################
# This function specifies the ct-sem model based on the types of interaction
# (i.e., no interaction, full interaction, or partial interaction of parameters)
# Args:
# 1) interaction (str): type of interaction ("none", "full", "partial")
# 2) TIpred (vector or str): time-independent moderators (covariates)
# 3) TDpred (vector or str): time-dependent predictors (control variables)

# Note, for partial interaction, it now includes drift effects, the latent
# process means at the first time point, and continuous-time intercepts to
# capture the main effects of the moderators. <= Could easily revise it based on needs
###################

model_spec <- function(interaction,
                       TIpred = c("group", "sex", "age",
                                  "race", "education"),
                       TDpred = NULL){
  # Define TDPREDEFFECT matrix if TDpred exists
  # If TDpred is NULL, this stays NULL
  td_effect_matrix <- NULL
  if(!is.null(TDpred)){
    # Rows = n.latent (2), Cols = n.TDpred (length of TDpred)
    td_effect_matrix <- matrix(
      paste0("td_", rep(c("NDVI", "Stress"), length(TDpred)), "_", 
             rep(TDpred, each = 2)), 
      nrow = 2, ncol = length(TDpred)
    )
  }

  model <- ctModel(
    type = 'ct', # Bayesian fitting with ctStanFit (using Stan)

    LAMBDA = diag(2),

    n.manifest = 2,
    n.latent = 2,

    # Manifest means (set to 0 since we scaled)
    MANIFESTMEANS = matrix(c(0, 0), nrow = 2),

    # Continuous intercepts (including random effects)
    CINT = matrix(c("cint_ndvi", "cint_stress"), nrow = 2),

    # DRIFT matrix with group interactions
    DRIFT = matrix(c(
      "ar_ndvi",         # NDVI autoregressive
      "cl_ndvi_stress",  # NDVI → Stress
      "cl_stress_ndvi",  # Stress → NDVI
      "ar_stress"        # Stress autoregressive
    ), nrow = 2, ncol = 2, byrow = TRUE),

    # Innovation diffusion matrix (lower triangular)
    DIFFUSION = matrix(c(
      "diffusion_ndvi", "0",
      "diffusion_corr", "diffusion_stress"
    ), nrow = 2, ncol = 2, byrow = TRUE),

    # Measurement error
    MANIFESTVAR = matrix(c(
      "resid_ndvi", "0",
      "0", "resid_stress"
    ), nrow = 2, byrow = TRUE),

    # Variable names
    manifestNames = c("iv", "dv"),
    latentNames = c("NDVI", "Stress"),

    # Time-independent predictors
    TIpredNames = TIpred,
    
    # Time-dependent predictors
    TDpredNames = TDpred,
    TDPREDEFFECT = td_effect_matrix,

    tipredDefault = (interaction == "full") # tipredDefault for full interaction of TI predictor effect
  )

  # Customize manual time independent predictor assignment for moderation effects
  if (interaction == "partial"){
    model$pars[, 9:ncol(model$pars)] = FALSE
    # Only include moderation on drift effects, the latent process means at the 
    # first time point, and continuous-time intercepts
    model$pars[c(1:2, 8:9, 21:22), 9:ncol(model$pars)] = TRUE
  }
  
  return(model)
}

###################
# This function extracts and formats time-dependent (TD) predictor effects from
# a ctsem model fit.
#
# Args:
# 1) fit (stanfit object): ctsem stanfit result
# 2) td_predictors (vector): names of TD predictors in the order they were specified
# 3) latent_names (vector): names of latent variables (default: c("NDVI", "Stress"))
# 4) output (str): output type (option: "table" or "text")
#
# Returns: A table showing how each TD predictor affects each latent variable
###################
extract_td_effects <- function(fit,
                               td_predictors,
                               latent_names = c("NDVI", "Stress"),
                               output = "table") {
  library(dplyr)
  library(knitr)

  # Get parameter matrices from summary
  params_matrices <- summary(fit)$parmatrices

  # Filter for TDPREDEFFECT
  td_effects <- params_matrices[params_matrices$matrix == "TDPREDEFFECT", ]

  if (nrow(td_effects) == 0) {
    warning("No TDPREDEFFECT found in model. Did you include TD predictors?")
    return(NULL)
  }

  # Create formatted table
  td_table <- data.frame(
    Latent = latent_names[as.numeric(td_effects$row)],
    TD_Predictor = td_predictors[as.numeric(td_effects$col)],
    Median = round(td_effects$`50%`, 3),
    SD = round(td_effects$sd, 3),
    CI_lower = round(td_effects$`2.5%`, 3),
    CI_upper = round(td_effects$`97.5%`, 3)
  )

  # Add significance column
  td_table$Significance <- ifelse(td_table$CI_lower * td_table$CI_upper > 0, "Sig", "Non-Sig")

  # Output the results
  if (output == "table") {
    return(td_table)
  } else {
    return(kable(td_table,
                 caption = "Time-Dependent Predictor Effects on Latent Variables",
                 row.names = FALSE))
  }
}

###################
# This function outputs a table of parameter estimates for the variable of interest, including
# posterior median, standard deviations (SD) and credibility intervals (CI) for
# means of estimated population distributions. 

# Args: 
# 1) fit (stanfit object): ctsem stanfit result
# 2) variable (str): variable of parameter esimates (options: "NDVI", "stress")
# 3) moderators (vector): moderators (covariates) in ctsem analysis
# 4) output (str): output type (option: "table" or "text")

###################
create_ctsem_table_1 <- function(fit, variable,
                                 moderators = c("group", "sex", "age", 
                                                "race", "education"),
                                output = "table") {
  library(dplyr)
  library(knitr)
  
  params_matrices <- summary(fit)$parmatrices
  tipreds <- summary(fit)$tipreds
  
  # Parameter estimates for latent process means at the first time point
  # and continuous-time intercepts for NDVI and stress
  t0_params <- params_matrices[params_matrices$matrix == "T0MEANS", ]
  cint_params <- params_matrices[params_matrices$matrix == "CINT", ]
  
  # Note: Currently not adding MAINIFEST because the manifest is set to zero 
  # (since we scale the variables)
  
  # Create structured output for the table
  table_1 <- data.frame(
    Parameter = character(),
    Median = numeric(),
    SD = numeric(),
    CI_lower = numeric(),
    CI_upper = numeric(),
    stringsAsFactors = FALSE
  )
  
  if (variable == "NDVI") index <- 1
  else index <- 2
  
  # Add T0 means
  table_1 <- rbind(table_1, data.frame(
    Parameter = "T0 mean",
    Median = t0_params[index, "50%"],
    SD = t0_params[index, "sd"],
    CI_lower = t0_params[index, "2.5%"],
    CI_upper = t0_params[index, "97.5%"]
  ))
  
  # Add CINT
  table_1 <- rbind(table_1, data.frame(
    Parameter = "CINT mean",
    Median = cint_params[index, "50%"],
    SD = cint_params[index, "sd"],
    CI_lower = cint_params[index, "2.5%"],
    CI_upper = cint_params[index, "97.5%"]
  ))
  
  # Add TI predictor effects (what we only need here is continuous-time
  # intercepts to capture the main effects of each moderator on variable of interest) 
  for(i in seq(moderators)){
    # Determine index position in `tipreds` dataframe
    if (variable == "NDVI") index <- 6 * i - 1
    else index <- 6 * i
    
    table_1 <- rbind(table_1, data.frame(
      Parameter = paste("Main effect of", moderators[i]),
      Median = tipreds[index, "50%"],
      SD = tipreds[index, "sd"],
      CI_lower = tipreds[index, "2.5%"],
      CI_upper = tipreds[index, "97.5%"]
    ))
  }
  
  # Add a column of significance based on CI_lower and CI_upper
  table_1$Significance <- ifelse(table_1$CI_lower * table_1$CI_upper > 0, "Sig", "Non-Sig")
  
  # Round all numeric columns to two digits
  table_1[c('Median', 'SD', 'CI_lower', 'CI_upper')] <- lapply(table_1[c('Median', 'SD', 'CI_lower', 'CI_upper')], round, 2)
  
  # Output the results
  if (output == "table"){
    table_1_output <- table_1
  } else {
    table_1_output <- kable(
      table_1,
      caption = paste0(
        "CTSEM Results: Posterior median, standard deviations (SD) and credibility 
      intervals (CI) for means of estimated", variable, "population distributions.")
    )
  }
  
  return(table_1_output)
}

###################
# This function outputs a table of autoregressive and cross-lagged effects 
# (both continuous and discrete), as well as diffusion matrices of NDVI and 
# stress.

# Args: 
# 1) fit (stanfit object): ctsem stanfit result
# 2) moderators (vector): moderators (covariates) in ctsem analysis
# 3) output (str): output type (option: "table" or "text")
###################

create_ctsem_table_2 <- function(fit,
                                 moderators = c("group", "sex", "age", 
                                                "race", "education"),
                                 output = "table") {
  
  library(dplyr)
  library(knitr)
  params_matrices <- summary(fit)$parmatrices
  
  # Create structured output for the second table
  table_2 <- data.frame(
    Type = character(),
    Parameter = character(),
    Median = numeric(),
    SD = numeric(),
    CI_lower = numeric(),
    CI_upper = numeric()
  )
  
  # Index drift parameters (both continuous and discrete) and diffusion parameters
  ct_drift_params <- params_matrices[params_matrices$matrix == "DRIFT", ]
  dt_drift_params <- params_matrices[params_matrices$matrix == "dtDRIFT", ]
  diffusion_covar_params <- params_matrices[params_matrices$matrix == "DIFFUSIONcov", ]
  
  if (nrow(ct_drift_params) == nrow(dt_drift_params) & nrow(dt_drift_params) == nrow(diffusion_covar_params)){
    for (i in 1:nrow(ct_drift_params)) {
      row_num <- ct_drift_params[i, "row"]
      col_num <- ct_drift_params[i, "col"]
      
      # Translate into NDVI or stress
      if (row_num == "1") row_var <- "NDVI" else row_var <- "Stress"
      if (col_num == "1") col_var <- "NDVI" else col_var <- "Stress"
      
      param_name <- paste(row_var, "➔", col_var)
      
      # Add continuous drift parameters
      table_2 <- rbind(table_2, data.frame(
        Type = "Continuous drift parameters",
        Parameter = param_name,
        Median = ct_drift_params[i, "50%"],
        SD = ct_drift_params[i, "sd"],
        CI_lower = ct_drift_params[i, "2.5%"],
        CI_upper = ct_drift_params[i, "97.5%"]
      ))
      
      # Add discrete drift parameters
      table_2 <- rbind(table_2, data.frame(
        Type = "Discrete drift parameters",
        Parameter = param_name,
        Median = dt_drift_params[i, "50%"],
        SD = dt_drift_params[i, "sd"],
        CI_lower = dt_drift_params[i, "2.5%"],
        CI_upper = dt_drift_params[i, "97.5%"]
      ))
      
      # Add diffusion parameters
      table_2 <- rbind(table_2, data.frame(
        Type = "Diffusion covariance",
        Parameter = param_name,
        Median = diffusion_covar_params[i, "50%"],
        SD = diffusion_covar_params[i, "sd"],
        CI_lower = diffusion_covar_params[i, "2.5%"],
        CI_upper = diffusion_covar_params[i, "97.5%"]
      ))
    }
  }
  
  # Reorder rows based on `Type` column
  table_2 <- table_2[order(
    factor(table_2$Type, 
           levels = c("Continuous drift parameters", 
                      "Discrete drift parameters", 
                      "Diffusion covariance"))
  ), ]

  # Add significance column
  table_2$Significance <- ifelse(table_2$CI_lower * table_2$CI_upper > 0, "Sig", "Non-Sig")
  
  # Round all numeric columns to two digits
  table_2[c('Median', 'SD', 'CI_lower', 'CI_upper')] <- lapply(table_2[c('Median', 'SD', 'CI_lower', 'CI_upper')], round, 2)
  
  # Output the results
  if (output == "table"){
    table_2_output <- table_2
  } else {
    table_2_output <- kable(
      table_2,
      caption = paste0(
        "CTSEM Results: Autoregressive and cross-lagged effects between NDVI and 
        Stress")
    )
  }
  
  return(table_2_output)
}

# ##############
# This function outputs long-term, asymptotic correlation between the 
# two processes in their random fluctuations.

# Args:
# 1) fit (stanfit object): ctsem stanfit result

# Note: some may calculate contemporaneous correlation between the two processes
# but we decided not to include that in our study as we were more interested in
# the long-term, steady-state relationship between the two processes. 
# If one would like to do contemporaneous correlation, they could do something 
# like the following:

# # Extract diffusion variance and covariance
# diffusion_var1 <- summary_fit$parmatrices$dtDIFFUSION[1,1]  # Process 1 variance
# diffusion_var2 <- summary_fit$parmatrices$dtDIFFUSION[2,2]  # Process 2 variance  
# diffusion_cov <- summary_fit$parmatrices$dtDIFFUSION[1,2]   # Covariance
# 
# # Calculate correlation
# contemporaneous_cor <- diffusion_cov / sqrt(diffusion_var1 * diffusion_var2)
# print(contemporaneous_cor)
# ##############
asymptotic_correlation <- function(fit){
  params_matrices <- summary(fit)$parmatrices
  
  # Extract mean, SD, and CI bounds
  asym_data <- params_matrices[params_matrices$matrix %in% "asymDIFFUSIONcov", ]
  
  # Create covariance matrices for mean and bounds
  asymDIFFUSIONcov_mean <- matrix(asym_data$Mean, nrow=2, ncol=2, byrow=TRUE)
  asymDIFFUSIONcov_lower <- matrix(asym_data$`2.5%`, nrow=2, ncol=2, byrow=TRUE) 
  asymDIFFUSIONcov_upper <- matrix(asym_data$`97.5%`, nrow=2, ncol=2, byrow=TRUE)
  
  # Convert to correlations
  cor_mean <- cov2cor(asymDIFFUSIONcov_mean)[1,2]  # Extract off-diagonal
  cor_lower <- cov2cor(asymDIFFUSIONcov_lower)[1,2]
  cor_upper <- cov2cor(asymDIFFUSIONcov_upper)[1,2]
  
  # Return formatted results
  cat("Asymptotic correlation between Greenspace Exposure and Perceived Stress:\n")
  cat(sprintf("r = %.2f, 95%% CI [%.2f, %.2f]\n", cor_mean, cor_lower, cor_upper))
} 

# ##############
# This function outputs moderated difference of cross-lagged effects in text.
# Essentially, the interpretation should be something like A ➔ B relationship 
# was stronger in Group X than in Group Y (Est. = [], SD = [], 95 % CI [; ]).

# Args:
# 1) fit (stanfit object): ctsem stanfit result
# 2) moderator (vector or string): grouping variables of interest (default is "group")
# ##############

cl_effect_moderate <- function(fit, moderator="group"){
  
  tipreds <- summary(fit)$tipreds
  
  # Create structured output for the table
  table <- data.frame(
    Parameter = character(),
    Median = numeric(),
    SD = numeric(),
    CI_lower = numeric(),
    CI_upper = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(i in seq(moderator)){
    # Determine index position in `tipreds` dataframe
    indices = c(6 * i - 3, 6 * i - 2)
    
    # Add NDVI ➔ Stress cross-lagged effects
    table <- rbind(table, data.frame(
      Parameter = paste("Moderated NDVI ➔ Stress cross-lagged effects of", moderator[i]),
      Median = tipreds[indices[1], "50%"],
      SD = tipreds[indices[1], "sd"],
      CI_lower = tipreds[indices[1], "2.5%"],
      CI_upper = tipreds[indices[1], "97.5%"]
    ))
    
    # Add Stress ➔ NDVI cross-lagged effects
    table <- rbind(table, data.frame(
      Parameter = paste("Moderated Stress ➔ NDVI cross-lagged effects of", moderator[i]),
      Median = tipreds[indices[2], "50%"],
      SD = tipreds[indices[2], "sd"],
      CI_lower = tipreds[indices[2], "2.5%"],
      CI_upper = tipreds[indices[2], "97.5%"]
    ))
  }
  
  table$Significance <- ifelse(table$CI_lower * table$CI_upper > 0, "Sig", "Non-Sig")
  
  output <- kable(table, digits = 2, caption = "Moderated Cross-Lagged Effects")
  
  return(output)
}

##############
# This function generates the Median and 95% CI of parameter estimates of two groups.

# Args:
# 1) fit (stanfit object): ctsem stanfit result
# 2) method (str): which method to derive parameters for each group 
#   (options: `ctStanTIpredeffects` and `ctExtract`) => `ctStanTIpredeffects` is a 
# convenient wrapper (from ctsem package) that automatically computes and 
# summarizes group-specific effects, while `ctExtract` gives you raw posterior 
# samples for custom Bayesian analyses (using rstan) but requires manual 
# calculation of the group effects. The bottom line is that the two methods returned
# almost identical results of parameter estimates. 
# 3) whichpars (vector or string): which matrices in ctsem to compute effects 
# of time independent predictors (default is the two cross-lagged effects)
# 4) whichTIpreds (vector or int): which of the tipreds in the fit object to calculate effects
# 5) timeinterval (numeric): positive numeric indicating time interval to use 
# for discrete time parameter matrices (default is NULL, indicating continuous time effects).
# 6) group_labels (str): labels for the two groups (defualt is CC and CHR)
# 7) probs (vector or numeric): quantile probabilities (defualt is 0.025, 0.5, 0.975)

##############
analyze_group_effects <- function(fit, method,
                                  whichpars = c('DRIFT[1,2]', 'DRIFT[2,1]'), 
                                  param_names = c('NDVI ➔ Stress', 'Stress ➔ NDVI'),
                                  whichTIpreds = 1, # group variable
                                  timeinterval = NULL,
                                  group_labels = c("CC", "CHR"),
                                  probs = c(0.025, 0.5, 0.975)) {
  
  if (method == 'ctStanTIpredeffects'){
    # Run ctStanTIpredeffects
    if(is.null(timeinterval)) {
      group_effects <- ctStanTIpredeffects(fit, 
                                           whichpars = whichpars,
                                           whichTIpreds = whichTIpreds,
                                           probs = probs)
    } else {
      group_effects <- ctStanTIpredeffects(fit, 
                                           whichpars = whichpars,
                                           whichTIpreds = whichTIpreds,
                                           timeinterval = timeinterval,
                                           probs = probs)
    }
    
    # Extract key group values (endpoints: -0.5 and +0.5)
    group_values <- group_effects$x[,1]
    min_idx <- which.min(group_values)  # Group 0 (-0.5)
    max_idx <- which.max(group_values)  # Group 1 (+0.5)
    
    # Create results table
    results_table <- data.frame()
    
    for(i in 1:length(param_names)) {
      # CC estimates
      group0_median <- group_effects$y[min_idx, i, "Quantile0.5"]
      group0_lower <- group_effects$y[min_idx, i, "Quantile0.025"]
      group0_upper <- group_effects$y[min_idx, i, "Quantile0.975"]
      
      # CHR estimates  
      group1_median <- group_effects$y[max_idx, i, "Quantile0.5"]
      group1_lower <- group_effects$y[max_idx, i, "Quantile0.025"]
      group1_upper <- group_effects$y[max_idx, i, "Quantile0.975"]
      
      # Add to table
      results_table <- rbind(results_table, 
                             data.frame(
                               Parameter = param_names[i],
                               Group = group_labels[1],
                               Median = group0_median,
                               CI_Lower = group0_lower,
                               CI_Upper = group0_upper,
                               Significant = ifelse(group0_lower * group0_upper > 0, "Sig", "Non-Sig"),
                               stringsAsFactors = FALSE
                             ))
      
      results_table <- rbind(results_table,
                             data.frame(
                               Parameter = param_names[i], 
                               Group = group_labels[2],
                               Median = group1_median,
                               CI_Lower = group1_lower,
                               CI_Upper = group1_upper,
                               Significant = ifelse(group1_lower * group1_upper > 0, "Sig", "Non-Sig"),
                               stringsAsFactors = FALSE
                             ))
    }
    # Run ctExtract (could be extended to add group differences as well)
  } else {
    # Extract posterior samples (load or save)
    extracted_params_path <- "output/extracted_params.rds"
    
    if (file.exists(extracted_params_path)) {
      cat("Loading extracted parameters from:", extracted_params_path, "\n")
      extracted_params <- readRDS(extracted_params_path)
    } else {
      cat("Extracting parameters from model...\n")
      
      extracted_params <- ctExtract(fit)
      saveRDS(extracted_params, extracted_params_path)
      cat("Saved extracted parameters to:", extracted_params_path, "\n")
    }
    
    # Extract baseline effects from pop_DRIFT (2x2 matrix)
    greenspace_to_stress_baseline <- extracted_params$pop_DRIFT[, 1, 2]  # DRIFT[1,2]
    stress_to_greenspace_baseline <- extracted_params$pop_DRIFT[, 2, 1]  # DRIFT[2,1]
    
    # Extract interaction effects from TIPREDEFFECT 
    # Based on popsetup in fit$setup: positions 4 and 5 for cross-lagged effects
    greenspace_to_stress_interaction <- extracted_params$TIPREDEFFECT[, 4, whichTIpreds]  # cl_ndvi_stress
    stress_to_greenspace_interaction <- extracted_params$TIPREDEFFECT[, 5, whichTIpreds]  # cl_stress_ndvi
    
    # Calculate group-specific effects with centered coding (-0.5, +0.5):
    # CC (group = -0.5): baseline + (-0.5) * interaction
    CC_green_to_stress <- greenspace_to_stress_baseline + (-0.5) * greenspace_to_stress_interaction
    CC_stress_to_green <- stress_to_greenspace_baseline + (-0.5) * stress_to_greenspace_interaction
    
    # CHR (group = +0.5): baseline + (+0.5) * interaction  
    chr_green_to_stress <- greenspace_to_stress_baseline + (0.5) * greenspace_to_stress_interaction
    chr_stress_to_green <- stress_to_greenspace_baseline + (0.5) * stress_to_greenspace_interaction
    
    # Combine samples for analysis
    CC_samples <- cbind(CC_green_to_stress, CC_stress_to_green)
    chr_samples <- cbind(chr_green_to_stress, chr_stress_to_green)
    
    # Calculate summary statistics
    results_table <- data.frame()
    
    for(i in 1:length(param_names)) {
      # CC estimates
      CC_median <- median(CC_samples[,i])
      CC_sd <- sd(CC_samples[,i])
      CC_lower <- quantile(CC_samples[,i], probs[1])
      CC_upper <- quantile(CC_samples[,i], probs[3])
      CC_sig <- ifelse(CC_lower * CC_upper > 0, "Sig", "Non-Sig")
      
      cat(
        "CC NDVI → Stress\n",
        "Median:", CC_median, "\n",
        "SD:", CC_sd, "\n",
        "CI Lower:", CC_lower, "\n",
        "CI Upper:", CC_upper, "\n",
        "Significance:", CC_sig, "\n\n"
      )
      
      results_table <- rbind(results_table, 
                             data.frame(
                               Parameter = param_names[i],
                               Group = group_labels[1],
                               Median = round(CC_median, 3),
                               SD = round(CC_sd, 3),
                               CI_Lower = round(CC_lower, 3),
                               CI_Upper = round(CC_upper, 3),
                               Significant = CC_sig
                             ))
      
      # CHR estimates  
      chr_median <- median(chr_samples[,i])
      chr_sd <- sd(chr_samples[,i])
      chr_lower <- quantile(chr_samples[,i], probs[1])
      chr_upper <- quantile(chr_samples[,i], probs[3])
      chr_sig <- ifelse(chr_lower * chr_upper > 0, "Sig", "Non-Sig")
      
      cat(
        "CHR NDVI → Stress\n",
        "Median:", chr_median, "\n",
        "SD:", chr_sd, "\n",
        "CI Lower:", chr_lower, "\n",
        "CI Upper:", chr_upper, "\n",
        "Significance:", chr_sig, "\n\n"
      )
      
      results_table <- rbind(results_table,
                             data.frame(
                               Parameter = param_names[i], 
                               Group = group_labels[2],
                               Median = round(chr_median, 3),
                               SD = round(chr_sd, 3),
                               CI_Lower = round(chr_lower, 3),
                               CI_Upper = round(chr_upper, 3),
                               Significant = chr_sig
                             ))
      
      rownames(results_table) <- NULL
    }
  }
  
  # Create output
  library(knitr)
  table_output <- kable(results_table, digits = 2, row.names = FALSE,
                        caption = ifelse(is.null(timeinterval), 
                                         "Group-Specific Continuous-Time Effects from CTSEM",
                                         paste0("Group-Specific Discrete-Time Effects (", timeinterval, " time unit) from CTSEM")))
    
    return(table_output)
}


##############
# This function plots discrete-time autoregressive and/or cross-lagged effects 
# of NDVI and stress.

# Args:
# 1) fit (stanfit object): ctsem stanfit result
# 2) indices (str): specifying type(s) of effect to plot (options: "all", "AR", 
#   "CR"). By default, plot both autoregressive and cross-lagged effects.
##############
plot_discrete_effects <- function(fit, indices = "all"){
  library(ggplot2)
  
  if (indices == "all") {
    effects <- "autoregressive and cross-lagged effects"
  } else if (indices == "AR") {
    effects <- "autoregressive effects"
  } else {
    effects <- "cross-lagged effects"
  }
  
  discrete_pars <- ctStanDiscretePars(fit, plot = FALSE) 
  
  p <- ctStanDiscreteParsPlot(
    discrete_pars,
    indices = indices,
    latentNames = c("GE", "PS"), 
    xlab = "Time interval (days)",
    ylab = "Drift effect estimate",
    title = paste("Discrete-time", effects, "\nbetween greenspace exposure (GE) and perceived stress (PS)\n")
  ) + 
    # Replace dots with arrows in legend labels
    scale_color_discrete(
      labels = function(x) gsub("\\.", " →", x)  
    ) +
    scale_fill_discrete(
      labels = function(x) gsub("\\.", " →", x)
    ) +
    theme(
      legend.key.height = unit(0.6, "cm"),  # This controls spacing for multi-line labels
      legend.text = element_text(lineheight = 1.2),  # Line spacing within each label
      panel.border = element_blank(),           # Remove panel border
      axis.line.x = element_line(color = "black"),  # Add x-axis line
      axis.line.y = element_line(color = "black"),  # Add y-axis line
      panel.grid.major = element_line(color = "grey90"),  # Keep major grid lines
      panel.grid.minor = element_blank()        # Remove minor grid lines
    )
  
  # Save the plot
  ggsave(paste0("output/discrete_time_plot_", indices,".png"))
  
  return(p)
}

##############
# This function plots discrete-time effects of NDVI and stress for the two diagnostic group.

# Args:
# 1) fit (stanfit object): ctsem stanfit result
# 2) whichTIpreds (vector or int): which of the tipreds in the fit object to 
# calculate effects. By default, set to 1 (the index position of diagnositic 
# group when specifying the TIpreds). 
# 3) nsamples (int): number of samples from the stanfit to use for analysi and 
# plotting. The default value is the same as ctStanDiscretePars, which is 200. 
# Higher values will increase smoothness / accuracy, at cost of plotting speed. 
# Values greater than the total number of samples will be set to total samples.
# 4) indices (str): specifying type(s) of effect to plot (options: "all", "AR", 
#   "CR"). By default, plot both autoregressive and cross-lagged effects.
# 5) plot_style (str): the style of each group's respective discrete-time plot
#   (options: "facet" - one plot with two facets, "separate" - two separate 
#   plots). By default, plot one plot with two facets. 
##############
plot_discrete_effects_group <- function(fit, whichTIpreds = 1, nsamples = 200, 
                                        indices = "all", plot_style = 'facet'){
  
    # Extract population-level parameters
    extracted_params <- ctExtract(fit, nsamples = nsamples, subjects = "popmean")
    
    niter = dim(extracted_params$pop_DRIFT)[1]
    if(nsamples > niter) nsamples <- niter
    
    sample_indices <- sample(1:niter, nsamples)
    
    group_results <- list()
    discrete_results <- list()
    
    # Define group values and labels
    # CC (group = -0.5): baseline + (-0.5) * interaction
    # CHR (group = 0.5): baseline + (0.5) * interaction
    group_values <- c(-0.5, 0.5)
    group_labels <- c("CC", "CHR")
    
    for(group_idx in 1:2) {
      group_val <- group_values[group_idx]
      
      # Population-level baseline matrices
      baseline_drift <- extracted_params$pop_DRIFT[sample_indices, , , drop = FALSE]
      baseline_diffcov <- extracted_params$pop_DIFFUSIONcov[sample_indices, , , drop = FALSE]
      baseline_asymdiffcov <- extracted_params$pop_asymDIFFUSIONcov[sample_indices, , , drop = FALSE]
      
      # TI interaction effects  
      tipred_effects <- extracted_params$TIPREDEFFECT[sample_indices, , whichTIpreds, drop = FALSE]
      
      # Create group-specific drift matrix
      group_drift <- baseline_drift
      group_drift[, 1, 2] <- baseline_drift[, 1, 2] + group_val * tipred_effects[, 4, ]  # GE -> PS
      group_drift[, 2, 1] <- baseline_drift[, 2, 1] + group_val * tipred_effects[, 5, ]  # PS -> GE
      
      # Package for ctStanDiscreteParsDrift (add subject dimension)
      ctpars <- list(
        DRIFT = array(group_drift, dim = c(nsamples, 1, 2, 2)),
        DIFFUSIONcov = array(baseline_diffcov, dim = c(nsamples, 1, 2, 2)),
        asymDIFFUSIONcov = array(baseline_asymdiffcov, dim = c(nsamples, 1, 2, 2))
      )
      
      group_results[[group_idx]] <- ctpars
      
      # Pass to ctStanDiscreteParsDrift (see function below)
      # (All following the default options in ctStanDiscretePars)
      discrete_effects <- ctStanDiscreteParsDrift(
        ctpars = ctpars,
        times = seq(from = 0, to = 10, by = 0.1),
        observational = FALSE,
        standardise = FALSE,
        cov = FALSE,
        discreteInput = fit$ctstanmodel$continuoustime == FALSE,
        quiet = TRUE
      )
      
      # Add proper dimnames for plotting compatibility
      dimnames(discrete_effects) <- list(
        Sample = sample_indices,
        Subject = "popmean",
        `Time interval` = seq(from = 0, to = 10, by = 0.1),
        row = fit$ctstanmodel$latentNames,
        col = fit$ctstanmodel$latentNames
      )
      
      # Set required attributes for plotting
      attributes(discrete_effects)$observational <- FALSE
      attributes(discrete_effects)$cov <- FALSE
      
      discrete_results[[group_idx]] <- discrete_effects
    }
      
    # Do the plotting
    if (indices == "all") {
      effects <- "autoregressive and cross-lagged effects"
    } else if (indices == "AR") {
      effects <- "autoregressive effects"
    } else {
      effects <- "cross-lagged effects"
    }
    
    if (plot_style == "facet"){
      # FACETED VERSION
      # Combine discrete results before calling ctStanDiscreteParsPlot
      combined_discrete <- abind::abind(
        discrete_results[[1]],
        discrete_results[[2]],
        along = 2  # Along subject dimension
      )
      
      # Set proper dimnames for 5D array [samples, subjects, times, row, col]
      dimnames(combined_discrete) <- list(
        Sample = sample_indices,
        Subject = group_labels,  # CC, CHR
        `Time interval` = seq(from = 0, to = 10, by = 0.1),
        row = fit$ctstanmodel$latentNames,
        col = fit$ctstanmodel$latentNames
      )
      
      # Set required attributes
      attributes(combined_discrete)$observational <- FALSE
      attributes(combined_discrete)$cov <- FALSE
      
      # Create the faceted overlay plot
      faceted_plot <- ctStanDiscreteParsPlot(
        x = combined_discrete,
        indices = indices,
        latentNames = c("GE", "PS"),
        splitSubjects = TRUE,
        facets = "Subject",  # This creates separate panels for each group
        xlab = "Time interval (days)",
        ylab = "Drift effect estimate", 
        title = paste("Discrete-time", effects, "between greenspace exposure (GE) and perceived stress (PS)\n")
      ) + 
        
        # Custom styling
        scale_color_discrete(
          name = "Effect",
          labels = function(x) gsub("\\.", " →", x)
        ) +
        
        scale_fill_discrete(
          name = "Effect", 
          labels = function(x) gsub("\\.", " →", x)
        ) +
        
        # Customize facet labels to show full group names
        facet_wrap(
          ~ factor(Subject,
                   levels = c("CHR", "CC"),
                   labels = c("Clinical High-Risk (CHR)", "Community Control (CC)")),
          scales = "fixed") + # Allow different y-scales if effects differ greatly
    
        theme(
          legend.key.height = unit(0.6, "cm"),
          legend.text = element_text(lineheight = 1.2),
          panel.border = element_blank(),
          axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black"),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank(),
          panel.spacing.x = unit(0.5, "cm"), # Horizontal spacing between two facets
          plot.title = element_text(hjust = 0.5, size = 12),
          
          # Facet styling
          strip.background.x = element_rect(colour = "black", fill = "white", size = 1.5, linetype = "solid"),
          strip.text = element_text(size = 11, face = "bold"),
          strip.text.x = element_text(margin = margin(t = 0.3, b = 0.3, unit = "cm")),
          
          # Legend positioning
          legend.position = "right",
          legend.box = "horizontal"
        )
      
      # Print the plot
      print(faceted_plot)
      
      # Save the faceted plot
      ggsave(paste0("output/discrete_time_plot_", indices, "_faceted.png"),
             plot = faceted_plot, width = 12, height = 6, dpi = 300)
    } else {
        # SEPARATED VERSION
        for(group_idx in 1:2) {
          p <- ctStanDiscreteParsPlot(
            x = discrete_results[[group_idx]],
            indices = "all",
            latentNames = c("GE", "PS"),
            xlab = "Time interval (days)",
            ylab = "Drift effect estimate",
            title = paste("Discrete-time", effects, "between greenspace exposure (GE) and perceived stress (PS) for", group_labels[group_idx], "\n")
          )
  
          # Add custom styling to match your original function
          p <- p +
            scale_color_discrete(labels = function(x) gsub("\\.", " →", x)) +
            scale_fill_discrete(labels = function(x) gsub("\\.", " →", x)) +
            theme(
              legend.key.height = unit(0.6, "cm"),
              legend.text = element_text(lineheight = 1.2),
              panel.border = element_blank(),
              axis.line.x = element_line(color = "black"),
              axis.line.y = element_line(color = "black"),
              panel.grid.major = element_line(color = "grey90"),
              panel.grid.minor = element_blank()
            )
          
          # Print the plot
          print(p)
  
          # Save individual plots
          ggsave(paste0("output/discrete_time_plot_", indices, "_", group_labels[group_idx], ".png"),
                 plot = p, width = 10, height = 6, dpi = 300)
        }
    }
}

##############
# This function is written by the authors of the package in their source code for 
# ctDiscretePars (https://rdrr.io/cran/ctsem/src/R/ctDiscretePars.R). 
##############
ctStanDiscreteParsDrift<-function(ctpars,times, observational,  standardise,cov=FALSE,
                                  types='dtDRIFT',discreteInput=FALSE, quiet=FALSE){
  
  nl=dim(ctpars$DRIFT)[3]
  
  if(!quiet) message('Computing temporal regression coefficients for ', dim(ctpars$DRIFT)[1],' samples')
  
  lapply(names(ctpars),function(x){ #add in extra dim if only 3 dims (e.g. when not individually varying)
    dm=dim(ctpars[[x]])
    if(length(dm)==3){
      ctpars[[x]] <<- array(ctpars[[x]],dim=c(dm[1],1,dm[2:3]))
    }
  })
  
  nsubs <- lapply(ctpars,function(x) dim(x)[2])
  
  
  if('dtDRIFT' %in% types){ 
    ctpars$dtDRIFT <- array(NA, dim=c(dim(ctpars$DRIFT)[1],max(unlist(nsubs)),length(times),dim(ctpars$DRIFT)[3:4]))
    
    mpow <- function(m,n){
      if(n==0) return(diag(1,nrow(m))) else{
        if(n>1){
          mo <-m
          for(i in 2:n){
            m <- m %*% mo
          }
        }
        return(m)
      }}
    
    for(i in 1:dim(ctpars$DRIFT)[1]){
      for(j in 1:dim(ctpars$DRIFT)[2]){
        for(ti in 1:length(times)){
          if(!discreteInput) ctpars$dtDRIFT[i,j,ti,,] <- expm::expm(as.matrix(ctpars$DRIFT[i,min(j,nsubs$DRIFT),,] * times[ti]))
          if(discreteInput) ctpars$dtDRIFT[i,j,ti,,] <- mpow(as.matrix(ctpars$DRIFT[i,min(j,nsubs$DRIFT),,]),times[ti])
          if(standardise) {
            if(any(diag(ctpars$asymDIFFUSIONcov[i,min(j,nsubs$asymDIFFUSIONcov),,]) < 0)) stop(
              "Asymptotic diffusion matrix has negative diagonals -- I don't know what non stationary standardization looks like")
            ctpars$dtDRIFT[i,j,ti,,] <- ctpars$dtDRIFT[i,j,ti,,] * 
              matrix(rep(sqrt(diag(ctpars$asymDIFFUSIONcov[i,min(j,nsubs$asymDIFFUSIONcov),,])+1e-10),each=nl) / 
                       rep((sqrt(diag(ctpars$asymDIFFUSIONcov[i,min(j,nsubs$asymDIFFUSIONcov),,]))),times=nl),nl)
          }
          if(observational){
            Qcor<-cov2cor(matrix(ctpars$DIFFUSIONcov[i,min(j,nsubs$DIFFUSIONcov),,],nl,nl)+diag(1e-8,nl)) 
            Qcor <- Qcor #* sign(Qcor) #why was this squared before?
            # browser()
            ctpars$dtDRIFT[i,j,ti,,]  <- ctpars$dtDRIFT[i,j,ti,,]  %*% Qcor
          }
          if(cov) ctpars$dtDRIFT[i,j,ti,,]  <- tcrossprod(ctpars$dtDRIFT[i,j,ti,,] )
        }
      }
    }
  } #end dtdrift
  
  
  return(ctpars$dtDRIFT)
}

# create_ctsem_table_2 => Extract matrices of interest
# analyze_group_effects => Extract matrices (parameters) of interest for each group
# ctStanDiscreteParsDrift => Get discrete-time effects based on provided matrices

