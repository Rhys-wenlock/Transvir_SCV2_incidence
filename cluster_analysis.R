require(lmtest)
require(broom)
require(tidyverse)
require(stringr)
require(stats)
require(lme4)
require(Matrix)

# Univariate --------------------------------------------------------------

#requires dataset = Total_single_cluster_data

#Merge datasets
total_single_cluster_data <- total_single_cluster_data %>%
  left_join(full_demog_data %>% select(Participant_ID, Total_household_no, age_cat, 
                                       Children_5_17, Children_2_4, Children_6m_2, Children_under_6m) %>% 
              rename(Index_ID = Participant_ID), by = "Index_ID") %>% 
  mutate(inf = case_when(
    PCR_result=="Positive" ~1, 
    PCR_result=="Negative"~0, 
    TRUE~NA
  )) %>% mutate(hh_size = as.numeric(Total_household_no)) %>%
  mutate(total_child = Children_5_17 + Children_2_4 + Children_6m_2 + Children_under_6m) %>%
  select(-Children_5_17, -Children_2_4, -Children_6m_2, -Children_under_6m) %>%
  mutate(total_child = as.numeric(total_child)) %>%
  mutate(age_cat = as.factor(age_cat)) %>%
  filter(!is.na(inf))


# List of columns
cols_of_interest <- c("age_cat", "period",
                      "hh_size", "total_child")


# Loop over the columns of interest

results_list <- list()

for (current_col in cols_of_interest) {
  
  # Fit the model
  current_model <- glm(
    formula = as.formula(paste0("inf ~ ", current_col)),
    family = binomial(),
    data = total_single_cluster_data
  )
  
  # Fit the null model without the current column
  null_model <- glm(
    formula = inf ~ 1,
    family = binomial(),
    data = total_single_cluster_data
  )
  
  # Conduct likelihood ratio test
  lrt_result <- anova(null_model, current_model, test = "LRT")
  lrt_pvalue <- lrt_result$`Pr(>Chi)`[2]  # Extracting the p-value of the test
  
  # Calculate hazard ratios
  model_results <- data.frame(
    variable_level = row.names(coef(summary(current_model))),
    or = exp(coef(current_model)),
    lower_ci_or = exp(confint(current_model)[, 1]),
    upper_ci_or = exp(confint(current_model)[, 2]),
    p_value = lrt_pvalue  # Add p-value to the model results dataframe
  )
  
  # Remove the column prefix for clarity
  model_results$variable_level <- str_replace(
    model_results$variable_level,
    paste0(current_col, ""),
    ""
  )
  
  # ... (previous code)
  
  if(is.numeric(total_single_cluster_data[[current_col]])) { 
    # Continuous variable
    current_results_table <- total_single_cluster_data %>%
      summarize(
        variable_level = as.character(current_col),
        num_exposed = paste0(median(!!sym(current_col)), " (", quantile(!!sym(current_col), 0.25), "-", quantile(!!sym(current_col), 0.75), ")"),
        prop_exposed = n() / 468,
        num_inf = sum(inf == 1), 
        prop_inf = (num_inf / n()) * 100
      )
  } else {
    # Categorical variable
    current_results_table <- total_single_cluster_data %>%
      group_by(!!sym(current_col)) %>%
      summarize(
        num_exposed = as.character(n()), # Convert to character here
        prop_exposed = n() / 468,
        num_inf = sum(inf == 1), 
        prop_inf = (num_inf / n()) * 100
      ) %>%
      mutate(variable_level = as.character(!!sym(current_col))) # Convert to character here too
  }
  
  # ... (rest of the code)
  
  
  # Common operations regardless of variable type
  current_results_table <- current_results_table %>%
    left_join(model_results, by = "variable_level") %>%
    mutate(OR = ifelse(is.na(or), NA, paste0(round(or, 2), "(", round(lower_ci_or, 2), "-", round(upper_ci_or, 2), ")"))) %>%
    select(variable_level, num_exposed, prop_exposed, num_inf, prop_inf, OR, p_value)  # Add p_value to the results table
  
  # Add the identifier column
  current_results_table$variable <- current_col
  
  # Append the table to the list
  results_list[[current_col]] <- current_results_table
}


for (current_col in cols_of_interest) {
  
  # Check if the current column is numeric
  if(is.numeric(total_single_cluster_data[[current_col]])) {
    # Continuous variable
    formula_str <- paste0("inf ~ ", current_col, " + (1|Household_ID)")
  } else {
    # Categorical variable
    formula_str <- paste0("inf ~ ", current_col, " + (1|Household_ID)")
  }
  
  # Fit the model
  current_model <- glmer(
    formula = as.formula(formula_str),
    family = binomial(),
    data = total_single_cluster_data
  )
  
  # Fit the null model without the current column
  null_model <- glmer(
    formula = inf ~ 1 + (1|Household_ID),
    family = binomial(),
    data = total_single_cluster_data
  )
  
  # Conduct likelihood ratio test
  lrt_result <- anova(null_model, current_model, test = "LRT")
  lrt_pvalue <- lrt_result$`Pr(>Chi)`[2]  # Extracting the p-value of the test
  
  # ...
  
  # Calculate hazard ratios
  fixed_effects <- fixef(current_model)
  
  if(is.factor(total_single_cluster_data[[current_col]])) {
    # Get the reference level for the factor
    ref_level <- levels(total_single_cluster_data[[current_col]])[1]
    
    # Exclude the reference level for the factor
    factor_levels <- setdiff(levels(total_single_cluster_data[[current_col]]), ref_level)
    
    # Extract the coefficients for the factor levels
    factor_coefs <- fixed_effects[names(fixed_effects) %in% paste0(current_col, factor_levels)]
    
    if(length(factor_coefs) != length(factor_levels)) {
      stop("Mismatch between factor levels and extracted coefficients.")
    }
    
    model_results <- data.frame(
      variable_level = factor_levels,
      or = exp(factor_coefs),
      lower_ci_or = exp(sapply(factor_coefs, function(x) confint(current_model, parm = names(x))[1])),
      upper_ci_or = exp(sapply(factor_coefs, function(x) confint(current_model, parm = names(x))[2])),
      p_value = rep(lrt_pvalue, length(factor_levels))
    )
    
  } else {
    predictor_coef <- fixed_effects[current_col]
    
    if(length(predictor_coef) == 0) {
      stop(paste("Coefficient for", current_col, "not found in the model."))
    }
    
    model_results <- data.frame(
      variable_level = current_col,
      or = exp(predictor_coef),
      lower_ci_or = exp(confint(current_model, parm = current_col)[1]),
      upper_ci_or = exp(confint(current_model, parm = current_col)[2]),
      p_value = lrt_pvalue
    )
  }
  
  # ...
  
  
  
  
  # Remove the column prefix for clarity
  model_results$variable_level <- str_replace(
    model_results$variable_level,
    paste0(current_col, ""),
    ""
  )
  if(is.numeric(total_single_cluster_data[[current_col]])) { 
    # Continuous variable
    current_results_table <- total_single_cluster_data %>%
      summarize(
        variable_level = as.character(current_col),
        num_exposed = paste0(median(!!sym(current_col)), " (", quantile(!!sym(current_col), 0.25), "-", quantile(!!sym(current_col), 0.75), ")"),
        prop_exposed = n() / 468,
        num_inf = sum(inf == 1), 
        prop_inf = (num_inf / n()) * 100
      )
  } else {
    # Categorical variable
    current_results_table <- total_single_cluster_data %>%
      group_by(!!sym(current_col)) %>%
      summarize(
        num_exposed = as.character(n()), # Convert to character here
        prop_exposed = n() / 468,
        num_inf = sum(inf == 1), 
        prop_inf = (num_inf / n()) * 100
      ) %>%
      mutate(variable_level = as.character(!!sym(current_col))) # Convert to character here too
  }
  
  # ... (rest of the code)
  
  
  # Common operations regardless of variable type
  current_results_table <- current_results_table %>%
    left_join(model_results, by = "variable_level") %>%
    mutate(OR = ifelse(is.na(or), NA, paste0(round(or, 2), "(", round(lower_ci_or, 2), "-", round(upper_ci_or, 2), ")"))) %>%
    select(variable_level, num_exposed, prop_exposed, num_inf, prop_inf, OR, p_value)  # Add p_value to the results table
  
  # Add the identifier column
  current_results_table$variable <- current_col
  
  # Append the table to the list
  results_list[[current_col]] <- current_results_table
}

null_model <- glmer(inf ~1 + (1|Household_ID), family=binomial(), total_single_cluster_data)

#Age_cat univariate
age_model <- glmer(inf ~ age_cat + (1|Household_ID), family=binomial(), total_single_cluster_data)
age_fef <- exp(fixef(age_model))
age_cint <- exp(confint(age_model))

results_list$age_cat$OR[2] <- 
  paste0(round(age_fef[2],2), " (", round(age_cint[3],2), "-", round(age_cint[8],2), ")")
results_list$age_cat$OR[3] <- 
  paste0(round(age_fef[3],2), " (", round(age_cint[4],2), "-", round(age_cint[9],2), ")")
results_list$age_cat$OR[4] <- 
  paste0(round(age_fef[4],2), " (", round(age_cint[5],2), "-", round(age_cint[10],2), ")")

p_val <- lrtest(age_model, null_model)
results_list$age_cat$p_value <- p_val$`Pr(>Chisq)`[2]

#ADDITIONAL Household Size Code
hh_size_model <- glm(inf ~ hh_size, family = binomial, data=total_single_cluster_data)

# Get the coefficient of the predictor (excluding the intercept)
pred_coef <- exp(coef(hh_size_model)[-1]) 
pred_conf_ints <- exp(confint(hh_size_model)[-1, ])

# Construct the OR string
OR_values <- paste0(round(pred_coef,2), " (", round(pred_conf_ints[1],2), "-", round(pred_conf_ints[2],2), ")")

# Assign the OR value
results_list$hh_size$OR <- OR_values

P_val <- lrtest(hh_size_model)
P_val <- P_val$`Pr(>Chisq)`[2]
results_list$hh_size$p_value <- P_val

#ADDITIONAL Child Size Code
child_size_model <- glm(inf ~ total_child, family = binomial, data=total_single_cluster_data)

# Get the coefficient of the predictor (excluding the intercept)
pred_coef <- exp(coef(child_size_model)[-1]) 
pred_conf_ints <- exp(confint(child_size_model)[-1, ])

# Construct the OR string
OR_values <- paste0(round(pred_coef,2), " (", round(pred_conf_ints[1],2), "-", round(pred_conf_ints[2],2), ")")

# Assign the OR value
results_list$total_child$OR <- OR_values

P_val <- lrtest(child_size_model)
P_val <- P_val$`Pr(>Chisq)`[2]
results_list$total_child$p_value <- P_val


# Bind all tables together into a single dataframe
results_table <- bind_rows(results_list)

