require(lmtest)
require(broom)
require(tidyverse)
require(matrix)
require(stringr)
require(stats)
require(lme4)
#requires lexis_inf_split

# Univariate --------------------------------------------------------------


#Set up variables
lexis_inf_split$prior_infection[lexis_inf_split$prior_infection==3] <- 2
lexis_inf_split$prior_infection[lexis_inf_split$prior_infection==4] <- 2
lexis_inf_split$prior_infection <- factor(lexis_inf_split$prior_infection, levels=c("0", "1", "2"))
lexis_inf_split$HIV <- as.factor(lexis_inf_split$HIV)
lexis_inf_split$HTN <- as.factor(lexis_inf_split$HTN)
lexis_inf_split$Diabetes <- as.factor(lexis_inf_split$Diabetes)
lexis_inf_split$Employed <- as.factor(lexis_inf_split$Employed)
lexis_inf_split$Steroid <- as.factor(lexis_inf_split$Steroid)
lexis_inf_split$Smoking <- as.factor(lexis_inf_split$Smoking)
lexis_inf_split$vaccinated <- as.factor(lexis_inf_split$vaccinated)
lexis_inf_split$Cancer <- as.factor(lexis_inf_split$Cancer)
symp_data <- lexis_inf_split[lexis_inf_split$lex.Xst==TRUE & lexis_inf_split$type=="PCR_pos",]
symp_data$age_cat <- factor(symp_data$age_cat, levels=c("18-49", "<5", "5-17", ">50"))
symp_data$prior_infection <- factor(symp_data$prior_infection, levels=c("0", "1", "2"))
symp_data <- symp_data %>%
  mutate(symp = case_when(
    symp==2 ~ 1, 
    symp==0 ~ 0, 
    TRUE ~ NA
  ))

# List of columns
cols_of_interest <- c("prior_infection", "period", "age_cat", 
                  "hh_size", "HIV", "Employed", "Steroid",
                  "Cancer", "Diabetes", "HTN", "Smoking",
                  "vaccinated")



# Loop over the columns of interest
results_list <- list()
for (current_col in cols_of_interest) {
  
  # Fit the model
  current_model <- glm(
    formula = as.formula(paste0("symp ~ ", current_col)),
    family = binomial(),
    data = symp_data[!is.na(symp_data$symp),]
  )
  
  # Fit the null model without the current column
  null_model <- glm(
    formula = symp ~ 1,
    family = binomial(),
    data = symp_data[!is.na(symp_data$symp),]
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
  
  # Create the results table
  current_results_table <- symp_data[!is.na(symp_data$symp),] %>%
    group_by(!!sym(current_col)) %>%
    summarize(
      num_infections = length(ID),
      prop_infections = num_infections / 263,
      num_symp = sum(symp == 1), 
      prop_symp = (num_symp / num_infections) * 100
    ) %>%
    rename(variable_level = !!sym(current_col)) %>%
    left_join(model_results, by = "variable_level") %>%
    mutate(OR = paste0(round(or, 2), "(", round(lower_ci_or, 2), "-", round(upper_ci_or, 2), ")")) %>%
    select(variable_level, num_infections, prop_infections, num_symp, prop_symp, OR, p_value)  # Add p_value to the results table
  
  # Add the identifier column
  current_results_table$variable <- current_col
  
  # Append the table to the list
  results_list[[current_col]] <- current_results_table
}

# Bind all tables together into a single dataframe
results_table <- bind_rows(results_list)
results_table

write.csv(results_table, "Symp_univariate.csv")

# Multivariate Modelling --------------------------------------------------

adjust_model <- glm(symp ~ age_cat + prior_infection + period, data=symp_data[symp_data$symp!=2,])

# Age ---------------------------------------------------------------------

#Produce null model
#test age cat
agenull_model <- glm(symp ~ prior_infection + period, data = symp_data[symp_data$symp!=2,])
#Crude model
age_model <- glm(symp ~age_cat, data = symp_data[symp_data$symp!=2,])

# Calculate p-values from likelihood ratio tests
p_val_crude_age <- lrtest(age_model)$Pr[2]
p_val_adjusted <- lrtest(adjust_model, agenull_model)$Pr[2]

# tidy() from the broom package extracts the necessary information from the models
tidy_age_model <- tidy(age_model, conf.int = TRUE) %>% 
  mutate(model = "Crude Model",
         p_value = p_val_crude_age)

tidy_adjust_model <- tidy(adjust_model, conf.int = TRUE) %>% 
  mutate(model = "Adjusted Model",
         p_value = p_val_adjusted)

# Join the results from both models
results_table_adjage <- rbind(tidy_age_model, tidy_adjust_model)

# Calculate hazard ratios and confidence intervals
results_table_adjage <- results_table_adjage %>%
  mutate(
    hr = exp(estimate),
    lower_ci = exp(conf.low),
    upper_ci = exp(conf.high),
    conf_int = paste0("(", round(lower_ci, 2), "-", round(upper_ci, 2), ")")
  )

results_table_adjage <- results_table_adjage[grepl("^age_cat", results_table_adjage$term), ]


# Select and rename necessary columns
results_table_adjage <- results_table_adjage %>%
  select(model, term, hr, lower_ci, upper_ci, conf_int, p_value) %>%
  rename(
    Model = model,
    `Age Category` = term,
    `Hazard Ratio` = hr,
    `Lower 95% CI` = lower_ci,
    `Upper 95% CI` = upper_ci,
    `95% CI` = conf_int,
    `p-value` = p_value
  ) %>%
  mutate(HR = paste0(round(`Hazard Ratio`, 2), "(", round(`Lower 95% CI`, 2), "-", round(`Upper 95% CI`, 2), ")")) %>%
  select(`Model`, `Age Category`,`HR`, `p-value`) %>%
  pivot_wider(
    names_from = Model,
    values_from = c(`HR`, `p-value`))


print(results_table_adjage)



# Prior Infection ---------------------------------------------------------

#Produce null model
#test prior inf cat
priorinfnull_model <- glm(symp ~ age_cat + period, data = symp_data[symp_data$symp!=2,])
#Crude model
priorinf_model <- glm(symp ~prior_infection, data = symp_data[symp_data$symp!=2,])
# Calculate p-values from likelihood ratio tests
p_val_crude_priorinf <- lrtest(priorinf_model)$Pr[2]
p_val_adjusted <- lrtest(adjust_model, priorinfnull_model)$Pr[2]

# tidy() from the broom package extracts the necessary information from the models
tidy_pinf_model <- tidy(priorinf_model, conf.int = TRUE) %>% 
  mutate(model = "Crude Model",
         p_value = p_val_crude_priorinf)

tidy_adjust_model <- tidy(adjust_model, conf.int = TRUE) %>% 
  mutate(model = "Adjusted model",
         p_value = p_val_adjusted)


# Join the results from both models
results_table_adjpinf <- rbind(tidy_pinf_model, tidy_adjust_model)

# Calculate hazard ratios and confidence intervals
results_table_adjpinf <- results_table_adjpinf %>%
  mutate(
    hr = exp(estimate),
    lower_ci = exp(conf.low),
    upper_ci = exp(conf.high),
    conf_int = paste0("(", round(lower_ci, 2), "-", round(upper_ci, 2), ")")
  )

results_table_adjpinf <- results_table_adjpinf[grepl("^prior_infection", results_table_adjpinf$term), ]


# Select and rename necessary columns
results_table_adjpinf <- results_table_adjpinf %>%
  select(model, term, hr, lower_ci, upper_ci, conf_int, p_value) %>%
  rename(
    Model = model,
    `Prior Infection` = term,
    `Hazard Ratio` = hr,
    `Lower 95% CI` = lower_ci,
    `Upper 95% CI` = upper_ci,
    `95% CI` = conf_int,
    `p-value` = p_value
  ) %>%
  mutate(HR = paste0(round(`Hazard Ratio`, 2), "(", round(`Lower 95% CI`, 2), "-", round(`Upper 95% CI`, 2), ")")) %>%
  select(`Model`, `Prior Infection`,`HR`, `p-value`) %>%
  pivot_wider(
    names_from = Model,
    values_from = c(`HR`, `p-value`))


print(results_table_adjpinf)


# Period ------------------------------------------------------------------


#Produce null model
#test prior inf cat
periodnull_model <- glm(symp ~ prior_infection + age_cat, data = symp_data[symp_data$symp!=2,])
#Crude model
period_model <- glm(symp ~period, data = symp_data[symp_data$symp!=2,])
# Calculate p-values from likelihood ratio tests
p_val_crude_period <- lrtest(period_model)$Pr[2]
p_val_adjusted <- lrtest(adjust_model, periodnull_model)$Pr[2]

# tidy() from the broom package extracts the necessary information from the models
tidy_period_model <- tidy(period_model, conf.int = TRUE) %>% 
  mutate(model = "Crude Model",
         p_value = p_val_crude_period)

tidy_adjust_model <- tidy(adjust_model, conf.int = TRUE) %>% 
  mutate(model = "Adjusted model",
         p_value = p_val_adjusted)


# Join the results from both models
results_table_adjperiod <- rbind(tidy_period_model, tidy_adjust_model)

# Calculate hazard ratios and confidence intervals
results_table_adjperiod <- results_table_adjperiod %>%
  mutate(
    hr = exp(estimate),
    lower_ci = exp(conf.low),
    upper_ci = exp(conf.high),
    conf_int = paste0("(", round(lower_ci, 2), "-", round(upper_ci, 2), ")")
  )

results_table_adjperiod <- results_table_adjperiod[grepl("^period", results_table_adjperiod$term), ]


# Select and rename necessary columns
results_table_adjperiod <- results_table_adjperiod %>%
  select(model, term, hr, lower_ci, upper_ci, conf_int, p_value) %>%
  rename(
    Model = model,
    `Period` = term,
    `Hazard Ratio` = hr,
    `Lower 95% CI` = lower_ci,
    `Upper 95% CI` = upper_ci,
    `95% CI` = conf_int,
    `p-value` = p_value
  ) %>%
  mutate(HR = paste0(round(`Hazard Ratio`, 2), "(", round(`Lower 95% CI`, 2), "-", round(`Upper 95% CI`, 2), ")")) %>%
  select(`Model`, `Period`,`HR`, `p-value`) %>%
  pivot_wider(
    names_from = Model,
    values_from = c(`HR`, `p-value`))


print(results_table_adjpinf)


adj_combined_results <- list(results_table_adjage, results_table_adjpinf, results_table_adjperiod) %>% 
  lapply(., function(x) setNames(x, c("Variable", "HR (Crude)", "HR (Adjusted)",
                                      "p-value (Crude)", "P-Value (Adjusted)"))) %>% 
  bind_rows()


write.csv(adj_combined_results, "adj_symp_risk.csv")
