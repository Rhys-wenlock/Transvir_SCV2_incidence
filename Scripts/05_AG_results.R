require(dplyr)
require(survival)
require(epiR)
require(lmtest)
require(broom)


#import lexis_inf_split from Lexis.code.R
lexis_inf_split <- read.csv(here::here("Output", "lexis_inf_split.csv"))

#set infection to be binary numerical
lexis_inf_split$infection <- ifelse(lexis_inf_split$lex.Xst==TRUE, 1,0)



# Age (Crude) -------------------------------------------------------------


# Create crude age_model
lexis_inf_split$age_cat <- factor(lexis_inf_split$age_cat, levels=c("18-49", "<5", "5-17", ">50"))
age_model <- coxph(Surv(lex.dur,lex.Xst)~age_cat+cluster(ID,Household_ID), data = lexis_inf_split)
summary(age_model)

# Calculate hazard ratios
model_results_age <- data.frame(
  age_cat = row.names(coef(summary(age_model))),
  hr = exp(coef(age_model)),
  lower_ci_hr = exp(confint(age_model)[, 1]),
  upper_ci_hr = exp(confint(age_model)[, 2])
)

# Clean up the age_cat variable in model_results to match with results_table
model_results_age$age_cat <- str_replace(model_results_age$age_cat, "age_cat", "")

# Create the results table
results_table_age <- lexis_inf_split %>%
  group_by(age_cat) %>%
  summarize(
    num_participants = n_distinct(ID),
    prop_participants = num_participants / 338,
    follow_up_time = sum(lex.dur),
    num_infections = sum(lex.Xst),
    incidence_rate = num_infections / (follow_up_time), # per 1,000 person-time
    lower_ci_incidence_rate = incidence_rate / exp(1.96 * sqrt(1 / num_infections)),
    upper_ci_incidence_rate = incidence_rate * exp(1.96 * sqrt(1 / num_infections))) %>%
  left_join(model_results_age, by = "age_cat") %>%
  mutate(
    incidence_rate_ci = paste0(round(incidence_rate, 2), " (", round(lower_ci_incidence_rate, 2), "-", round(upper_ci_incidence_rate, 2), ")")
  ) %>%
  mutate(HR = paste0(round(hr, 2), "(", round(lower_ci_hr, 2), "-", round(upper_ci_hr, 2), ")")) %>%
  select(age_cat, num_participants, prop_participants, num_infections, incidence_rate_ci, HR)
results_table_age$type <- "Age"
# Print the results table
print(results_table_age)


# Prior Infection (Crude) -------------------------------------------------

### Create crude age_model###
lexis_inf_split$prior_infection[lexis_inf_split$prior_infection==3] <- 2
lexis_inf_split$prior_infection[lexis_inf_split$prior_infection==4] <- 2
lexis_inf_split$prior_infection <- factor(lexis_inf_split$prior_infection, levels=c("0", "1", "2"))

priorinf_model <- coxph(Surv(lex.dur,lex.Xst)~prior_infection+cluster(ID,Household_ID), data = lexis_inf_split)
summary(priorinf_model)


# Calculate hazard ratios
model_results_pi <- data.frame(
  prior_infection = row.names(coef(summary(priorinf_model))),
  hr = exp(coef(priorinf_model)),
  lower_ci_hr = exp(confint(priorinf_model)[, 1]),
  upper_ci_hr = exp(confint(priorinf_model)[, 2])
)

# Clean up the age_cat variable in model_results to match with results_table
model_results_pi$prior_infection <- str_replace(model_results_pi$prior_infection, "prior_infection", "")

# Create the results table
results_table_pi <- lexis_inf_split %>%
  group_by(prior_infection) %>%
  summarize(
    num_participants = n_distinct(ID),
    prop_participants = num_participants / 338,
    follow_up_time = sum(lex.dur),
    num_infections = sum(lex.Xst),
    incidence_rate = num_infections / (follow_up_time), # per 1,000 person-time
    lower_ci_incidence_rate = incidence_rate / exp(1.96 * sqrt(1 / num_infections)),
    upper_ci_incidence_rate = incidence_rate * exp(1.96 * sqrt(1 / num_infections))) %>%
  left_join(model_results_pi, by = "prior_infection") %>%
  mutate(
    incidence_rate_ci = paste0(round(incidence_rate, 2), " (", round(lower_ci_incidence_rate, 2), "-", round(upper_ci_incidence_rate, 2), ")")
  ) %>%
  mutate(HR = paste0(round(hr, 2), "(", round(lower_ci_hr, 2), "-", round(upper_ci_hr, 2), ")")) %>%
  select(prior_infection, num_participants, prop_participants, num_infections, incidence_rate_ci, HR)
results_table_pi$type <- "PI"
# Print the results table
print(results_table_pi)




# HH Size (Crude) ---------------------------------------------------------

lexis_inf_split$hh_size <- factor(lexis_inf_split$hh_size, levels=c("[5,7)", "[7,10)", "[10,Inf)"))
hh_sizemodel <- coxph(Surv(lex.dur,lex.Xst)~hh_size+cluster(ID,Household_ID), data = lexis_inf_split)
summary(hh_sizemodel)

# Calculate hazard ratios
model_results_hh<- data.frame(
  hh_size = row.names(coef(summary(hh_sizemodel))),
  hr = exp(coef(hh_sizemodel)),
  lower_ci_hr = exp(confint(hh_sizemodel)[, 1]),
  upper_ci_hr = exp(confint(hh_sizemodel)[, 2])
)

# Clean up the age_cat variable in model_results to match with results_table
model_results_hh$hh_size <- str_replace(model_results_hh$hh_size, "hh_size", "")

# Create the results table
results_table_hh <- lexis_inf_split %>%
  group_by(hh_size) %>%
  summarize(
    num_participants = n_distinct(ID),
    prop_participants = num_participants / 338,
    follow_up_time = sum(lex.dur),
    num_infections = sum(lex.Xst),
    incidence_rate = num_infections / (follow_up_time), # per 1,000 person-time
    lower_ci_incidence_rate = incidence_rate / exp(1.96 * sqrt(1 / num_infections)),
    upper_ci_incidence_rate = incidence_rate * exp(1.96 * sqrt(1 / num_infections))) %>%
  left_join(model_results_hh, by = "hh_size") %>%
  mutate(
    incidence_rate_ci = paste0(round(incidence_rate, 2), " (", round(lower_ci_incidence_rate, 2), "-", round(upper_ci_incidence_rate, 2), ")")
  ) %>%
  mutate(HR = paste0(round(hr, 2), "(", round(lower_ci_hr, 2), "-", round(upper_ci_hr, 2), ")")) %>%
  select(hh_size, num_participants, prop_participants, num_infections, incidence_rate_ci, HR)
results_table_hh$type <- "HH"
# Print the results table
print(results_table_hh)




# #HIV Crude --------------------------------------------------------------
lexis_inf_split$HIV <- as.factor(lexis_inf_split$HIV)
HIV_model <- coxph(Surv(lex.dur,lex.Xst)~HIV+cluster(ID,Household_ID), data = lexis_inf_split)
summary(HIV_model)

# Calculate hazard ratios
model_results_HIV<- data.frame(
  HIV = row.names(coef(summary(HIV_model))),
  hr = exp(coef(HIV_model)),
  lower_ci_hr = exp(confint(HIV_model)[, 1]),
  upper_ci_hr = exp(confint(HIV_model)[, 2])
)

# Clean up the age_cat variable in model_results to match with results_table
model_results_HIV$HIV <- str_replace(model_results_HIV$HIV, "HIV", "")

# Create the results table
results_table_HIV <- lexis_inf_split %>%
  group_by(HIV) %>%
  summarize(
    num_participants = n_distinct(ID),
    prop_participants = num_participants / 338,
    follow_up_time = sum(lex.dur),
    num_infections = sum(lex.Xst),
    incidence_rate = num_infections / (follow_up_time), # per 1,000 person-time
    lower_ci_incidence_rate = incidence_rate / exp(1.96 * sqrt(1 / num_infections)),
    upper_ci_incidence_rate = incidence_rate * exp(1.96 * sqrt(1 / num_infections))) %>%
  left_join(model_results_HIV, by = "HIV") %>%
  mutate(
    incidence_rate_ci = paste0(round(incidence_rate, 2), " (", round(lower_ci_incidence_rate, 2), "-", round(upper_ci_incidence_rate, 2), ")")
  ) %>%
  mutate(HR = paste0(round(hr, 2), "(", round(lower_ci_hr, 2), "-", round(upper_ci_hr, 2), ")")) %>%
  select(HIV, num_participants, prop_participants, num_infections, incidence_rate_ci, HR)
results_table_HIV$type <- "HIV"
# Print the results table
print(results_table_HIV)


# #Employed Crude ---------------------------------------------------------

lexis_inf_split$Employed <- as.factor(lexis_inf_split$Employed)
Employed_model <- coxph(Surv(lex.dur,lex.Xst)~Employed+cluster(ID,Household_ID), data = lexis_inf_split)
summary(Employed_model)

# Calculate hazard ratios
model_results_Employed<- data.frame(
  Employed = row.names(coef(summary(Employed_model))),
  hr = exp(coef(Employed_model)),
  lower_ci_hr = exp(confint(Employed_model)[, 1]),
  upper_ci_hr = exp(confint(Employed_model)[, 2])
)

# Clean up the age_cat variable in model_results to match with results_table
model_results_Employed$Employed <- str_replace(model_results_Employed$Employed, "Employed", "")

# Create the results table
results_table_Employed <- lexis_inf_split %>%
  group_by(Employed) %>%
  summarize(
    num_participants = n_distinct(ID),
    prop_participants = num_participants / 338,
    follow_up_time = sum(lex.dur),
    num_infections = sum(lex.Xst),
    incidence_rate = num_infections / (follow_up_time), # per 1,000 person-time
    lower_ci_incidence_rate = incidence_rate / exp(1.96 * sqrt(1 / num_infections)),
    upper_ci_incidence_rate = incidence_rate * exp(1.96 * sqrt(1 / num_infections))) %>%
  left_join(model_results_Employed, by = "Employed") %>%
  mutate(
    incidence_rate_ci = paste0(round(incidence_rate, 2), " (", round(lower_ci_incidence_rate, 2), "-", round(upper_ci_incidence_rate, 2), ")")
  ) %>%
  mutate(HR = paste0(round(hr, 2), "(", round(lower_ci_hr, 2), "-", round(upper_ci_hr, 2), ")")) %>%
  select(Employed, num_participants, prop_participants, num_infections, incidence_rate_ci, HR)
results_table_Employed$type <- "Employed"
# Print the results table
print(results_table_Employed)


# #Steroid Crude ----------------------------------------------------------
lexis_inf_split$Steroid <- as.factor(lexis_inf_split$Steroid)
Steroid_model <- coxph(Surv(lex.dur,lex.Xst)~Steroid+cluster(ID,Household_ID), data = lexis_inf_split)
summary(Steroid_model)

# Calculate hazard ratios
model_results_Steroid<- data.frame(
  Steroid = row.names(coef(summary(Steroid_model))),
  hr = exp(coef(Steroid_model)),
  lower_ci_hr = exp(confint(Steroid_model)[, 1]),
  upper_ci_hr = exp(confint(Steroid_model)[, 2])
)

# Clean up the age_cat variable in model_results to match with results_table
model_results_Steroid$Steroid <- str_replace(model_results_Steroid$Steroid, "Steroid", "")

# Create the results table
results_table_Steroid <- lexis_inf_split %>%
  group_by(Steroid) %>%
  summarize(
    num_participants = n_distinct(ID),
    prop_participants = num_participants / 338,
    follow_up_time = sum(lex.dur),
    num_infections = sum(lex.Xst),
    incidence_rate = num_infections / (follow_up_time), # per 1,000 person-time
    lower_ci_incidence_rate = incidence_rate / exp(1.96 * sqrt(1 / num_infections)),
    upper_ci_incidence_rate = incidence_rate * exp(1.96 * sqrt(1 / num_infections))) %>%
  left_join(model_results_Steroid, by = "Steroid") %>%
  mutate(
    incidence_rate_ci = paste0(round(incidence_rate, 2), " (", round(lower_ci_incidence_rate, 2), "-", round(upper_ci_incidence_rate, 2), ")")
  ) %>%
  mutate(HR = paste0(round(hr, 2), "(", round(lower_ci_hr, 2), "-", round(upper_ci_hr, 2), ")")) %>%
  select(Steroid, num_participants, prop_participants, num_infections, incidence_rate_ci, HR)
results_table_Steroid$type <- "Steroid"
# Print the results table
print(results_table_Steroid)


# #Cancer crude -----------------------------------------------------------
lexis_inf_split$Cancer <- as.factor(lexis_inf_split$Cancer)
Cancer_model <- coxph(Surv(lex.dur,lex.Xst)~Cancer+cluster(ID,Household_ID), data = lexis_inf_split)
summary(Cancer_model)

# Calculate hazard ratios
model_results_Cancer<- data.frame(
  Cancer = row.names(coef(summary(Cancer_model))),
  hr = exp(coef(Cancer_model)),
  lower_ci_hr = exp(confint(Cancer_model)[, 1]),
  upper_ci_hr = exp(confint(Cancer_model)[, 2])
)

# Clean up the age_cat variable in model_results to match with results_table
model_results_Cancer$Cancer <- str_replace(model_results_Cancer$Cancer, "Cancer", "")

# Create the results table
results_table_Cancer <- lexis_inf_split %>%
  group_by(Cancer) %>%
  summarize(
    num_participants = n_distinct(ID),
    prop_participants = num_participants / 338,
    follow_up_time = sum(lex.dur),
    num_infections = sum(lex.Xst),
    incidence_rate = num_infections / (follow_up_time), # per 1,000 person-time
    lower_ci_incidence_rate = incidence_rate / exp(1.96 * sqrt(1 / num_infections)),
    upper_ci_incidence_rate = incidence_rate * exp(1.96 * sqrt(1 / num_infections))) %>%
  left_join(model_results_Cancer, by = "Cancer") %>%
  mutate(
    incidence_rate_ci = paste0(round(incidence_rate, 2), " (", round(lower_ci_incidence_rate, 2), "-", round(upper_ci_incidence_rate, 2), ")")
  ) %>%
  mutate(HR = paste0(round(hr, 2), "(", round(lower_ci_hr, 2), "-", round(upper_ci_hr, 2), ")")) %>%
  select(Cancer, num_participants, prop_participants, num_infections, incidence_rate_ci, HR)
results_table_Cancer$type <- "Cancer"
# Print the results table
print(results_table_Cancer)


# #Diabetes Crude ---------------------------------------------------------

lexis_inf_split$Diabetes <- as.factor(lexis_inf_split$Diabetes)
Diabetes_model <- coxph(Surv(lex.dur,lex.Xst)~Diabetes+cluster(ID,Household_ID), data = lexis_inf_split)
summary(Diabetes_model)

# Calculate hazard ratios
model_results_Diabetes<- data.frame(
  Diabetes = row.names(coef(summary(Diabetes_model))),
  hr = exp(coef(Diabetes_model)),
  lower_ci_hr = exp(confint(Diabetes_model)[, 1]),
  upper_ci_hr = exp(confint(Diabetes_model)[, 2])
)

# Clean up the age_cat variable in model_results to match with results_table
model_results_Diabetes$Diabetes <- str_replace(model_results_Diabetes$Diabetes, "Diabetes", "")

# Create the results table
results_table_Diabetes <- lexis_inf_split %>%
  group_by(Diabetes) %>%
  summarize(
    num_participants = n_distinct(ID),
    prop_participants = num_participants / 338,
    follow_up_time = sum(lex.dur),
    num_infections = sum(lex.Xst),
    incidence_rate = num_infections / (follow_up_time), # per 1,000 person-time
    lower_ci_incidence_rate = incidence_rate / exp(1.96 * sqrt(1 / num_infections)),
    upper_ci_incidence_rate = incidence_rate * exp(1.96 * sqrt(1 / num_infections))) %>%
  left_join(model_results_Diabetes, by = "Diabetes") %>%
  mutate(
    incidence_rate_ci = paste0(round(incidence_rate, 2), " (", round(lower_ci_incidence_rate, 2), "-", round(upper_ci_incidence_rate, 2), ")")
  ) %>%
  mutate(HR = paste0(round(hr, 2), "(", round(lower_ci_hr, 2), "-", round(upper_ci_hr, 2), ")")) %>%
  select(Diabetes, num_participants, prop_participants, num_infections, incidence_rate_ci, HR)
results_table_Diabetes$type <- "Diabetes"
# Print the results table
print(results_table_Diabetes)


# #HTN Crude -------------------------------------------------------------

lexis_inf_split$HTN <- as.factor(lexis_inf_split$HTN)
HTN_model <- coxph(Surv(lex.dur,lex.Xst)~HTN+cluster(ID,Household_ID), data = lexis_inf_split)
summary(HTN_model)

# Calculate hazard ratios
model_results_HTN<- data.frame(
  HTN = row.names(coef(summary(HTN_model))),
  hr = exp(coef(HTN_model)),
  lower_ci_hr = exp(confint(HTN_model)[, 1]),
  upper_ci_hr = exp(confint(HTN_model)[, 2])
)

# Clean up the age_cat variable in model_results to match with results_table
model_results_HTN$HTN <- str_replace(model_results_HTN$HTN, "HTN", "")

# Create the results table
results_table_HTN <- lexis_inf_split %>%
  group_by(HTN) %>%
  summarize(
    num_participants = n_distinct(ID),
    prop_participants = num_participants / 338,
    follow_up_time = sum(lex.dur),
    num_infections = sum(lex.Xst),
    incidence_rate = num_infections / (follow_up_time), # per 1,000 person-time
    lower_ci_incidence_rate = incidence_rate / exp(1.96 * sqrt(1 / num_infections)),
    upper_ci_incidence_rate = incidence_rate * exp(1.96 * sqrt(1 / num_infections))) %>%
  left_join(model_results_HTN, by = "HTN") %>%
  mutate(
    incidence_rate_ci = paste0(round(incidence_rate, 2), " (", round(lower_ci_incidence_rate, 2), "-", round(upper_ci_incidence_rate, 2), ")")
  ) %>%
  mutate(HR = paste0(round(hr, 2), "(", round(lower_ci_hr, 2), "-", round(upper_ci_hr, 2), ")")) %>%
  select(HTN, num_participants, prop_participants, num_infections, incidence_rate_ci, HR)
results_table_HTN$type <- "HTN"
# Print the results table
print(results_table_HTN)


# #Smoking Crude ---------------------------------------------------------

lexis_inf_split$Smoking <- as.factor(lexis_inf_split$Smoking)
Smoking_model <- coxph(Surv(lex.dur,lex.Xst)~Smoking+cluster(ID,Household_ID), data = lexis_inf_split)
summary(Smoking_model)

# Calculate hazard ratios
model_results_Smoking<- data.frame(
  Smoking = row.names(coef(summary(Smoking_model))),
  hr = exp(coef(Smoking_model)),
  lower_ci_hr = exp(confint(Smoking_model)[, 1]),
  upper_ci_hr = exp(confint(Smoking_model)[, 2])
)

# Clean up the age_cat variable in model_results to match with results_table
model_results_Smoking$Smoking <- str_replace(model_results_Smoking$Smoking, "Smoking", "")

# Create the results table
results_table_Smoking <- lexis_inf_split %>%
  group_by(Smoking) %>%
  summarize(
    num_participants = n_distinct(ID),
    prop_participants = num_participants / 338,
    follow_up_time = sum(lex.dur),
    num_infections = sum(lex.Xst),
    incidence_rate = num_infections / (follow_up_time), # per 1,000 person-time
    lower_ci_incidence_rate = incidence_rate / exp(1.96 * sqrt(1 / num_infections)),
    upper_ci_incidence_rate = incidence_rate * exp(1.96 * sqrt(1 / num_infections))) %>%
  left_join(model_results_Smoking, by = "Smoking") %>%
  mutate(
    incidence_rate_ci = paste0(round(incidence_rate, 2), " (", round(lower_ci_incidence_rate, 2), "-", round(upper_ci_incidence_rate, 2), ")")
  ) %>%
  mutate(HR = paste0(round(hr, 2), "(", round(lower_ci_hr, 2), "-", round(upper_ci_hr, 2), ")")) %>%
  select(Smoking, num_participants, prop_participants, num_infections, incidence_rate_ci, HR)
results_table_Smoking$type <- "Smoking"
# Print the results table
print(results_table_Smoking)


# Vaccinated Crude --------------------------------------------------------
lexis_inf_split$vaccinated <- as.factor(lexis_inf_split$vaccinated)
vaccinated_model <- coxph(Surv(lex.dur,lex.Xst)~vaccinated+cluster(ID,Household_ID), data = lexis_inf_split)
summary(vaccinated_model)

# Calculate hazard ratios
model_results_vaccinated<- data.frame(
  vaccinated = row.names(coef(summary(vaccinated_model))),
  hr = exp(coef(vaccinated_model)),
  lower_ci_hr = exp(confint(vaccinated_model)[, 1]),
  upper_ci_hr = exp(confint(vaccinated_model)[, 2])
)

# Clean up the age_cat variable in model_results to match with results_table
model_results_vaccinated$vaccinated <- str_replace(model_results_vaccinated$vaccinated, "vaccinated", "")

# Create the results table
results_table_vaccinated <- lexis_inf_split %>%
  group_by(vaccinated) %>%
  summarize(
    num_participants = n_distinct(ID),
    prop_participants = num_participants / 338,
    follow_up_time = sum(lex.dur),
    num_infections = sum(lex.Xst),
    incidence_rate = num_infections / (follow_up_time), # per 1,000 person-time
    lower_ci_incidence_rate = incidence_rate / exp(1.96 * sqrt(1 / num_infections)),
    upper_ci_incidence_rate = incidence_rate * exp(1.96 * sqrt(1 / num_infections))) %>%
  left_join(model_results_vaccinated, by = "vaccinated") %>%
  mutate(
    incidence_rate_ci = paste0(round(incidence_rate, 2), " (", round(lower_ci_incidence_rate, 2), "-", round(upper_ci_incidence_rate, 2), ")")
  ) %>%
  mutate(HR = paste0(round(hr, 2), "(", round(lower_ci_hr, 2), "-", round(upper_ci_hr, 2), ")")) %>%
  select(vaccinated, num_participants, prop_participants, num_infections, incidence_rate_ci, HR)
results_table_vaccinated$type <- "vaccinated"
# Print the results table
print(results_table_vaccinated)




df.combined <- list(results_table_age, results_table_pi, results_table_hh,
                    results_table_Cancer, results_table_Diabetes, results_table_Employed, 
                    results_table_HIV, results_table_HTN, results_table_Smoking, results_table_vaccinated) %>% 
  lapply(., function(x) setNames(x, c("Variable", "num_participants", "prop_participants",
                                      "num_infections", "incidence_rate_ci", "HR","type"))) %>% 
  bind_rows()

write.csv(df.combined, "crude_AGresults.csv")


# Adjusted Models ---------------------------------------------------------

#Produce adjusted models
lexis_inf_split$age_cat <- factor(lexis_inf_split$age_cat, levels=c("18-49", "<5", "5-17", ">50"))
lexis_inf_split$hh_size <- factor(lexis_inf_split$hh_size, levels=c("[5,7)", "[7,10)", "[10,Inf)"))
lexis_inf_split$prior_infection <- factor(lexis_inf_split$prior_infection, levels=c("0", "1", "2"))
adjust_model <- coxph(Surv(lex.dur,lex.Xst)~age_cat+prior_infection+vaccinated+hh_size+ strata(period)+cluster(ID,Household_ID), data = lexis_inf_split)
summary(adjust_model)


#test hh size
priorinf_age_model <- coxph(Surv(lex.dur,lex.Xst)~prior_infection+vaccinated+age_cat+strata(period)+cluster(ID,Household_ID), data = lexis_inf_split)
lrtest(adjust_model, priorinf_age_model)

adjusted_results <- data.frame(
  age_cat = row.names(coef(summary(adjust_model))),
  hr = exp(coef(adjust_model)),
  lower_ci_hr = exp(confint(adjust_model)[, 1]),
  upper_ci_hr = exp(confint(adjust_model)[, 2])
)

# Age (adjusted) ----------------------------------------------------------

#Produce null model
#test age cat
agenull_model <- coxph(Surv(lex.dur,lex.Xst)~prior_infection+vaccinated+hh_size+strata(period)+cluster(ID,Household_ID), data = lexis_inf_split)

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



# Prior Infection (Adjusted) ----------------------------------------------

#Produce null model
#test age cat
priorinfnull_model <- coxph(Surv(lex.dur,lex.Xst)~age_cat+vaccinated+hh_size+strata(period)+cluster(ID,Household_ID), data = lexis_inf_split)

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



# Household Size (Adjusted) -----------------------------------------------

#Produce null model
#test household size
hhsizenull_model <- coxph(Surv(lex.dur,lex.Xst)~age_cat+vaccinated+prior_infection+strata(period)+cluster(ID,Household_ID), data = lexis_inf_split)

# Calculate p-values from likelihood ratio tests
p_val_crude_hhsize <- lrtest(hh_sizemodel)$Pr[2]
p_val_adjusted <- lrtest(adjust_model, hhsizenull_model)$Pr[2]

# tidy() from the broom package extracts the necessary information from the models
tidy_hhsize_model <- tidy(hh_sizemodel, conf.int = TRUE) %>% 
  mutate(model = "Crude Model",
         p_value = p_val_crude_hhsize)

tidy_adjust_model <- tidy(adjust_model, conf.int = TRUE) %>% 
  mutate(model = "Adjusted model",
         p_value = p_val_adjusted)


# Join the results from both models
results_table_adjhhsize <- rbind(tidy_hhsize_model, tidy_adjust_model)

# Calculate hazard ratios and confidence intervals
results_table_adjhhsize <- results_table_adjhhsize %>%
  mutate(
    hr = exp(estimate),
    lower_ci = exp(conf.low),
    upper_ci = exp(conf.high),
    conf_int = paste0("(", round(lower_ci, 2), "-", round(upper_ci, 2), ")")
  )

results_table_adjhhsize <- results_table_adjhhsize[grepl("^hh_size", results_table_adjhhsize$term), ]


# Select and rename necessary columns
results_table_adjhhsize <- results_table_adjhhsize %>%
  select(model, term, hr, lower_ci, upper_ci, conf_int, p_value) %>%
  rename(
    Model = model,
    `Household Size` = term,
    `Hazard Ratio` = hr,
    `Lower 95% CI` = lower_ci,
    `Upper 95% CI` = upper_ci,
    `95% CI` = conf_int,
    `p-value` = p_value
  ) %>%
  mutate(HR = paste0(round(`Hazard Ratio`, 2), "(", round(`Lower 95% CI`, 2), "-", round(`Upper 95% CI`, 2), ")")) %>%
  select(`Model`, `Household Size`,`HR`, `p-value`) %>%
  pivot_wider(
    names_from = Model,
    values_from = c(`HR`, `p-value`))


print(results_table_adjhhsize)

adj_combined_results <- list(results_table_adjage, results_table_adjpinf, results_table_adjhhsize) %>% 
  lapply(., function(x) setNames(x, c("Variable", "HR (Crude)", "HR (Adjusted)",
                                      "p-value (Crude)", "P-Value (Adjusted)"))) %>% 
  bind_rows()

write.csv(adj_combined_results, "adj_AG_combined_results.csv")


# Stratified Results by Delta/Omicron -------------------------------------

#Age

adjust_model <- coxph(Surv(lex.dur,lex.Xst)~age_cat+prior_infection+vaccinated+hh_size+ period+cluster(ID,Household_ID), data = lexis_inf_split)
age_interactmodel <- coxph(Surv(lex.dur,lex.Xst)~prior_infection+vaccinated+hh_size+age_cat*period+cluster(ID,Household_ID), data = lexis_inf_split)


# Calculate p-values from likelihood ratio tests
p_val_interact_age <- lrtest(adjust_model, age_interactmodel)$Pr[2]


# tidy() from the broom package extracts the necessary information from the models
lexis_inf_split$period <- factor(lexis_inf_split$period, levels=c("Delta", "Omicron"))
age_interactmodel <- coxph(Surv(lex.dur,lex.Xst)~prior_infection+vaccinated+hh_size+age_cat*period+cluster(ID,Household_ID), data = lexis_inf_split)

tidy_age_model_delta <- tidy(age_interactmodel, conf.int = TRUE) %>% 
  mutate(model = "Delta", p_value = NA)

lexis_inf_split$period <- factor(lexis_inf_split$period, levels=c("Omicron", "Delta"))
age_interactmodel <- coxph(Surv(lex.dur,lex.Xst)~prior_infection+vaccinated+hh_size+age_cat*period+cluster(ID,Household_ID), data = lexis_inf_split)

tidy_age_model_omicron <- tidy(age_interactmodel, conf.int = TRUE) %>% 
  mutate(model = "Omicron", 
         p_value = p_val_interact_age)

results_table_interact_age <- rbind(tidy_age_model_delta, tidy_age_model_omicron)

results_table_interact_age <- results_table_interact_age[grepl("^age_cat", results_table_interact_age$term), ]

# Calculate hazard ratios and confidence intervals
results_table_interact_age <- results_table_interact_age %>%
  mutate(
    hr = exp(estimate),
    lower_ci = exp(conf.low),
    upper_ci = exp(conf.high),
    conf_int = paste0("(", round(lower_ci, 2), "-", round(upper_ci, 2), ")")
  )

# Select and rename necessary columns
results_table_interact_age <- results_table_interact_age %>%
  select(model, term, hr, lower_ci, upper_ci, conf_int, p_value) %>%
  rename(
    Model = model,
    `Variable` = term,
    `Hazard Ratio` = hr,
    `Lower 95% CI` = lower_ci,
    `Upper 95% CI` = upper_ci,
    `95% CI` = conf_int,
    `p-value` = p_value
  ) %>%
  mutate(HR = paste0(round(`Hazard Ratio`, 2), "(", round(`Lower 95% CI`, 2), "-", round(`Upper 95% CI`, 2), ")")) %>%
  select(`Model`, `Variable`,`HR`, `p-value`) %>%
  pivot_wider(
    names_from = Model,
    values_from = c(`HR`, `p-value`))

results_table_interact_age <- results_table_interact_age[c(1:3),]

#Prior Infectin

adjust_model <- coxph(Surv(lex.dur,lex.Xst)~age_cat+prior_infection+vaccinated+hh_size+period+cluster(ID,Household_ID), data = lexis_inf_split)
pinf_interactmodel <- coxph(Surv(lex.dur,lex.Xst)~vaccinated+hh_size+age_cat+prior_infection*period+cluster(ID,Household_ID), data = lexis_inf_split)


# Calculate p-values from likelihood ratio tests
p_val_interact_pinf <- lrtest(adjust_model, pinf_interactmodel)$Pr[2]


# tidy() from the broom package extracts the necessary information from the models
lexis_inf_split$period <- factor(lexis_inf_split$period, levels=c("Delta", "Omicron"))
pinf_interactmodel <- coxph(Surv(lex.dur,lex.Xst)~vaccinated+hh_size+age_cat+prior_infection*period+cluster(ID,Household_ID), data = lexis_inf_split)

tidy_pinf_model_delta <- tidy(pinf_interactmodel, conf.int = TRUE) %>% 
  mutate(model = "Delta", p_value = NA)

lexis_inf_split$period <- factor(lexis_inf_split$period, levels=c("Omicron", "Delta"))
pinf_interactmodel <- coxph(Surv(lex.dur,lex.Xst)~vaccinated+hh_size+age_cat+prior_infection*period+cluster(ID,Household_ID), data = lexis_inf_split)

tidy_pinf_model_omicron <- tidy(pinf_interactmodel, conf.int = TRUE) %>% 
  mutate(model = "Omicron", 
         p_value = p_val_interact_pinf)

results_table_interact_pinf <- rbind(tidy_pinf_model_delta, tidy_pinf_model_omicron)

results_table_interact_pinf <- results_table_interact_pinf[grepl("^prior_infection", results_table_interact_pinf$term), ]

# Calculate hazard ratios and confidence intervals
results_table_interact_pinf <- results_table_interact_pinf %>%
  mutate(
    hr = exp(estimate),
    lower_ci = exp(conf.low),
    upper_ci = exp(conf.high),
    conf_int = paste0("(", round(lower_ci, 2), "-", round(upper_ci, 2), ")")
  )

# Select and rename necessary columns
results_table_interact_pinf <- results_table_interact_pinf %>%
  select(model, term, hr, lower_ci, upper_ci, conf_int, p_value) %>%
  rename(
    Model = model,
    `Variable` = term,
    `Hazard Ratio` = hr,
    `Lower 95% CI` = lower_ci,
    `Upper 95% CI` = upper_ci,
    `95% CI` = conf_int,
    `p-value` = p_value
  ) %>%
  mutate(HR = paste0(round(`Hazard Ratio`, 2), "(", round(`Lower 95% CI`, 2), "-", round(`Upper 95% CI`, 2), ")")) %>%
  select(`Model`, `Variable`,`HR`, `p-value`) %>%
  pivot_wider(
    names_from = Model,
    values_from = c(`HR`, `p-value`))

results_table_interact_pinf <- results_table_interact_pinf[c(1:2),]

results_interact_combined <- rbind(results_table_interact_age, results_table_interact_pinf)

write.csv(results_interact_combined, "results_AG_interact_combined.csv")


# Supplementary waning analysis -------------------------------------------

#Produce adjusted models
lexis_inf_split$age_cat <- factor(lexis_inf_split$age_cat, levels=c("18-49", "<5", "5-17", ">50"))
lexis_inf_split$hh_size <- factor(lexis_inf_split$hh_size, levels=c("[5,7)", "[7,10)", "[10,Inf)"))
lexis_inf_split$prior_infection <- factor(lexis_inf_split$prior_infection, levels=c("0", "1", "2"))
adjust_model <- coxph(Surv(lex.dur,lex.Xst)~age_cat+combined_variable+vaccinated+hh_size+ strata(period)+cluster(ID,Household_ID), data = lexis_inf_split)
summary(adjust_model)



#Produce null model
#test age cat
combinednull_model <- coxph(Surv(lex.dur,lex.Xst)~age_cat+vaccinated+hh_size+strata(period)+cluster(ID,Household_ID), data = lexis_inf_split)

# Calculate p-values from likelihood ratio tests
p_val_crude_priorinf <- lrtest(priorinf_model)$Pr[2]
p_val_adjusted <- lrtest(adjust_model, combinednull_model)$Pr[2]

# tidy() from the broom package extracts the necessary information from the models

tidy_adjust_model <- tidy(adjust_model, conf.int = TRUE) %>% 
  mutate(model = "Adjusted model",
         p_value = p_val_adjusted)


# Join the results from both models
results_table_adjcomb <- tidy_adjust_model

# Calculate hazard ratios and confidence intervals
results_table_adjcomb <- results_table_adjcomb %>%
  mutate(
    hr = exp(estimate),
    lower_ci = exp(conf.low),
    upper_ci = exp(conf.high),
    conf_int = paste0("(", round(lower_ci, 2), "-", round(upper_ci, 2), ")")
  )

results_table_adjcomb <- results_table_adjcomb[grepl("^combined_variable", results_table_adjcomb$term), ]


# Select and rename necessary columns
results_table_adjcomb <- results_table_adjcomb %>%
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


print(results_table_adjcomb)

#Interaction analysis

adjust_model <- coxph(Surv(lex.dur,lex.Xst)~age_cat+combined_variable+vaccinated+hh_size+ period+cluster(ID,Household_ID), data = lexis_inf_split)
comb_interactmodel <- coxph(Surv(lex.dur,lex.Xst)~vaccinated+hh_size+age_cat+combined_variable*period+cluster(ID,Household_ID), data = lexis_inf_split)


# Calculate p-values from likelihood ratio tests
p_val_interact_comb <- lrtest(adjust_model, comb_interactmodel)$Pr[2]


# tidy() from the broom package extracts the necessary information from the models
lexis_inf_split$period <- factor(lexis_inf_split$period, levels=c("Delta", "Omicron"))
comb_interactmodel <- coxph(Surv(lex.dur,lex.Xst)~vaccinated+hh_size+age_cat+combined_variable*period+cluster(ID,Household_ID), data = lexis_inf_split)

tidy_comb_model_delta <- tidy(comb_interactmodel, conf.int = TRUE) %>% 
  mutate(model = "Delta", p_value = NA)

lexis_inf_split$period <- factor(lexis_inf_split$period, levels=c("Omicron", "Delta"))
comb_interactmodel <- coxph(Surv(lex.dur,lex.Xst)~vaccinated+hh_size+age_cat+combined_variable*period+cluster(ID,Household_ID), data = lexis_inf_split)

tidy_comb_model_omicron <- tidy(comb_interactmodel, conf.int = TRUE) %>% 
  mutate(model = "Omicron", 
         p_value = p_val_interact_pinf)

results_table_interact_comb <- rbind(tidy_comb_model_delta, tidy_comb_model_omicron)

results_table_interact_comb <- results_table_interact_comb[grepl("^combined_variable", results_table_interact_comb$term), ]

# Calculate hazard ratios and confidence intervals
results_table_interact_comb <- results_table_interact_comb %>%
  mutate(
    hr = exp(estimate),
    lower_ci = exp(conf.low),
    upper_ci = exp(conf.high),
    conf_int = paste0("(", round(lower_ci, 2), "-", round(upper_ci, 2), ")")
  )

# Select and rename necessary columns
results_table_interact_comb <- results_table_interact_comb %>%
  select(model, term, hr, lower_ci, upper_ci, conf_int, p_value) %>%
  rename(
    Model = model,
    `Variable` = term,
    `Hazard Ratio` = hr,
    `Lower 95% CI` = lower_ci,
    `Upper 95% CI` = upper_ci,
    `95% CI` = conf_int,
    `p-value` = p_value
  ) %>%
  mutate(HR = paste0(round(`Hazard Ratio`, 2), "(", round(`Lower 95% CI`, 2), "-", round(`Upper 95% CI`, 2), ")")) %>%
  select(`Model`, `Variable`,`HR`, `p-value`) %>%
  pivot_wider(
    names_from = Model,
    values_from = c(`HR`, `p-value`))

results_table_interact_comb <- results_table_interact_comb[c(1:4),]

supp_table_comb <- cbind(results_table_adjcomb, results_table_interact_comb)

write.csv(supp_table_comb, "supp_table_comb.csv")



