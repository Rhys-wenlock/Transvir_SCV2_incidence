require(lmtest)
require(broom)
require(tidyverse)
require(lme4)

#Run 02_cluster_code_final before this
#keep stored in R environment

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
  filter(!is.na(inf)) %>%
  mutate(period = factor(period, levels=c("Pre-Delta", "Delta", "Omicron")))


summary_single_cluster_data <- total_single_cluster_data %>% 
  group_by(cluster_ID) %>%
  summarize(
    num_participants = n(),
    num_positive = sum(ifelse(PCR_result == "Positive", 1, 0), na.rm = TRUE),
    index_id = first(Index_ID),
    index_ct = first(Index_Ct),
    symp = first(symp)
  )

sum(summary_single_cluster_data$num_participants)
summary(summary_single_cluster_data$num_participants)
sum(summary_single_cluster_data$num_positive)
table(summary_single_cluster_data$num_positive)

# Summary Table -----------------------------------------------------------


cols_of_interest <- c("age_cat", "period",
                      "hh_size", "total_child", "Index_Ct", "symp", "sero")
results_list <- list()

for (current_col in cols_of_interest) {
  
  if(is.numeric(total_single_cluster_data[[current_col]])) { 
    # Continuous variable
    current_results_table <- total_single_cluster_data %>%
      summarize(
        variable_level = as.character(current_col),
        num_exposed = paste0(median(!!sym(current_col)), " (", quantile(!!sym(current_col), 0.25), "-", quantile(!!sym(current_col), 0.75), ")"),
        prop_exposed = n() / 468,
        num_inf = sum(inf == 1), 
        prop_inf = (num_inf / n()) * 100, 
        OR = NA, 
        p_value = NA
      )
  } else {
    # Categorical variable
    current_results_table <- total_single_cluster_data %>%
      group_by(!!sym(current_col)) %>%
      summarize(
        num_exposed = as.character(n()),
        prop_exposed = n() / 468,
        num_inf = sum(inf == 1), 
        prop_inf = (num_inf / n()) * 100, 
        OR = NA, 
        p_value = NA
      ) %>%
      mutate(variable_level = as.character(!!sym(current_col)))
  }
  
  # Add the identifier column
  current_results_table$variable <- current_col
  
  # Append the table to the list
  results_list[[current_col]] <- current_results_table
}

# To view results for a specific column, for example "age_cat", use:
# results_list[["age_cat"]]



null_model <- glmer(inf ~1 + (1|Household_ID), family=binomial(), total_single_cluster_data)

# Univariate age ----------------------------------------------------------
#age_cat univariate
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




# Univariate Period -------------------------------------------------------
period_model <- glmer(inf ~ period + (1|Household_ID), family=binomial(), total_single_cluster_data)
period_fef <- exp(fixef(period_model))
period_cint <- exp(confint(period_model))

results_list$period$OR[2] <- 
  paste0(round(period_fef[2],2), " (", round(period_cint[3],2), "-", round(period_cint[7],2), ")")
results_list$period$OR[3] <- 
  paste0(round(period_fef[3],2), " (", round(period_cint[4],2), "-", round(period_cint[8],2), ")")

p_val <- lrtest(period_model, null_model)
results_list$period$p_value <- p_val$`Pr(>Chisq)`[2]


# Univariate hh_size ------------------------------------------------------
hh_size_model <- glmer(inf ~ hh_size + (1|Household_ID), family=binomial(), total_single_cluster_data)
hh_size_fef <- exp(fixef(hh_size_model))
hh_size_cint <- exp(confint(hh_size_model))

results_list$hh_size$OR[1] <- 
  paste0(round(hh_size_fef[2],2), " (", round(hh_size_cint[3],2), "-", round(hh_size_cint[6],2), ")")

p_val <- lrtest(hh_size_model, null_model)
results_list$hh_size$p_value <- p_val$`Pr(>Chisq)`[2]


# Univariate HH Children --------------------------------------------------
hh_child_model <- glmer(inf ~ total_child + (1|Household_ID), family=binomial(), total_single_cluster_data)
hh_child_fef <- exp(fixef(hh_child_model))

se <- sqrt(diag(vcov(hh_child_model)))
hh_child_cint <- data.frame(
  Estimate = hh_child_fef,
  `2.5 %` = hh_child_fef - 1.96*se,
  `97.5 %` = hh_child_fef + 1.96*se
)


results_list$total_child$OR[1] <- 
  paste0(round(hh_child_fef[2],2), " (", round(hh_child_cint[2,2],2), "-", round(hh_child_cint[2,3],2), ")")

p_val <- lrtest(hh_child_model, null_model)
results_list$total_child$p_value <- p_val$`Pr(>Chisq)`[2]



# Univariate Ct -----------------------------------------------------------

Ct_model <- glmer(inf ~ Index_Ct + (1|Household_ID), family=binomial(), total_single_cluster_data)
Ct_fef <- exp(fixef(Ct_model))
Ct_cint <- exp(confint(Ct_model))

results_list$Index_Ct$OR[1] <- 
  paste0(round(Ct_fef[2],2), " (", round(Ct_cint[3],2), "-", round(Ct_cint[6],2), ")")

p_val <- lrtest(Ct_model, null_model)
results_list$Index_Ct$p_value <- p_val$`Pr(>Chisq)`[2]


# Univariate symptom ------------------------------------------------------
symp_model <- glmer(inf ~ symp + (1|Household_ID), family=binomial(), total_single_cluster_data)
symp_fef <- exp(fixef(symp_model))
symp_cint <- exp(confint(symp_model))

results_list$symp$OR[2] <- 
  paste0(round(symp_fef[2],2), " (", round(symp_cint[3],2), "-", round(symp_cint[7],2), ")")
results_list$symp$OR[3] <- 
  paste0(round(symp_fef[3],2), " (", round(symp_cint[4],2), "-", round(symp_cint[8],2), ")")

p_val <- lrtest(symp_model, null_model)
results_list$symp$p_value <- p_val$`Pr(>Chisq)`[2]


# Univariate sero ---------------------------------------------------------
sero_model <- glmer(inf ~ sero + (1|Household_ID), family=binomial(), total_single_cluster_data)
sero_fef <- exp(fixef(sero_model))
sero_cint <- exp(confint(sero_model))

results_list$sero$OR[2] <- 
  paste0(round(sero_fef[2],2), " (", round(sero_cint[3],2), "-", round(sero_cint[6],2), ")")

p_val <- lrtest(sero_model, null_model)
results_list$sero$p_value <- p_val$`Pr(>Chisq)`[2]

results_table <- bind_rows(results_list)
write.csv(results_table, "HCIR_analysis.csv")


# Multivariate Modelling --------------------------------------------------

adjust_model <- glmer(inf ~ sero +sero_contact+ (1|Household_ID), family=binomial(), total_single_cluster_data)
adjust_fef <- exp(fixef(adjust_model))
adjust_cint <- exp(confint(adjust_model))

test_model <- glmer(inf ~ sero+ (1|Household_ID), family=binomial(), total_single_cluster_data[!is.na(total_single_cluster_data$sero_contact),])
lrtest(adjust_model, test_model)  
  
  