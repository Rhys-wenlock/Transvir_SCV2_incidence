#This produces a summary table of the number of infections, number per period, number of 1st/2nd/3rd infections for all sensitivity analysis
#Supplementary Table 2

require(tidyverse)

#Import dataframes
all_inf_data <- read_csv(here::here("Data","all_inf_data_base.csv"))
all_inf_data1 <- read_csv(here::here("Data","all_inf_data_1.csv"))
all_inf_data2 <- read_csv(here::here("Data","all_inf_data_2.csv"))
all_inf_data3 <- read_csv(here::here("Data","all_inf_data_3.csv"))
all_inf_data4 <- read_csv(here::here("Data","all_inf_data_4.csv"))
all_inf_data5 <- read_csv(here::here("Data","all_inf_data_5.csv"))
all_inf_data6 <- read_csv(here::here("Data","all_inf_data_6.csv"))
all_inf_data7 <- read_csv(here::here("Data","all_inf_data_7.csv"))
all_inf_data8 <- read_csv(here::here("Data","all_inf_data_8.csv"))

# Assuming you have 8 dataframes: df1, df2, ..., df8
list_of_dataframes <- list(all_inf_data, all_inf_data1, all_inf_data2,
                           all_inf_data3, all_inf_data4, all_inf_data5,
                           all_inf_data6, all_inf_data7, all_inf_data8)
names(list_of_dataframes) <- c("Base", "Scenario 1 (IQRS1N + NV)", "Scenario 2 (IQR1SN + IQR1NV)",
                               "Scenario 3 (MedianS + MedianNV)", "Scenario 4 (MedianN + MedianNV)", 
                               "Scenario 5 (MixtureS + MixtureNV)", "Scenario 6 (MixtureN + MixtureNV)",
                               "Scenario 7 (90 days, MedianSN + MedianNV)", "Scenario 8 (PCR-positive only)")  # Dataset names

main_results <- data.frame()

for (i in 1:length(list_of_dataframes)) {
  all_inf_data <- list_of_dataframes[[i]]
  
  visit_number <- length(visit_data_fu$Participant_ID[!is.na(visit_data_fu$Date) & !is.na(visit_data_fu$PCR_result)])
  no_pos_visits <- sum(all_inf_data$pcr_num, na.rm = TRUE)
  no_infected <- length(unique(all_inf_data$Participant_ID[all_inf_data$total_inf>=1]))
  prop_infected <- (no_infected/338)*100
  no_infections <- all_inf_data %>%
    group_by(Participant_ID) %>%
    summarize(infections = max(total_inf))
  no_infections <- sum(no_infections$infections)
  no_1st <- length(unique(all_inf_data$Participant_ID[all_inf_data$inf_num==1]))
  no_2nd <- length(unique(all_inf_data$Participant_ID[all_inf_data$inf_num==2]))
  no_3rd <- length(unique(all_inf_data$Participant_ID[all_inf_data$inf_num==3]))
  no_1PCR <- length(unique(all_inf_data$Participant_ID[all_inf_data$pcr_num==1]))
  no_2PCR <- length(unique(all_inf_data$Participant_ID[all_inf_data$pcr_num==2]))
  no_3PCR <- length(unique(all_inf_data$Participant_ID[all_inf_data$pcr_num==3]))
  
  inf_predelta <- length(unique(all_inf_data$Participant_ID[all_inf_data$Inf_date<="2021-07-07"]))
  inf_delta <- length(unique(all_inf_data$Participant_ID[all_inf_data$Inf_date>"2021-07-07" & 
                                                           all_inf_data$Inf_date<="2021-12-04"]))
  inf_omicron <- length(unique(all_inf_data$Participant_ID[all_inf_data$Inf_date>"2021-12-04"]))
  FU_predelta <- length(unique(visit_data_fu$Participant_ID[visit_data$Date<="2021-07-07"]))
  FU_delta <- length(unique(visit_data_fu$Participant_ID[visit_data$Date>"2021-07-07" & 
                                                           visit_data$Date <="2021-12-04"]))
  FU_omicron <- length(unique(visit_data_fu$Participant_ID[visit_data$Date>"2021-12-04"]))
  AR_predelta <- (inf_predelta/FU_predelta)*100
  AR_delta <- (inf_delta/FU_delta)*100
  AR_omicron <- (inf_omicron/FU_omicron)*100
  inf_length <- summary(all_inf_data$length_infection)
  all_inf_data$length_cut <- cut(all_inf_data$length_infection, breaks=c(-Inf,7,14,21,Inf),)
  table(all_inf_data$length_cut, all_inf_data$type)
  pcr_neg <- length(all_inf_data$Participant_ID[all_inf_data$type=="PCR_neg"])
  
  cbind(visit_number, no_pos_visits, no_infected, prop_infected, no_infections, no_1st, no_2nd, no_3rd, 
        no_1PCR, no_2PCR, no_3PCR, inf_predelta, inf_delta, inf_omicron, FU_predelta, 
        FU_delta, FU_omicron, AR_predelta, AR_delta, AR_omicron, pcr_neg)
  
  # Consolidate results into a dataframe
  results_df <- as.data.frame(cbind(visit_number, no_pos_visits, no_infected, prop_infected, 
                                    no_infections, no_1st, no_2nd, no_3rd, 
                                    no_1PCR, no_2PCR, no_3PCR, inf_predelta, inf_delta, inf_omicron, 
                                    FU_predelta, FU_delta, FU_omicron, AR_predelta, AR_delta, 
                                    AR_omicron, pcr_neg))
  
  results_df$dataset_name <- names(list_of_dataframes)[i]  # Add a column for dataset name
  
  # Append to the main results
  main_results <- rbind(main_results, results_df)
}

# main_results will have one row per input dataset with a column identifying the dataset.

main_results <- main_results %>% select(no_infected, prop_infected, no_infections, 
                       no_1st, no_2nd, no_3rd, inf_predelta, inf_delta,
                       inf_omicron, AR_predelta, AR_delta, AR_omicron, 
                       pcr_neg, dataset_name)
write.csv(main_results, "Summary_sensitivity_results.csv")
