##Visit Data Cleaning


require(readxl)
require(tidyverse)

#Overall visit cleaning

Aligned_PCR_Symp_12m_data <- read_excel(here::here("Data", "Aligned_Symp_PCRV12.xlsx"))
df <- Aligned_PCR_Symp_12m_data
df <- df %>% group_by(df$`Participant ID`) %>%
  distinct(`Date of Interview`, .keep_all = TRUE)
colnames(df)[7] <- "Ct_val"
colnames(df)[6] <- "PCR_result"
colnames(df)[5] <- "Date"
colnames(df)[1] <- "Participant_ID"
df <- df[,-3]

visit_data <- df
visit_data$Participant_ID[visit_data$Participant_ID=="44-517H"] <- "44-517F"


#Demographic Data Cleaning

#this includes vaccine, serology data for all visits

Transvir_baseline_metadata_immunology_data_V3rw <- read_excel(here::here("Data", "Transvir_baseline_metadata_immunology_data_V3rw.xls"))
demog_data <- Transvir_baseline_metadata_immunology_data_V3rw
demog_data <- demog_data %>% 
  mutate(age_cat = cut(Age, breaks = c(0, 5, 18, 50, Inf), labels = c("<5", "5-17", "18-49", ">50"), right = FALSE)) %>%
  select(Participant_ID, age_cat, Sex, Household_ID, Total_household_no, 
         Adults_over_50, Adults_18_49, Children_5_17, Children_2_4, Children_6m_2, Children_under_6m, Total_rooms)
demog_data$Participant_ID[demog_data$Participant_ID=="44-517H"] <- "44-517F"

TransVir_covariates <- read_excel(here::here("Data","TransVir_covariates.xlsx"))

demog_data <- left_join(demog_data, TransVir_covariates, by="Participant_ID")
demog_data$`Lung Disease`[is.na(demog_data$`Lung Disease`)] <- 0
demog_data$Diabetes[is.na(demog_data$Diabetes)] <- 0
demog_data$Steroid[is.na(demog_data$Steroid)] <- 0
demog_data$Cancer[is.na(demog_data$Cancer)] <- 0
demog_data$HIV[is.na(demog_data$HIV)] <- 0
demog_data$`Other immunsupression`[is.na(demog_data$`Other immunsupression`)] <- 0
demog_data$`Lung Disease`[is.na(demog_data$`Lung Disease`)] <- 0
demog_data$HTN[is.na(demog_data$HTN)] <- 0
demog_data$Smoking[is.na(demog_data$Smoking)] <- 0
demog_data$`Cardiovascular Disease`[is.na(demog_data$`Cardiovascular Disease`)] <- 0
demog_data$Employed[is.na(demog_data$Employed)] <- 0


vaccine_data <- read.csv(here::here("Data", "vaccine_data.csv"))
vaccine_data <- vaccine_data %>%
  select(Participant_ID, vaccinated, vax_date)

vaccine_data$Participant_ID[vaccine_data$Participant_ID=="44-517H"] <- "44-517F"
#19-219C had a vaccine date after V3, assumed 1-year typo, corrected
vaccine_data$vax_date[vaccine_data$Participant_ID=="19-219C"] <- "2021-11-20"

#Serology data cleaning
Transvir_complete_ELISA_data <- read.csv(here::here("Data","Transvir_complete_ELISA_data.csv"))
                                              
All_visit_serol2 <- Transvir_complete_ELISA_data[,c(-1,-4:-5,-7,-9,-11:-25)]
Serol_v1 <- All_visit_serol2[All_visit_serol2$Visits=="V1",]
colnames(Serol_v1) <- c("Participant ID", "Visits", "V1_date", "V1_spike", "V1_ncp")

Serol_v2 <- All_visit_serol2[All_visit_serol2$Visits=="V2",]
colnames(Serol_v2) <- c("Participant ID", "Visits", "V2_date", "V2_spike", "V2_ncp")

Serol_v3 <- All_visit_serol2[All_visit_serol2$Visits=="V3",]
colnames(Serol_v3) <- c("Participant ID", "Visits", "V3_date", "V3_spike", "V3_ncp")

All_visit_serol <- Serol_v1 %>%
  left_join(Serol_v2, by="Participant ID") %>%
  left_join(Serol_v3, by="Participant ID") %>%
  select(-Visits.x, -Visits.y, -Visits)

colnames(All_visit_serol)[1] <- "Participant_ID"
All_visit_serol$Participant_ID[All_visit_serol$Participant_ID=="44-517H"] <- "44-517F"


full_demog_data <- demog_data %>%
  left_join(vaccine_data, by="Participant_ID")





