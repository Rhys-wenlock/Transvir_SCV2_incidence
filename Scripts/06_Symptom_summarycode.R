require(tidyverse)


#Extract Full Symptom data
#I.e. all symptoms

#Merge visit data with full list of PCR-positive infections
#requires all_inf_data from incidence_data base
all_inf_data <- read.csv(here::here("Output", "all_inf_data_base.csv"))
visit_data <- read.csv(here::here("Output", "visit_data.csv"))
  
#requires visit data from data_cleaning.R

#Select dates of PCR-positive infections
Inf_dates <- all_inf_data %>% filter(type=="PCR_pos") %>%
  select(Participant_ID, Inf_date)

Full_symp_data <- visit_data_fu %>% 
  left_join(Inf_dates, by="Participant_ID") %>%
  mutate(diff = difftime(as.Date(Date), as.Date(Inf_date), units="days")) %>%
  filter(diff<10 & diff>=0)

colnames(Full_symp_data)[10] <- "Fever"
colnames(Full_symp_data)[11] <- "Cough"
colnames(Full_symp_data)[12] <- "SOB"
colnames(Full_symp_data)[13] <- "Anosmia"
colnames(Full_symp_data)[14] <- "Sore_throat"
colnames(Full_symp_data)[15] <- "Nasal_cong"
colnames(Full_symp_data)[16] <- "CP"
colnames(Full_symp_data)[17] <- "Myalgia"
colnames(Full_symp_data)[18] <- "Headache"
colnames(Full_symp_data)[19] <- "Vom"
colnames(Full_symp_data)[20] <- "Diar"

Full_symp_data <- Full_symp_data %>%
  mutate(Fever = case_when(
    Fever == "Unchecked" ~ 0, 
    Fever == "Checked" ~ 1, 
    TRUE ~NA
  ), 
  Cough = case_when(
    Cough == "Unchecked" ~ 0, 
    Cough == "Checked" ~ 1, 
    TRUE ~NA), 
  SOB = case_when(
    SOB == "Unchecked" ~ 0, 
    SOB == "Checked" ~ 1, 
    TRUE ~NA), 
  Anosmia = case_when(
    Anosmia == "Unchecked" ~ 0, 
    Anosmia == "Checked" ~ 1, 
    TRUE ~NA), 
  Sore_throat = case_when(
      Sore_throat == "Unchecked" ~ 0, 
      Sore_throat == "Checked" ~ 1, 
      TRUE ~NA), 
  Nasal_cong = case_when(
    Nasal_cong == "Unchecked" ~ 0, 
    Nasal_cong == "Checked" ~ 1, 
    TRUE ~NA), 
  CP = case_when(
    CP == "Unchecked" ~ 0, 
    CP == "Checked" ~ 1, 
    TRUE ~NA),
  Myalgia = case_when(
    Myalgia == "Unchecked" ~ 0, 
    Myalgia == "Checked" ~ 1, 
    TRUE ~NA), 
  Diar = case_when(
    Diar == "Unchecked" ~ 0, 
    Diar == "Checked" ~ 1, 
    TRUE ~NA), 
  Headache = case_when(
    Headache == "Unchecked" ~ 0, 
    Headache == "Checked" ~ 1, 
    TRUE ~NA), 
  Vom = case_when(
    Vom == "Unchecked" ~ 0, 
    Vom == "Checked" ~ 1, 
    TRUE ~NA)) %>%
  group_by(Participant_ID, Inf_date) %>%
  mutate(
    Fever=max(Fever),
    Cough = max(Cough), 
    SOB = max(SOB),
    Anosmia = max(Anosmia), 
    Sore_throat = max(Sore_throat), 
    Nasal_cong = max(Nasal_cong), 
    Headache=max(Headache), 
    Vom = max(Vom), 
    Diar=max(Diar), 
    CP=max(CP), 
    Myalgia = max(Myalgia)
  ) %>% filter(row_number()==1) %>%
  mutate(Sum_symp = Fever + Cough + SOB + CP + Anosmia + Headache + Myalgia + Sore_throat + Nasal_cong + Vom + Diar)

#Number of symptoms per infection episode
table(Full_symp_data$Sum_symp)


#Produce summary table
symptom_data <- Full_symp_data %>% ungroup() %>% 
  select(Fever, Cough, Headache, Myalgia, Anosmia, CP, SOB, Nasal_cong, Sore_throat, Vom, Diar)  # replace with your symptom column names
symptom_counts <- colSums(symptom_data, na.rm=TRUE)
summary_table <- data.frame(
  symptom = names(symptom_counts),
  count = as.numeric(symptom_counts)
)

write.csv(summary_table, "Symp_summarytable.csv")



