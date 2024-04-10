
require(readr)
require(tidyverse)
require(lubridate)
require(ggplot2)
require(vcdExtra)

#requires: 
##visit_data from data_cleaning.R
#full_demog_data from data_cleaning.R
#all_visit_serol from data_cleaning.R



# Identifying all infections ----------------------------------------------


# Load data
#visit_data

#remove participants with <5 visits FU
count_df <- visit_data %>%
  filter(`Event Name` %in% c("Weekly Home Visit", "6 Monthly Clinic Visit", "Initial Clinic Visit")) %>%
  group_by(Participant_ID, `Event Name`) %>%
  summarize(Cumulative_Count = n(), .groups = "drop") %>%
  group_by(Participant_ID) %>%
  summarize(Max_Count = max(Cumulative_Count), .groups = "drop")

# Joining the count data back to the original data
visit_data <- left_join(visit_data, count_df, by = "Participant_ID")

visit_data_fu <- visit_data %>%
  filter(Max_Count>=5) %>%
  select(-Max_Count)

# Select positive infection records
pos_data <- visit_data_fu %>% 
  filter(PCR_result=="Positive", 
         !is.na(Date),
         !is.na(`Participant_ID`))


pos_data$type <- "PCR_pos"
pos_data <- pos_data[!is.na(pos_data$`Participant_ID`),]
#pos_data <- bind_rows(pos_data, PCR_neg_inf) #add PCR_neg_inf when we want to
pos_data$PCR_result <- "Positive"


# Get data grouped by ID and arranged by Interview Date
test <- pos_data %>% 
  group_by(`Participant_ID`) %>%
  arrange(`Date`) %>% 
  mutate(num = row_number())

# Store first positive date for each ID
first_pos_date <- test %>% filter(num==1) %>% select(`Participant_ID`, `Date`) %>%
  rename(First_positive_date = `Date`)

# Add First_positive_date to test visit_data and calculate date_diff
test <- test %>%
  left_join(first_pos_date, by= "Participant_ID") %>%
  mutate(date_diff = as.numeric(difftime(`Date`, First_positive_date, units="days")))

# Get first infection data and calculate stats
first_inf <- test %>% filter(date_diff<=28)

# Produce first infection data (infection length, PCR number)
first_inf_data <- first_inf %>%
  group_by(Participant_ID) %>%
  summarize(length_infection = max(date_diff),
            pcr_num = max(num)) %>%
  left_join(first_inf, by="Participant_ID") %>%
  mutate(date_diff = as.numeric(difftime(`Date`, First_positive_date, units="days"))) %>%
  distinct(`Participant_ID`, .keep_all = TRUE) %>%
  rename(Inf_date = First_positive_date) %>%
  select(Participant_ID, `Date`, Inf_date, PCR_result, Ct_val, pcr_num, length_infection, type)

#Produce first infection symptom data
first_inf_symp_data <- visit_data %>%
  left_join(first_pos_date, by="Participant_ID") %>% 
  mutate(date_diff = as.numeric(difftime(`Date`, First_positive_date, units="days"))) %>%
  mutate(symp = case_when(
    is.na(`1 - Have you had ILI symptoms in the last week?`) ~ 1,
    `1 - Have you had ILI symptoms in the last week?`=="Yes" ~ 2,
    TRUE ~ 0)) %>%
  filter(date_diff<10 & date_diff>=0) %>%
  arrange(desc(symp)) %>%
  distinct(`Participant_ID`, .keep_all = TRUE) %>%
  select(Participant_ID, symp)

first_inf_symp_data <- first_inf_symp_data[,-1]

#merge first infection data
first_inf_data <- first_inf_data %>%
  left_join(first_inf_symp_data, by="Participant_ID")%>%
  mutate(inf_num=1)

##SECOND INFECTION##

# Get data grouped by ID and arranged by Interview Date
test <- pos_data %>% 
  left_join(first_pos_date, by="Participant_ID") %>%
  mutate(date_diff = as.numeric(difftime(`Date`, First_positive_date, units="days"))) %>%
  filter(date_diff>28) %>%
  group_by(`Participant_ID`) %>%
  arrange(`Date`) %>% 
  mutate(num = row_number()) 


# Store second positive date for each ID
second_pos_date <- test %>% filter(num==1) %>% select(`Participant_ID`, `Date`) %>%
  rename(second_positive_date = `Date`)

# Add second_positive_date to test visit_data and calculate date_diff
test <- test %>%
  left_join(second_pos_date, by= "Participant_ID") %>%
  mutate(date_diff = as.numeric(difftime(`Date`, second_positive_date, units="days")))

# Get second infection data and calculate stats
second_inf <- test %>% filter(date_diff<=28)

# Produce second infection data (infection length, PCR number)
second_inf_data <- second_inf %>%
  group_by(`Participant_ID`) %>%
  summarize(length_infection = max(date_diff),
            pcr_num = max(num)) %>%
  left_join(second_inf, by="Participant_ID") %>%
  mutate(date_diff = as.numeric(difftime(`Date`, second_positive_date, units="days"))) %>%
  distinct(`Participant_ID`, .keep_all = TRUE) %>%
  rename(Inf_date = second_positive_date) %>%
  select(`Participant_ID`, `Date`, Inf_date, PCR_result, Ct_val, pcr_num, length_infection, type)

#Produce second infection symptom data
second_inf_symp_data <- visit_data %>%
  left_join(second_pos_date, by="Participant_ID") %>% 
  mutate(date_diff = as.numeric(difftime(`Date`, second_positive_date, units="days"))) %>%
  mutate(symp = case_when(
    is.na(`1 - Have you had ILI symptoms in the last week?`) ~ 1,
    `1 - Have you had ILI symptoms in the last week?`=="Yes" ~ 2,
    TRUE ~ 0)) %>%
  filter(date_diff<10 & date_diff>=0) %>%
  arrange(desc(symp)) %>%
  distinct(`Participant_ID`, .keep_all = TRUE) %>%
  select(`Participant_ID`, symp)

second_inf_symp_data <- second_inf_symp_data[,-1]


#merge second infection data
second_inf_data <- second_inf_data %>%
  left_join(second_inf_symp_data, by="Participant_ID") %>%
  mutate(inf_num=2)


##third INFECTION##

# Get data grouped by ID and arranged by Interview Date
test <- pos_data %>% 
  left_join(second_pos_date, by="Participant_ID") %>%
  mutate(date_diff = as.numeric(difftime(`Date`, second_positive_date, units="days"))) %>%
  filter(date_diff>28) %>%
  group_by(`Participant_ID`) %>%
  arrange(`Date`) %>% 
  mutate(num = row_number()) 


# Store third positive date for each ID
third_pos_date <- test %>% filter(num==1) %>% select(`Participant_ID`, `Date`) %>%
  rename(third_positive_date = `Date`)

# Add third_positive_date to test visit_data and calculate date_diff
test <- test %>%
  left_join(third_pos_date, by= "Participant_ID") %>%
  mutate(date_diff = as.numeric(difftime(`Date`, third_positive_date, units="days")))

# Get third infection data and calculate stats
third_inf <- test %>% filter(date_diff<=28)

# Produce third infection data (infection length, PCR number)
third_inf_data <- third_inf %>%
  group_by(`Participant_ID`) %>%
  summarize(length_infection = max(date_diff),
            pcr_num = max(num)) %>%
  left_join(third_inf, by="Participant_ID") %>%
  mutate(date_diff = as.numeric(difftime(`Date`, third_positive_date, units="days"))) %>%
  distinct(`Participant_ID`, .keep_all = TRUE) %>%
  rename(Inf_date = third_positive_date) %>%
  select(`Participant_ID`, `Date`,  Inf_date, PCR_result, Ct_val, pcr_num, length_infection, type)

#Produce third infection symptom data
third_inf_symp_data <- visit_data %>%
  left_join(third_pos_date, by="Participant_ID") %>% 
  mutate(date_diff = as.numeric(difftime(`Date`, third_positive_date, units="days"))) %>%
  mutate(symp = case_when(
    is.na(`1 - Have you had ILI symptoms in the last week?`) ~ 1,
    `1 - Have you had ILI symptoms in the last week?`=="Yes" ~ 2,
    TRUE ~ 0)) %>%
  filter(date_diff<10 & date_diff>=0) %>%
  arrange(desc(symp)) %>%
  distinct(`Participant_ID`, .keep_all = TRUE) %>%
  select(`Participant_ID`, symp)

third_inf_symp_data <- third_inf_symp_data[,-1]


#merge third infection data
third_inf_data <- third_inf_data %>%
  left_join(third_inf_symp_data, by="Participant_ID")%>%
  mutate(inf_num=3)


#Produce no infection data 

no_inf_data <- visit_data_fu %>% 
  group_by(`Participant_ID`) %>% 
  filter(row_number() == first(which(visit_data$PCR_result=="Negative"))) %>%
  mutate(Inf_date = NA,
         symp = NA,
         length_infection = NA,
         pcr_num = NA,
         inf_num = 0, 
         type = "None",
         PCR_result ="Negative") %>%
  select(`Participant_ID`, `Date`, Inf_date ,PCR_result, Ct_val, pcr_num, length_infection,inf_num, symp, type)


#Merge all infection episode data together

all_inf_data <- rbind(first_inf_data, second_inf_data, third_inf_data, no_inf_data)

#Produce infection numbers/total infection number
all_inf_data <- all_inf_data %>%
  group_by(`Participant_ID`) %>%
  mutate(total_inf = max(inf_num))

#Produce dataframe of infection dates per participant
inf_dates <- first_pos_date %>%
  left_join(second_pos_date, by="Participant_ID") %>%
  left_join(third_pos_date, by="Participant_ID") 


# Setting FU start and end dates ------------------------------------------


# Create last_date dataframe
last_date <- visit_data %>%
  group_by(`Participant_ID`) %>%
  summarize(max_date = max(as.Date(`Date`), na.rm=TRUE))

all_inf_data <- all_inf_data %>%
  left_join(inf_dates, by="Participant_ID") %>%
  left_join(All_visit_serol, by="Participant_ID")  %>%
  left_join(last_date, by="Participant_ID")






all_inf_data <- all_inf_data %>% 
  mutate(across(c(V1_date, First_positive_date, second_positive_date, third_positive_date, Inf_date), as.Date)) %>%
  mutate(start_date = case_when(
    total_inf == 0 ~ V1_date,
    inf_num == 1 ~ V1_date,
    inf_num == 0 & total_inf == 1 ~ First_positive_date + 28,
    inf_num == 2 ~ First_positive_date + 28,
    inf_num == 0 & total_inf == 2 ~ second_positive_date + 28,
    inf_num == 3 ~ second_positive_date + 28,
    inf_num == 0 & total_inf == 3 ~ third_positive_date + 28,
    TRUE ~ NA
  )
  ) %>%
  mutate(
    end_date = case_when(
      total_inf == 0 ~ max_date,
      inf_num == 1 ~ Inf_date,
      inf_num == 0 & total_inf == 1 ~ max_date,
      inf_num == 2 ~ Inf_date,
      inf_num == 0 & total_inf == 2 ~ max_date,
      inf_num == 3 ~ Inf_date,
      inf_num == 0 & total_inf == 3 ~ max_date,
      TRUE ~ as.Date(NA)  # To handle other cases if any
    )
  ) %>%
  select(-max_date)
# Convert end_date to Date type
all_inf_data <- all_inf_data %>% mutate(end_date = as.Date(end_date, origin="1970-01-01"))

#Produce code for combined prior infection variable
#requires infections_by_week_start from 02

all_inf_data$begin_date <- as.Date("2020-03-15")

sample_weeks <- function(begin_date, start_date, infections_by_week_start) {
  sample_date <- tryCatch(
    {
      sampled_week <- infections_by_week_start %>%
        filter(week >= as.Date(begin_date) & week <= as.Date(start_date)) %>%
        mutate(prob = infections / sum(infections)) %>%
        sample_n(1, replace = TRUE, weight = prob) %>%
        pull(week)
      
      return(sampled_week)
    },
    error = function(e) {
      print(paste("Error:", e$message))
      return(NA)
    }
  )
  
  return(sample_date)
}

all_inf_data <- all_inf_data %>%
  rowwise() %>%
  mutate(latest_infection = case_when(
    inf_num == 0 & is.na(First_positive_date) & is.na(second_positive_date) & is.na(third_positive_date) & V1_spike =="Pos" ~ as.Date(sample_weeks(begin_date, start_date, infections_by_week_start)),
    inf_num == 0 ~ pmax(First_positive_date, second_positive_date, third_positive_date, na.rm=TRUE),
    inf_num == 1 & V1_spike == "Pos" ~ as.Date(sample_weeks(begin_date, start_date, infections_by_week_start)),
    inf_num == 1 & V1_spike == "Neg" ~ as.Date(NA), 
    inf_num == 2 ~ First_positive_date, 
    inf_num == 3 ~ second_positive_date, 
    inf_num == 4 ~ third_positive_date, 
    TRUE ~ NA
  ))


write.csv(all_inf_data, "all_inf_data_8.csv")



