###Set up base assumption dataset for incidence modelling
#Median FC boost in SN/N in vacc

require(readr)
require(tidyverse)
require(lubridate)
require(ggplot2)
require(vcdExtra)

#requires: 
##visit_data from data_cleaning.R
#full_demog_data from data_cleaning.R
#all_visit_serol from data_cleaning.R
#sample weeks function from PCR_negative_dates.R

# PCR Negative Infections -------------------------------------------------


#load infections_by_week_start

#select the PCR negative cohort
inf_vaccine_serol_FC_countinf <- read_csv(here::here("Data","inf_vaccine_serol_medianS_medianNvacc.csv"))
PCR_neg_inf <- inf_vaccine_serol_FC_countinf[inf_vaccine_serol_FC_countinf$PCR_neg_booster_0_6=="Y"|inf_vaccine_serol_FC_countinf$PCR_neg_booster_6_12=="Y"|
                                               inf_vaccine_serol_FC_countinf$PCR_serol_0_6=="PCR_neg_0_6_seroconv"|inf_vaccine_serol_FC_countinf$PCR_serol_6_12=="PCR_neg_6_12_seroconv",]


# Subset and select columns
PCR_neg_inf <- PCR_neg_inf[!is.na(PCR_neg_inf$Participant_ID), ]
PCR_neg_inf <- PCR_neg_inf[, c(1, 19:21, 33:34, 39:40)]

# Create binary variables based on conditions
PCR_neg_inf$inf0_6 <- as.integer(PCR_neg_inf$PCR_serol_0_6 == "PCR_neg_0_6_seroconv" | PCR_neg_inf$PCR_neg_booster_0_6 == "Y")
PCR_neg_inf$inf6_12 <- as.integer(PCR_neg_inf$PCR_serol_6_12 == "PCR_neg_6_12_seroconv" | PCR_neg_inf$PCR_neg_booster_6_12 == "Y")

# Subset data for inf0_6
PCR_neg_inf06 <- PCR_neg_inf[PCR_neg_inf$inf0_6 == 1 & !is.na(PCR_neg_inf$inf0_6), ]

# Generate random dates
PCR_neg_inf06$start_date <- as.Date(PCR_neg_inf06$V1_date)
PCR_neg_inf06$end_date <- as.Date(PCR_neg_inf06$V2_date)

# Function to sample weeks (from PCR_neg_dates_code.R)
sample_weeks <- function(start_date, end_date, size = 1) {
  infections_by_week_start %>%
    filter(week >= start_date & week <= end_date) %>%
    mutate(prob = infections / sum(infections)) %>%
    { sample(.$week, size, replace = TRUE, prob = .$prob) }
}

# Sample weeks for each participant
PCR_neg_inf06 <- PCR_neg_inf06 %>%
  mutate(sampled_week = map2(start_date, end_date, sample_weeks)) %>%
  unnest(sampled_week) %>%
  mutate(sampled_week = as.Date(sampled_week, origin = "1970-01-01"))

# Subset data for inf6_12
PCR_neg_inf612 <- PCR_neg_inf[PCR_neg_inf$inf6_12 == 1 & !is.na(PCR_neg_inf$inf6_12), ]
PCR_neg_inf612$start_date <- as.Date(PCR_neg_inf612$V2_date)
PCR_neg_inf612$end_date <- as.Date(PCR_neg_inf612$V3_date)

# Sample weeks for each participant
PCR_neg_inf612 <- PCR_neg_inf612 %>%
  mutate(sampled_week = map2(start_date, end_date, sample_weeks)) %>%
  unnest(sampled_week) %>%
  mutate(sampled_week = as.Date(sampled_week, origin = "1970-01-01"))

# Combine dataframes
PCR_neg_inf <- rbind(PCR_neg_inf06, PCR_neg_inf612)
PCR_neg_inf <- PCR_neg_inf[, c(1, 13)]
colnames(PCR_neg_inf) <- c("Participant_ID", "Date")
PCR_neg_inf$type <- "PCR_neg"


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
pos_data <- bind_rows(pos_data, PCR_neg_inf) #add PCR_neg_inf when we want to
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
    is.na(`1 - Have you had ILI symptoms in the last week?`) ~ 2,
    `1 - Have you had ILI symptoms in the last week?`=="Yes" ~ 1,
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
    is.na(`1 - Have you had ILI symptoms in the last week?`) ~ 2,
    `1 - Have you had ILI symptoms in the last week?`=="Yes" ~ 1,
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
    is.na(`1 - Have you had ILI symptoms in the last week?`) ~ 2,
    `1 - Have you had ILI symptoms in the last week?`=="Yes" ~ 1,
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


write.csv(all_inf_data, "all_inf_data_3.csv")



# Producing PCR Negative Plot ---------------------------------------------

#load gambia_case_linelist

colnames(gambia_case_linelist)[1] <- "Type"
gambia_case_linelist$date <- as.Date(gambia_case_linelist$date)

PCR_neg_freq <- PCR_neg_inf[,c(3,2)]
colnames(PCR_neg_freq)[1] <- "Type"
colnames(PCR_neg_freq)[2] <- "date"
PCR_neg_freq$date <- as.Date(PCR_neg_freq$date)


epi_curve <- rbind(gambia_case_linelist, PCR_neg_freq)
epi_curve$date <- as.Date(epi_curve$date, origin="1970-01-01")
weekly_breaks <- seq.Date(from = as.Date("2021-03-02"),
                          to = as.Date("2022-05-31"),
                          by = "week")

supp_figure_1d <- ggplot(epi_curve) + 
  
  geom_histogram(
    mapping = aes(
      x = date, color = Type),    # arguments inside aes() apply by group
    stat = "count",
    # arguments outside aes() apply to all data
    
    # histogram breaks
    breaks = weekly_breaks, # pre-defined date vector (see earlier in this page)
    closed = "left" # count cases from start of breakpoint
  ) + 
  
  facet_wrap(
    Type~., ncol=1, scales = "free_y") + ylab("") + xlab("Week") + 
  theme(
    # Remove panel border
    panel.border = element_blank(), panel.background = element_blank(), panel.grid = element_line(colour = "grey", size = 0.1),
    axis.line = element_line(colour = "grey"), legend.position = "bottom", legend.title = element_blank(),
    axis.title.x = element_text(vjust=-0.5)
  ) + scale_color_manual(values=cols)

ggsave("supp_fig_1d.png", supp_figure_1d, width = 10, height = 6, dpi=600)
