require(dplyr)
require(lubridate)
require(tidvyerse)

#REQUIRES: 
#run Data_cleaning.R


#set serial interval
time_window <- 14


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


# Cluster 1 ---------------------------------------------------------------

# Join and basic transformations
pos_HH_data <- pos_data %>%
  left_join(full_demog_data, by="Participant_ID") %>%
  mutate(Household_ID = as.factor(Household_ID)) #join demography data to positive visit data

PCR_data <- visit_data_fu %>%
  left_join(full_demog_data, by="Participant_ID") %>%
  mutate(Household_ID = as.factor(Household_ID)) #join demography data to all visit data

# Calculate index data
index_data1 <- pos_data %>%
  left_join(full_demog_data, by="Participant_ID") %>% #join pos data to demog data
  mutate(Household_ID = as.factor(Household_ID)) %>%
  group_by(Household_ID) %>%
  summarize(first_case_date = min(Date)) %>% #identify the first case in the household
  left_join(pos_HH_data, by="Household_ID") %>% #join all other positive 
  filter(Date == first_case_date) %>% #select only participants whos date is equal to first case date
  mutate(
    first_end_date = as.Date(first_case_date) + time_window,
    index_interval = interval(first_case_date, first_end_date)
  )  %>% 
  select(Household_ID, first_case_date, Participant_ID, Ct_val, first_end_date, index_interval) 

#End result = participants whos first positive date = households first positive date, with an interval of exposure

# Primary data
primary_1_data <- PCR_data %>% #using all visit data 
  left_join(select(index_data1, Household_ID, index_interval), by="Household_ID") %>% #left join the household ID and index interval from index1
  filter(Date %within% index_interval & !is.na(Participant_ID)) %>% #keep only dates within index interval
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>% #order PCRs so that positive are first
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1) %>% #keep the first result, ensuring that positives within interval are kept before negatives
  anti_join(index_data1, by="Participant_ID") %>% #remove index cases from it 
  select(Household_ID, Participant_ID, PCR_result, Ct_val, Date)

# Positive primary data and primary dates
pos_primary_1_data <- filter(primary_1_data, PCR_result == "Positive") #from primary exposures, keep only those positive

primary_1_dates <- pos_primary_1_data %>%
  group_by(Household_ID) %>%
  summarize(primary1_max_date = max(Date)) %>% #identify last positive date amongst primaries within household
  mutate(
    primary_1_enddate = as.Date(primary1_max_date) + time_window,
    primary1_interval = interval(primary1_max_date, primary_1_enddate)
  ) #set up the interval from which household participants can be exposed to primary (non-index) cases

# Secondary data
secondary_1_data <- PCR_data %>%
  left_join(primary_1_dates, by="Household_ID") %>%
  filter(Date %within% primary1_interval & !is.na(Participant_ID)) %>% #identify participants with visits in the exposure interval of primaries
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>%
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1)  %>%#keep one result for each participant exposed to a primary, positive >negative
  select(Household_ID, Participant_ID, PCR_result, Ct_val, Date) 
#those exposed to a primary and positive are secondaries

# Positive secondary data and secondary dates
pos_secondary_1_data <- filter(secondary_1_data, PCR_result == "Positive") #keep only positive secondaries

secondary_1_dates <- pos_secondary_1_data %>%
  group_by(Household_ID) %>%
  summarize(secondary1_max_date = max(Date)) %>%
  mutate(
    secondary_1_enddate = as.Date(secondary1_max_date) + time_window,
    secondary1_interval = interval(secondary1_max_date, secondary_1_enddate)
  ) #identify latest positive date amongst secondaries and produce exposure interval

# Tertiary data
tertiary_1_data <- PCR_data %>%
  left_join(secondary_1_dates, by="Household_ID") %>%
  filter(Date %within% secondary1_interval & !is.na(Participant_ID)) %>%
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>%
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1) %>%#identify any tertiary cases, exposed to secondaries
  select(Household_ID, Participant_ID, PCR_result, Ct_val, Date)  

#remove unnecessary columns
prim_1_data <- primary_1_data
sec_1_data <- secondary_1_data
tert_1_data <- tertiary_1_data
ind_1_data <- index_data1

#Keep one row per exposed participant, keeping positive > Neg >NA
cluster_1_data <- prim_1_data %>%
  bind_rows(sec_1_data, tert_1_data) %>%
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative",NA))) %>% 
  group_by(Participant_ID) %>% 
  arrange(PCR_result, Date) %>%
  filter(row_number()==1) 

ind_1_data <- ind_1_data  %>% group_by(Household_ID) %>% mutate(count = n_distinct(Participant_ID)) #count number of index's per household 
single_index_HH_1 <- ind_1_data %>% 
  filter(count==1) %>%
  select(-count) %>%
  mutate(
    Index_ID = Participant_ID, 
    Date = first_case_date
  )

cluster_1_data_sing <- single_index_HH_1 %>%
  left_join(cluster_1_data, by="Household_ID") %>%
  select(-Participant_ID.x, -Ct_val.y) %>%
  rename(
    Participant_ID = Participant_ID.y,
    Index_Ct = Ct_val.x,
    Date = Date.y
  )
#results in single index clusters with index data included 


#separate positives to calculate end dates
all_pos <- pos_HH_data %>%
  group_by(Participant_ID) %>%
  arrange(Date) %>%
  mutate(num=1:n()) %>% 
  select(Participant_ID, Household_ID, Date)


#need to define an end date to the cluster
end_cluster1 <- cluster_1_data_sing %>%
  select(Household_ID, Participant_ID, PCR_result, Date) %>%
  filter(PCR_result == "Positive") %>%
  bind_rows(single_index_HH_1) %>%
  select(Household_ID, Participant_ID, Date) %>%
  left_join(all_pos, by=c("Participant_ID", "Household_ID")) %>%
  mutate(date_diff = difftime(Date.y, Date.x, units="days")) %>%
  filter(date_diff <=28) %>%
  group_by(Participant_ID) %>%
  mutate(cluster_end_date = max(Date.y)) %>%
  group_by(Household_ID) %>%
  summarize(cluster_end_date = max(cluster_end_date))


#Symptom data 
PCR_index_1 <- PCR_data %>%
  left_join(single_index_HH_1, by="Participant_ID") %>%
  filter(!is.na(first_case_date)) %>%
  mutate(date_diff = difftime(as.Date(Date.x), as.Date(first_case_date), units="days")) %>%
  filter(date_diff <10 & date_diff >=0) %>%
  filter(!is.na(Participant_ID)) %>%
  mutate(symp = case_when(
    is.na(`1 - Have you had ILI symptoms in the last week?`) ~ "Missing",
    `1 - Have you had ILI symptoms in the last week?` == "Yes" ~ "1",
    TRUE ~ "0"
  )) %>% arrange(factor(symp, levels=c(1, "Missing", 0))) %>%
                   filter(row_number()==1) %>%
  select(Participant_ID, symp)

single_index_HH_1 <- single_index_HH_1 %>% 
  left_join(PCR_index_1, by="Participant_ID") %>%
  select(Household_ID, Index_ID, Ct_val, symp, Date)  %>%
  rename(
    Index_date = Date
  )


cluster_1_data_sing <- single_index_HH_1 %>%
  left_join(cluster_1_data, by="Household_ID") %>%
  select(-Ct_val.y) %>%
  rename(
    Index_Ct = Ct_val.x
  ) %>%
  mutate(Cluster_no = 1)

cluster_1_data_sing_prim <- single_index_HH_1 %>%
  left_join(prim_1_data, by="Household_ID") %>%
  select(-Ct_val.y) %>%
  rename(
    Index_Ct = Ct_val.x
  ) %>%
  mutate(Cluster_no = 1)

# Cluster 2 ---------------------------------------------------------------

pos_HH_data <- pos_data %>%
  left_join(full_demog_data, by="Participant_ID") %>%
  mutate(Household_ID = as.factor(Household_ID))

PCR_data <- visit_data_fu %>%
  left_join(full_demog_data, by="Participant_ID") %>%
  mutate(Household_ID = as.factor(Household_ID))

# Calculate index data
index_data2 <- pos_data %>%
  left_join(full_demog_data, by="Participant_ID") %>% #join pos data to demog data
  mutate(Household_ID = as.factor(Household_ID))%>%
  left_join(end_cluster1, by="Household_ID") %>%
  filter(Date > cluster_end_date) %>%
  group_by(Household_ID) %>%
  summarize(first_case_date = min(Date)) %>% #identify the first case in the household
  left_join(pos_HH_data, by="Household_ID") %>% #join all other positive 
  filter(Date == first_case_date) %>% #select only participants whos date is equal to first case date
  mutate(
    first_end_date = as.Date(first_case_date) + time_window,
    index_interval = interval(first_case_date, first_end_date)
  )  %>% 
  select(Household_ID, first_case_date, Participant_ID, Ct_val, first_end_date, index_interval) 


# Primary data
primary_2_data <- PCR_data %>% #using all visit data 
  left_join(select(index_data2, Household_ID, index_interval), by="Household_ID") %>% #left join the household ID and index interval from index2
  filter(Date %within% index_interval & !is.na(Participant_ID)) %>% #keep only dates within index interval
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>% #order PCRs so that positive are first
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1) %>% #keep the first result, ensuring that positives within interval are kept before negatives
  anti_join(index_data2, by="Participant_ID") %>% #remove index cases from it 
  select(Household_ID, Participant_ID, PCR_result, Ct_val, Date)

# Positive primary data and primary dates
pos_primary_2_data <- filter(primary_2_data, PCR_result == "Positive") #from primary exposures, keep only those positive

primary_2_dates <- pos_primary_2_data %>%
  group_by(Household_ID) %>%
  summarize(primary2_max_date = max(Date)) %>% #identify last positive date amongst primaries within household
  mutate(
    primary_2_enddate = as.Date(primary2_max_date) + time_window,
    primary2_interval = interval(primary2_max_date, primary_2_enddate)
  ) #set up the interval from which household participants can be exposed to primary (non-index) cases

# Secondary data
secondary_2_data <- PCR_data %>%
  left_join(primary_2_dates, by="Household_ID") %>%
  filter(Date %within% primary2_interval & !is.na(Participant_ID)) %>% #identify participants with visits in the exposure interval of primaries
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>%
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1)  %>%#keep one result for each participant exposed to a primary, positive >negative
  select(Household_ID, Participant_ID, PCR_result, Ct_val, Date) 
#those exposed to a primary and positive are secondaries

# Positive secondary data and secondary dates
pos_secondary_2_data <- filter(secondary_2_data, PCR_result == "Positive") #keep only positive secondaries

secondary_2_dates <- pos_secondary_2_data %>%
  group_by(Household_ID) %>%
  summarize(secondary2_max_date = max(Date)) %>%
  mutate(
    secondary_2_enddate = as.Date(secondary2_max_date) + time_window,
    secondary2_interval = interval(secondary2_max_date, secondary_2_enddate)
  ) #identify latest positive date amongst secondaries and produce exposure interval

# Tertiary data
tertiary_2_data <- PCR_data %>%
  left_join(secondary_2_dates, by="Household_ID") %>%
  filter(Date %within% secondary2_interval & !is.na(Participant_ID)) %>%
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>%
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1) %>%#identify any tertiary cases, exposed to secondaries
  select(Household_ID, Participant_ID, PCR_result, Ct_val, Date)  

#remove unnecessary columns
prim_2_data <- primary_2_data
sec_2_data <- secondary_2_data
tert_2_data <- tertiary_2_data
ind_2_data <- index_data2

#Keep one row per exposed participant, keeping positive > Neg >NA
cluster_2_data <- prim_2_data %>%
  bind_rows(sec_2_data, tert_2_data) %>%
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative",NA))) %>% 
  group_by(Participant_ID) %>% 
  arrange(PCR_result, Date) %>%
  filter(row_number()==1) 

ind_2_data <- ind_2_data  %>% group_by(Household_ID) %>% mutate(count = n_distinct(Participant_ID)) #count number of index's per household 
single_index_HH_2 <- ind_2_data %>% 
  filter(count==1) %>%
  select(-count) %>%
  mutate(
    Index_ID = Participant_ID, 
    Date = first_case_date
  )

cluster_2_data_sing <- single_index_HH_2 %>%
  left_join(cluster_2_data, by="Household_ID") %>%
  select(-Participant_ID.x, -Ct_val.y) %>%
  rename(
    Participant_ID = Participant_ID.y,
    Index_Ct = Ct_val.x,
    Date = Date.y
  )
#results in single index clusters with index data included 


#separate positives to calculate end dates
all_pos <- pos_HH_data %>%
  group_by(Participant_ID) %>%
  arrange(Date) %>%
  mutate(num=1:n()) %>% 
  select(Participant_ID, Household_ID, Date)


#need to define an end date to the cluster
end_cluster2 <- cluster_2_data_sing %>%
  select(Household_ID, Participant_ID, PCR_result, Date) %>%
  filter(PCR_result == "Positive") %>%
  bind_rows(single_index_HH_2) %>%
  select(Household_ID, Participant_ID, Date) %>%
  left_join(all_pos, by=c("Participant_ID", "Household_ID")) %>%
  mutate(date_diff = difftime(Date.y, Date.x, units="days")) %>%
  filter(date_diff <=28) %>%
  group_by(Participant_ID) %>%
  mutate(cluster_end_date = max(Date.y)) %>%
  group_by(Household_ID) %>%
  summarize(cluster_end_date = max(cluster_end_date))


#Symptom data 
PCR_index_2 <- PCR_data %>%
  left_join(single_index_HH_2, by="Participant_ID") %>%
  filter(!is.na(first_case_date)) %>%
  mutate(date_diff = difftime(as.Date(Date.x), as.Date(first_case_date), units="days")) %>%
  filter(date_diff <10 & date_diff >=0) %>%
  filter(!is.na(Participant_ID)) %>%
  mutate(symp = case_when(
    is.na(`1 - Have you had ILI symptoms in the last week?`) ~ "Missing",
    `1 - Have you had ILI symptoms in the last week?` == "Yes" ~ "1",
    TRUE ~ "0"
  )) %>% arrange(factor(symp, levels=c(1, "Missing", 0))) %>%
  filter(row_number()==1) %>%
  select(Participant_ID, symp)

single_index_HH_2 <- single_index_HH_2 %>% 
  left_join(PCR_index_2, by="Participant_ID") %>%
  select(Household_ID, Index_ID, Ct_val, symp, Date)  %>%
  rename(
    Index_date = Date
  )



cluster_2_data_sing <- single_index_HH_2 %>%
  left_join(cluster_2_data, by="Household_ID") %>%
  select(-Ct_val.y) %>%
  rename(
    Index_Ct = Ct_val.x
  ) %>%
  mutate(Cluster_no = 2)

cluster_2_data_sing_prim <- single_index_HH_2 %>%
  left_join(prim_2_data, by="Household_ID") %>%
  select(-Ct_val.y) %>%
  rename(
    Index_Ct = Ct_val.x
  ) %>%
  mutate(Cluster_no = 2)




# Cluster 3 ---------------------------------------------------------------

pos_HH_data <- pos_data %>%
  left_join(full_demog_data, by="Participant_ID") %>%
  mutate(Household_ID = as.factor(Household_ID))

PCR_data <- visit_data_fu %>%
  left_join(full_demog_data, by="Participant_ID") %>%
  mutate(Household_ID = as.factor(Household_ID))

# Calculate index data
index_data3 <- pos_data %>%
  left_join(full_demog_data, by="Participant_ID") %>% #join pos data to demog data
  mutate(Household_ID = as.factor(Household_ID))%>%
  left_join(end_cluster2, by="Household_ID") %>%
  filter(Date > cluster_end_date) %>%
  group_by(Household_ID) %>%
  summarize(first_case_date = min(Date)) %>% #identify the first case in the household
  left_join(pos_HH_data, by="Household_ID") %>% #join all other positive 
  filter(Date == first_case_date) %>% #select only participants whos date is equal to first case date
  mutate(
    first_end_date = as.Date(first_case_date) + time_window,
    index_interval = interval(first_case_date, first_end_date)
  )  %>% 
  select(Household_ID, first_case_date, Participant_ID, Ct_val, first_end_date, index_interval) 


# Primary data
primary_3_data <- PCR_data %>% #using all visit data 
  left_join(select(index_data3, Household_ID, index_interval), by="Household_ID") %>% #left join the household ID and index interval from index3
  filter(Date %within% index_interval & !is.na(Participant_ID)) %>% #keep only dates within index interval
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>% #order PCRs so that positive are first
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1) %>% #keep the first result, ensuring that positives within interval are kept before negatives
  anti_join(index_data3, by="Participant_ID") %>% #remove index cases from it 
  select(Household_ID, Participant_ID, PCR_result, Ct_val, Date)

# Positive primary data and primary dates
pos_primary_3_data <- filter(primary_3_data, PCR_result == "Positive") #from primary exposures, keep only those positive

primary_3_dates <- pos_primary_3_data %>%
  group_by(Household_ID) %>%
  summarize(primary3_max_date = max(Date)) %>% #identify last positive date amongst primaries within household
  mutate(
    primary_3_enddate = as.Date(primary3_max_date) + time_window,
    primary3_interval = interval(primary3_max_date, primary_3_enddate)
  ) #set up the interval from which household participants can be exposed to primary (non-index) cases

# Secondary data
secondary_3_data <- PCR_data %>%
  left_join(primary_3_dates, by="Household_ID") %>%
  filter(Date %within% primary3_interval & !is.na(Participant_ID)) %>% #identify participants with visits in the exposure interval of primaries
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>%
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1)  %>%#keep one result for each participant exposed to a primary, positive >negative
  select(Household_ID, Participant_ID, PCR_result, Ct_val, Date) 
#those exposed to a primary and positive are secondaries

# Positive secondary data and secondary dates
pos_secondary_3_data <- filter(secondary_3_data, PCR_result == "Positive") #keep only positive secondaries

secondary_3_dates <- pos_secondary_3_data %>%
  group_by(Household_ID) %>%
  summarize(secondary3_max_date = max(Date)) %>%
  mutate(
    secondary_3_enddate = as.Date(secondary3_max_date) + time_window,
    secondary3_interval = interval(secondary3_max_date, secondary_3_enddate)
  ) #identify latest positive date amongst secondaries and produce exposure interval

# Tertiary data
tertiary_3_data <- PCR_data %>%
  left_join(secondary_3_dates, by="Household_ID") %>%
  filter(Date %within% secondary3_interval & !is.na(Participant_ID)) %>%
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>%
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1) %>%#identify any tertiary cases, exposed to secondaries
  select(Household_ID, Participant_ID, PCR_result, Ct_val, Date)  

#remove unnecessary columns
prim_3_data <- primary_3_data
sec_3_data <- secondary_3_data
tert_3_data <- tertiary_3_data
ind_3_data <- index_data3

#Keep one row per exposed participant, keeping positive > Neg >NA
cluster_3_data <- prim_3_data %>%
  bind_rows(sec_3_data, tert_3_data) %>%
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative",NA))) %>% 
  group_by(Participant_ID) %>% 
  arrange(PCR_result, Date) %>%
  filter(row_number()==1) 

ind_3_data <- ind_3_data  %>% group_by(Household_ID) %>% mutate(count = n_distinct(Participant_ID)) #count number of index's per household 
single_index_HH_3 <- ind_3_data %>% 
  filter(count==1) %>%
  select(-count) %>%
  mutate(
    Index_ID = Participant_ID, 
    Date = first_case_date
  )

cluster_3_data_sing <- single_index_HH_3 %>%
  left_join(cluster_3_data, by="Household_ID") %>%
  select(-Participant_ID.x, -Ct_val.y) %>%
  rename(
    Participant_ID = Participant_ID.y,
    Index_Ct = Ct_val.x,
    Date = Date.y
  )
#results in single index clusters with index data included 


#separate positives to calculate end dates
all_pos <- pos_HH_data %>%
  group_by(Participant_ID) %>%
  arrange(Date) %>%
  mutate(num=1:n()) %>% 
  select(Participant_ID, Household_ID, Date)


#need to define an end date to the cluster
end_cluster3 <- cluster_3_data_sing %>%
  select(Household_ID, Participant_ID, PCR_result, Date) %>%
  filter(PCR_result == "Positive") %>%
  bind_rows(single_index_HH_3) %>%
  select(Household_ID, Participant_ID, Date) %>%
  left_join(all_pos, by=c("Participant_ID", "Household_ID")) %>%
  mutate(date_diff = difftime(Date.y, Date.x, units="days")) %>%
  filter(date_diff <=28) %>%
  group_by(Participant_ID) %>%
  mutate(cluster_end_date = max(Date.y)) %>%
  group_by(Household_ID) %>%
  summarize(cluster_end_date = max(cluster_end_date))


#Symptom data 
PCR_index_3 <- PCR_data %>%
  left_join(single_index_HH_3, by="Participant_ID") %>%
  filter(!is.na(first_case_date)) %>%
  mutate(date_diff = difftime(as.Date(Date.x), as.Date(first_case_date), units="days")) %>%
  filter(date_diff <10 & date_diff >=0) %>%
  filter(!is.na(Participant_ID)) %>%
  mutate(symp = case_when(
    is.na(`1 - Have you had ILI symptoms in the last week?`) ~ "Missing",
    `1 - Have you had ILI symptoms in the last week?` == "Yes" ~ "1",
    TRUE ~ "0"
  )) %>% arrange(factor(symp, levels=c(1, "Missing", 0))) %>%
  filter(row_number()==1) %>%
  select(Participant_ID, symp)

single_index_HH_3 <- single_index_HH_3 %>% 
  left_join(PCR_index_3, by="Participant_ID") %>%
  select(Household_ID, Index_ID, Ct_val, symp, Date)  %>%
  rename(
    Index_date = Date
  )



cluster_3_data_sing <- single_index_HH_3 %>%
  left_join(cluster_3_data, by="Household_ID") %>%
  select(-Ct_val.y) %>%
  rename(
    Index_Ct = Ct_val.x
  ) %>%
  mutate(Cluster_no = 3)

cluster_3_data_sing_prim <- single_index_HH_3 %>%
  left_join(prim_3_data, by="Household_ID") %>%
  select(-Ct_val.y) %>%
  rename(
    Index_Ct = Ct_val.x
  ) %>%
  mutate(Cluster_no = 3)





# Cluster 4 ---------------------------------------------------------------

pos_HH_data <- pos_data %>%
  left_join(full_demog_data, by="Participant_ID") %>%
  mutate(Household_ID = as.factor(Household_ID))

PCR_data <- visit_data_fu %>%
  left_join(full_demog_data, by="Participant_ID") %>%
  mutate(Household_ID = as.factor(Household_ID))

# Calculate index data
index_data4 <- pos_data %>%
  left_join(full_demog_data, by="Participant_ID") %>% #join pos data to demog data
  mutate(Household_ID = as.factor(Household_ID))%>%
  left_join(end_cluster3, by="Household_ID") %>%
  filter(Date > cluster_end_date) %>%
  group_by(Household_ID) %>%
  summarize(first_case_date = min(Date)) %>% #identify the first case in the household
  left_join(pos_HH_data, by="Household_ID") %>% #join all other positive 
  filter(Date == first_case_date) %>% #select only participants whos date is equal to first case date
  mutate(
    first_end_date = as.Date(first_case_date) + time_window,
    index_interval = interval(first_case_date, first_end_date)
  )  %>% 
  select(Household_ID, first_case_date, Participant_ID, Ct_val, first_end_date, index_interval) 


# Primary data
primary_4_data <- PCR_data %>% #using all visit data 
  left_join(select(index_data4, Household_ID, index_interval), by="Household_ID") %>% #left join the household ID and index interval from index4
  filter(Date %within% index_interval & !is.na(Participant_ID)) %>% #keep only dates within index interval
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>% #order PCRs so that positive are first
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1) %>% #keep the first result, ensuring that positives within interval are kept before negatives
  anti_join(index_data4, by="Participant_ID") %>% #remove index cases from it 
  select(Household_ID, Participant_ID, PCR_result, Ct_val, Date)

# Positive primary data and primary dates
pos_primary_4_data <- filter(primary_4_data, PCR_result == "Positive") #from primary exposures, keep only those positive

primary_4_dates <- pos_primary_4_data %>%
  group_by(Household_ID) %>%
  summarize(primary4_max_date = max(Date)) %>% #identify last positive date amongst primaries within household
  mutate(
    primary_4_enddate = as.Date(primary4_max_date) + time_window,
    primary4_interval = interval(primary4_max_date, primary_4_enddate)
  ) #set up the interval from which household participants can be exposed to primary (non-index) cases

# Secondary data
secondary_4_data <- PCR_data %>%
  left_join(primary_4_dates, by="Household_ID") %>%
  filter(Date %within% primary4_interval & !is.na(Participant_ID)) %>% #identify participants with visits in the exposure interval of primaries
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>%
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1)  %>%#keep one result for each participant exposed to a primary, positive >negative
  select(Household_ID, Participant_ID, PCR_result, Ct_val, Date) 
#those exposed to a primary and positive are secondaries

# Positive secondary data and secondary dates
pos_secondary_4_data <- filter(secondary_4_data, PCR_result == "Positive") #keep only positive secondaries

secondary_4_dates <- pos_secondary_4_data %>%
  group_by(Household_ID) %>%
  summarize(secondary4_max_date = max(Date)) %>%
  mutate(
    secondary_4_enddate = as.Date(secondary4_max_date) + time_window,
    secondary4_interval = interval(secondary4_max_date, secondary_4_enddate)
  ) #identify latest positive date amongst secondaries and produce exposure interval

# Tertiary data
tertiary_4_data <- PCR_data %>%
  left_join(secondary_4_dates, by="Household_ID") %>%
  filter(Date %within% secondary4_interval & !is.na(Participant_ID)) %>%
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>%
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1) %>%#identify any tertiary cases, exposed to secondaries
  select(Household_ID, Participant_ID, PCR_result, Ct_val, Date)  

#remove unnecessary columns
prim_4_data <- primary_4_data
sec_4_data <- secondary_4_data
tert_4_data <- tertiary_4_data
ind_4_data <- index_data4

#Keep one row per exposed participant, keeping positive > Neg >NA
cluster_4_data <- prim_4_data %>%
  bind_rows(sec_4_data, tert_4_data) %>%
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative",NA))) %>% 
  group_by(Participant_ID) %>% 
  arrange(PCR_result, Date) %>%
  filter(row_number()==1) 

ind_4_data <- ind_4_data  %>% group_by(Household_ID) %>% mutate(count = n_distinct(Participant_ID)) #count number of index's per household 
single_index_HH_4 <- ind_4_data %>% 
  filter(count==1) %>%
  select(-count) %>%
  mutate(
    Index_ID = Participant_ID, 
    Date = first_case_date
  )

cluster_4_data_sing <- single_index_HH_4 %>%
  left_join(cluster_4_data, by="Household_ID") %>%
  select(-Participant_ID.x, -Ct_val.y) %>%
  rename(
    Participant_ID = Participant_ID.y,
    Index_Ct = Ct_val.x,
    Date = Date.y
  )
#results in single index clusters with index data included 


#separate positives to calculate end dates
all_pos <- pos_HH_data %>%
  group_by(Participant_ID) %>%
  arrange(Date) %>%
  mutate(num=1:n()) %>% 
  select(Participant_ID, Household_ID, Date)


#need to define an end date to the cluster
end_cluster4 <- cluster_4_data_sing %>%
  select(Household_ID, Participant_ID, PCR_result, Date) %>%
  filter(PCR_result == "Positive") %>%
  bind_rows(single_index_HH_4) %>%
  select(Household_ID, Participant_ID, Date) %>%
  left_join(all_pos, by=c("Participant_ID", "Household_ID")) %>%
  mutate(date_diff = difftime(Date.y, Date.x, units="days")) %>%
  filter(date_diff <=28) %>%
  group_by(Participant_ID) %>%
  mutate(cluster_end_date = max(Date.y)) %>%
  group_by(Household_ID) %>%
  summarize(cluster_end_date = max(cluster_end_date))


#Symptom data 
PCR_index_4 <- PCR_data %>%
  left_join(single_index_HH_4, by="Participant_ID") %>%
  filter(!is.na(first_case_date)) %>%
  mutate(date_diff = difftime(as.Date(Date.x), as.Date(first_case_date), units="days")) %>%
  filter(date_diff <10 & date_diff >=0) %>%
  filter(!is.na(Participant_ID)) %>%
  mutate(symp = case_when(
    is.na(`1 - Have you had ILI symptoms in the last week?`) ~ "Missing",
    `1 - Have you had ILI symptoms in the last week?` == "Yes" ~ "1",
    TRUE ~ "0"
  )) %>% arrange(factor(symp, levels=c(1, "Missing", 0))) %>%
  filter(row_number()==1) %>%
  select(Participant_ID, symp)

single_index_HH_4 <- single_index_HH_4 %>% 
  left_join(PCR_index_4, by="Participant_ID") %>%
  select(Household_ID, Index_ID, Ct_val, symp, Date) %>%
  rename(
    Index_date = Date
  )



cluster_4_data_sing <- single_index_HH_4 %>%
  left_join(cluster_4_data, by="Household_ID") %>%
  select(-Ct_val.y) %>%
  rename(
    Index_Ct = Ct_val.x
  ) %>%
  mutate(Cluster_no = 4)

cluster_4_data_sing_prim <- single_index_HH_4 %>%
  left_join(prim_4_data, by="Household_ID") %>%
  select(-Ct_val.y) %>%
  rename(
    Index_Ct = Ct_val.x
  ) %>%
  mutate(Cluster_no = 4)



# Cluster 5 ---------------------------------------------------------------

pos_HH_data <- pos_data %>%
  left_join(full_demog_data, by="Participant_ID") %>%
  mutate(Household_ID = as.factor(Household_ID))

PCR_data <- visit_data_fu %>%
  left_join(full_demog_data, by="Participant_ID") %>%
  mutate(Household_ID = as.factor(Household_ID))

# Calculate index data
index_data5 <- pos_data %>%
  left_join(full_demog_data, by="Participant_ID") %>% #join pos data to demog data
  mutate(Household_ID = as.factor(Household_ID))%>%
  left_join(end_cluster4, by="Household_ID") %>%
  filter(Date > cluster_end_date) %>%
  group_by(Household_ID) %>%
  summarize(first_case_date = min(Date)) %>% #identify the first case in the household
  left_join(pos_HH_data, by="Household_ID") %>% #join all other positive 
  filter(Date == first_case_date) %>% #select only participants whos date is equal to first case date
  mutate(
    first_end_date = as.Date(first_case_date) + time_window,
    index_interval = interval(first_case_date, first_end_date)
  )  %>% 
  select(Household_ID, first_case_date, Participant_ID, Ct_val, first_end_date, index_interval) 


# Primary data
primary_5_data <- PCR_data %>% #using all visit data 
  left_join(select(index_data5, Household_ID, index_interval), by="Household_ID") %>% #left join the household ID and index interval from index5
  filter(Date %within% index_interval & !is.na(Participant_ID)) %>% #keep only dates within index interval
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>% #order PCRs so that positive are first
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1) %>% #keep the first result, ensuring that positives within interval are kept before negatives
  anti_join(index_data5, by="Participant_ID") %>% #remove index cases from it 
  select(Household_ID, Participant_ID, PCR_result, Ct_val, Date)

# Positive primary data and primary dates
pos_primary_5_data <- filter(primary_5_data, PCR_result == "Positive") #from primary exposures, keep only those positive

primary_5_dates <- pos_primary_5_data %>%
  group_by(Household_ID) %>%
  summarize(primary5_max_date = max(Date)) %>% #identify last positive date amongst primaries within household
  mutate(
    primary_5_enddate = as.Date(primary5_max_date) + time_window,
    primary5_interval = interval(primary5_max_date, primary_5_enddate)
  ) #set up the interval from which household participants can be exposed to primary (non-index) cases

# Secondary data
secondary_5_data <- PCR_data %>%
  left_join(primary_5_dates, by="Household_ID") %>%
  filter(Date %within% primary5_interval & !is.na(Participant_ID)) %>% #identify participants with visits in the exposure interval of primaries
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>%
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1)  %>%#keep one result for each participant exposed to a primary, positive >negative
  select(Household_ID, Participant_ID, PCR_result, Ct_val, Date) 
#those exposed to a primary and positive are secondaries

# Positive secondary data and secondary dates
pos_secondary_5_data <- filter(secondary_5_data, PCR_result == "Positive") #keep only positive secondaries

secondary_5_dates <- pos_secondary_5_data %>%
  group_by(Household_ID) %>%
  summarize(secondary5_max_date = max(Date)) %>%
  mutate(
    secondary_5_enddate = as.Date(secondary5_max_date) + time_window,
    secondary5_interval = interval(secondary5_max_date, secondary_5_enddate)
  ) #identify latest positive date amongst secondaries and produce exposure interval

# Tertiary data
tertiary_5_data <- PCR_data %>%
  left_join(secondary_5_dates, by="Household_ID") %>%
  filter(Date %within% secondary5_interval & !is.na(Participant_ID)) %>%
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>%
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1) %>%#identify any tertiary cases, exposed to secondaries
  select(Household_ID, Participant_ID, PCR_result, Ct_val, Date)  

#remove unnecessary columns
prim_5_data <- primary_5_data
sec_5_data <- secondary_5_data
tert_5_data <- tertiary_5_data
ind_5_data <- index_data5

#Keep one row per exposed participant, keeping positive > Neg >NA
cluster_5_data <- prim_5_data %>%
  bind_rows(sec_5_data, tert_5_data) %>%
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative",NA))) %>% 
  group_by(Participant_ID) %>% 
  arrange(PCR_result, Date) %>%
  filter(row_number()==1) 

ind_5_data <- ind_5_data  %>% group_by(Household_ID) %>% mutate(count = n_distinct(Participant_ID)) #count number of index's per household 
single_index_HH_5 <- ind_5_data %>% 
  filter(count==1) %>%
  select(-count) %>%
  mutate(
    Index_ID = Participant_ID, 
    Date = first_case_date
  )

cluster_5_data_sing <- single_index_HH_5 %>%
  left_join(cluster_5_data, by="Household_ID") %>%
  select(-Participant_ID.x, -Ct_val.y) %>%
  rename(
    Participant_ID = Participant_ID.y,
    Index_Ct = Ct_val.x,
    Date = Date.y
  )
#results in single index clusters with index data included 


#separate positives to calculate end dates
all_pos <- pos_HH_data %>%
  group_by(Participant_ID) %>%
  arrange(Date) %>%
  mutate(num=1:n()) %>% 
  select(Participant_ID, Household_ID, Date)


#need to define an end date to the cluster
end_cluster5 <- cluster_5_data_sing %>%
  select(Household_ID, Participant_ID, PCR_result, Date) %>%
  filter(PCR_result == "Positive") %>%
  bind_rows(single_index_HH_5) %>%
  select(Household_ID, Participant_ID, Date) %>%
  left_join(all_pos, by=c("Participant_ID", "Household_ID")) %>%
  mutate(date_diff = difftime(Date.y, Date.x, units="days")) %>%
  filter(date_diff <=28) %>%
  group_by(Participant_ID) %>%
  mutate(cluster_end_date = max(Date.y)) %>%
  group_by(Household_ID) %>%
  summarize(cluster_end_date = max(cluster_end_date))


#Symptom data 
PCR_index_5 <- PCR_data %>%
  left_join(single_index_HH_5, by="Participant_ID") %>%
  filter(!is.na(first_case_date)) %>%
  mutate(date_diff = difftime(as.Date(Date.x), as.Date(first_case_date), units="days")) %>%
  filter(date_diff <10 & date_diff >=0) %>%
  filter(!is.na(Participant_ID)) %>%
  mutate(symp = case_when(
    is.na(`1 - Have you had ILI symptoms in the last week?`) ~ "Missing",
    `1 - Have you had ILI symptoms in the last week?` == "Yes" ~ "1",
    TRUE ~ "0"
  )) %>% arrange(factor(symp, levels=c(1, "Missing", 0))) %>%
  filter(row_number()==1) %>%
  select(Participant_ID, symp)

single_index_HH_5 <- single_index_HH_5 %>% 
  left_join(PCR_index_5, by="Participant_ID") %>%
  select(Household_ID, Index_ID, Ct_val, symp, Date)  %>%
  rename(
    Index_date = Date
  )



cluster_5_data_sing <- single_index_HH_5 %>%
  left_join(cluster_5_data, by="Household_ID") %>%
  select(-Ct_val.y) %>%
  rename(
    Index_Ct = Ct_val.x
  ) %>%
  mutate(Cluster_no = 5)

cluster_5_data_sing_prim <- single_index_HH_5 %>%
  left_join(prim_5_data, by="Household_ID") %>%
  select(-Ct_val.y) %>%
  rename(
    Index_Ct = Ct_val.x
  ) %>%
  mutate(Cluster_no = 5)



# Cluster 6 ---------------------------------------------------------------

pos_HH_data <- pos_data %>%
  left_join(full_demog_data, by="Participant_ID") %>%
  mutate(Household_ID = as.factor(Household_ID))

PCR_data <- visit_data_fu %>%
  left_join(full_demog_data, by="Participant_ID") %>%
  mutate(Household_ID = as.factor(Household_ID))

# Calculate index data
index_data6 <- pos_data %>%
  left_join(full_demog_data, by="Participant_ID") %>% #join pos data to demog data
  mutate(Household_ID = as.factor(Household_ID))%>%
  left_join(end_cluster5, by="Household_ID") %>%
  filter(Date > cluster_end_date) %>%
  group_by(Household_ID) %>%
  summarize(first_case_date = min(Date)) %>% #identify the first case in the household
  left_join(pos_HH_data, by="Household_ID") %>% #join all other positive 
  filter(Date == first_case_date) %>% #select only participants whos date is equal to first case date
  mutate(
    first_end_date = as.Date(first_case_date) + time_window,
    index_interval = interval(first_case_date, first_end_date)
  )  %>% 
  select(Household_ID, first_case_date, Participant_ID, Ct_val, first_end_date, index_interval) 


# Primary data
primary_6_data <- PCR_data %>% #using all visit data 
  left_join(select(index_data6, Household_ID, index_interval), by="Household_ID") %>% #left join the household ID and index interval from index6
  filter(Date %within% index_interval & !is.na(Participant_ID)) %>% #keep only dates within index interval
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>% #order PCRs so that positive are first
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1) %>% #keep the first result, ensuring that positives within interval are kept before negatives
  anti_join(index_data6, by="Participant_ID") %>% #remove index cases from it 
  select(Household_ID, Participant_ID, PCR_result, Ct_val, Date)

# Positive primary data and primary dates
pos_primary_6_data <- filter(primary_6_data, PCR_result == "Positive") #from primary exposures, keep only those positive

primary_6_dates <- pos_primary_6_data %>%
  group_by(Household_ID) %>%
  summarize(primary6_max_date = max(Date)) %>% #identify last positive date amongst primaries within household
  mutate(
    primary_6_enddate = as.Date(primary6_max_date) + time_window,
    primary6_interval = interval(primary6_max_date, primary_6_enddate)
  ) #set up the interval from which household participants can be exposed to primary (non-index) cases

# Secondary data
secondary_6_data <- PCR_data %>%
  left_join(primary_6_dates, by="Household_ID") %>%
  filter(Date %within% primary6_interval & !is.na(Participant_ID)) %>% #identify participants with visits in the exposure interval of primaries
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>%
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1)  %>%#keep one result for each participant exposed to a primary, positive >negative
  select(Household_ID, Participant_ID, PCR_result, Ct_val, Date) 
#those exposed to a primary and positive are secondaries

# Positive secondary data and secondary dates
pos_secondary_6_data <- filter(secondary_6_data, PCR_result == "Positive") #keep only positive secondaries

secondary_6_dates <- pos_secondary_6_data %>%
  group_by(Household_ID) %>%
  summarize(secondary6_max_date = max(Date)) %>%
  mutate(
    secondary_6_enddate = as.Date(secondary6_max_date) + time_window,
    secondary6_interval = interval(secondary6_max_date, secondary_6_enddate)
  ) #identify latest positive date amongst secondaries and produce exposure interval

# Tertiary data
tertiary_6_data <- PCR_data %>%
  left_join(secondary_6_dates, by="Household_ID") %>%
  filter(Date %within% secondary6_interval & !is.na(Participant_ID)) %>%
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>%
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1) %>%#identify any tertiary cases, exposed to secondaries
  select(Household_ID, Participant_ID, PCR_result, Ct_val, Date)  

#remove unnecessary columns
prim_6_data <- primary_6_data
sec_6_data <- secondary_6_data
tert_6_data <- tertiary_6_data
ind_6_data <- index_data6

#Keep one row per exposed participant, keeping positive > Neg >NA
cluster_6_data <- prim_6_data %>%
  bind_rows(sec_6_data, tert_6_data) %>%
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative",NA))) %>% 
  group_by(Participant_ID) %>% 
  arrange(PCR_result, Date) %>%
  filter(row_number()==1) 

ind_6_data <- ind_6_data  %>% group_by(Household_ID) %>% mutate(count = n_distinct(Participant_ID)) #count number of index's per household 
single_index_HH_6 <- ind_6_data %>% 
  filter(count==1) %>%
  select(-count) %>%
  mutate(
    Index_ID = Participant_ID, 
    Date = first_case_date
  )

cluster_6_data_sing <- single_index_HH_6 %>%
  left_join(cluster_6_data, by="Household_ID") %>%
  select(-Participant_ID.x, -Ct_val.y) %>%
  rename(
    Participant_ID = Participant_ID.y,
    Index_Ct = Ct_val.x,
    Date = Date.y
  )
#results in single index clusters with index data included 


#separate positives to calculate end dates
all_pos <- pos_HH_data %>%
  group_by(Participant_ID) %>%
  arrange(Date) %>%
  mutate(num=1:n()) %>% 
  select(Participant_ID, Household_ID, Date)


#need to define an end date to the cluster
end_cluster6 <- cluster_6_data_sing %>%
  select(Household_ID, Participant_ID, PCR_result, Date) %>%
  filter(PCR_result == "Positive") %>%
  bind_rows(single_index_HH_6) %>%
  select(Household_ID, Participant_ID, Date) %>%
  left_join(all_pos, by=c("Participant_ID", "Household_ID")) %>%
  mutate(date_diff = difftime(Date.y, Date.x, units="days")) %>%
  filter(date_diff <=28) %>%
  group_by(Participant_ID) %>%
  mutate(cluster_end_date = max(Date.y)) %>%
  group_by(Household_ID) %>%
  summarize(cluster_end_date = max(cluster_end_date))


#Symptom data 
PCR_index_6 <- PCR_data %>%
  left_join(single_index_HH_6, by="Participant_ID") %>%
  filter(!is.na(first_case_date)) %>%
  mutate(date_diff = difftime(as.Date(Date.x), as.Date(first_case_date), units="days")) %>%
  filter(date_diff <10 & date_diff >=0) %>%
  filter(!is.na(Participant_ID)) %>%
  mutate(symp = case_when(
    is.na(`1 - Have you had ILI symptoms in the last week?`) ~ "Missing",
    `1 - Have you had ILI symptoms in the last week?` == "Yes" ~ "1",
    TRUE ~ "0"
  )) %>% arrange(factor(symp, levels=c(1, "Missing", 0))) %>%
  filter(row_number()==1) %>%
  select(Participant_ID, symp)

single_index_HH_6 <- single_index_HH_6 %>% 
  left_join(PCR_index_6, by="Participant_ID") %>%
  select(Household_ID, Index_ID, Ct_val, symp, Date) %>%
  rename(
    Index_date = Date
  )



cluster_6_data_sing <- single_index_HH_6 %>%
  left_join(cluster_6_data, by="Household_ID") %>%
  select(-Ct_val.y) %>%
  rename(
    Index_Ct = Ct_val.x
  ) %>%
  mutate(Cluster_no = 6)


cluster_6_data_sing_prim <- single_index_HH_6 %>%
  left_join(prim_6_data, by="Household_ID") %>%
  select(-Ct_val.y) %>%
  rename(
    Index_Ct = Ct_val.x
  ) %>%
  mutate(Cluster_no = 6)




# ANALYSIS ----------------------------------------------------------------

#Create cluster_no as variable
ind_1_data$cluster_no <- 1
ind_2_data$cluster_no <- 2
ind_3_data$cluster_no <- 3
ind_4_data$cluster_no <- 4
ind_5_data$cluster_no <- 5

#Create dataframe of all index cases
total_index_data <- bind_rows(
  ind_1_data, ind_2_data, 
  ind_3_data, ind_4_data, 
  ind_5_data
)
#Create Cluster ID
total_index_data <- total_index_data %>%
  group_by(Household_ID, first_case_date) %>%
  mutate(cluster_ID = cur_group_id()) %>%
  ungroup()

#Identify whether clusters were single or multi index
total_index_data <- total_index_data %>%
  group_by(cluster_ID) %>%
  filter(row_number()==1) %>%
  mutate(S_Ind = case_when(
    count == 1 ~ "1", 
    count != 1 ~ "0",
    TRUE ~ NA
  ))

summary(total_index_data$cluster_no) #Number of clusters per household
table(total_index_data$S_Ind) #number of single vs multi index
table(total_index_data$count, total_index_data$S_Ind) #split by number of indexes

#Single Index only dataframe
total_single_index_data <- single_index_HH_1 %>% bind_rows(
  single_index_HH_2, 
  single_index_HH_3, single_index_HH_4, 
  single_index_HH_5
) %>% rename(
  Participant_ID = Index_ID
) %>% 
  left_join(full_demog_data %>% select(Participant_ID, age_cat), by="Participant_ID") %>%
  rename(
    Index_ID = Participant_ID
  )

#Produce single cluster dataset for analysis
total_single_cluster_data <- cluster_1_data_sing %>% 
  bind_rows(
    cluster_2_data_sing, cluster_3_data_sing, cluster_4_data_sing, cluster_5_data_sing
  ) %>%
  group_by(Household_ID, Cluster_no) %>%
  mutate(cluster_ID = cur_group_id()) %>%
  ungroup() %>%
  mutate(period=case_when( #incorporate period as a variable
    Index_date <= "2021-07-07" ~ "Pre-Delta", 
    Index_date > "2021-07-07" & Index_date <= "2021-12-04" ~ "Delta", 
    Index_date > "2021-12-04" ~ "Omicron", 
    TRUE ~ NA
  ))



cluster_data_prim <- cluster_1_data_sing_prim %>%
  bind_rows(
    cluster_2_data_sing_prim, 
    cluster_3_data_sing_prim, 
    cluster_4_data_sing_prim, 
    cluster_5_data_sing_prim
  ) %>%
  group_by(Household_ID, Cluster_no) %>%
  mutate(cluster_ID = cur_group_id()) %>%
  ungroup() %>%
  mutate(period=case_when( #incorporate period as a variable
    Index_date <= "2021-07-07" ~ "Pre-Delta", 
    Index_date > "2021-07-07" & Index_date <= "2021-12-04" ~ "Delta", 
    Index_date > "2021-12-04" ~ "Omicron", 
    TRUE ~ NA
  ))


#Addition of serology data for index cases

#import all_visit_serol

index_serol <- All_visit_serol
colnames(index_serol)[1] <- "Index_ID"

contact_serol <- All_visit_serol
colnames(contact_serol)[1] <- "Participant_ID"

cluster_data_prim <- cluster_data_prim %>%
  left_join(index_serol, by="Index_ID") %>%
  mutate(sero = case_when(
    V1_spike == "Pos" ~ "Pos",
    Index_date < V2_date & V1_spike == "Neg" ~ "Neg",
    Index_date > V2_date & V2_spike == "Neg" ~ "Neg", 
    Index_date > V2_date & V2_spike =="Pos" ~ "Pos", 
    TRUE ~ NA
  )) %>% 
  select(-V1_date, -V1_spike, -V1_ncp, 
         -V2_date, -V2_spike, -V2_ncp, 
         -V3_date, -V3_spike, -V3_ncp) %>%
  left_join(contact_serol, by="Participant_ID") %>%
  mutate(sero_contact = case_when(
    V1_spike == "Pos" ~ "Pos",
    Index_date < V2_date & V1_spike == "Neg" ~ "Neg",
    Index_date > V2_date & V2_spike == "Neg" ~ "Neg", 
    Index_date > V2_date & V2_spike =="Pos" ~ "Pos", 
    TRUE ~ NA
  )) %>% 
  select(-V1_date, -V1_spike, -V1_ncp, 
         -V2_date, -V2_spike, -V2_ncp, 
         -V3_date, -V3_spike, -V3_ncp)

total_single_cluster_data <- total_single_cluster_data %>%
  left_join(index_serol, by="Index_ID") %>%
  mutate(sero = case_when(
    V1_spike == "Pos" ~ "Pos",
    Index_date < V2_date & V1_spike == "Neg" ~ "Neg",
    Index_date > V2_date & V2_spike == "Neg" ~ "Neg", 
    Index_date > V2_date & V2_spike =="Pos" ~ "Pos", 
    TRUE ~ NA
  )) %>% 
  select(-V1_date, -V1_spike, -V1_ncp, 
         -V2_date, -V2_spike, -V2_ncp, 
         -V3_date, -V3_spike, -V3_ncp) %>%
  left_join(contact_serol, by="Participant_ID") %>%
  mutate(sero_contact = case_when(
    V1_spike == "Pos" ~ "Pos",
    Index_date < V2_date & V1_spike == "Neg" ~ "Neg",
    Index_date > V2_date & V2_spike == "Neg" ~ "Neg", 
    Index_date > V2_date & V2_spike =="Pos" ~ "Pos", 
    TRUE ~ NA
  )) %>% 
  select(-V1_date, -V1_spike, -V1_ncp, 
         -V2_date, -V2_spike, -V2_ncp, 
         -V3_date, -V3_spike, -V3_ncp)
  
