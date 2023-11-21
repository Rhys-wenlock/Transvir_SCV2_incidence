library(dplyr)
library(lubridate)

time_window <- 14

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
  ) #create index time interval for exposure

#End result = participants whos first positive date = households first positive date, with an interval of exposure

# Primary data
primary_1_data <- PCR_data %>% #using all visit data 
  left_join(select(index_data1, Household_ID, index_interval), by="Household_ID") %>% #left join the household ID and index interval from index1
  filter(Date %within% index_interval & !is.na(Participant_ID)) %>% #keep only dates within index interval
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>% #order PCRs so that positive are first
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1) %>% #keep the first result, ensuring that positives within interval are kept before negatives
  anti_join(index_data1, by="Participant_ID") #remove index cases from it

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
  filter(row_number() == 1) #keep one result for each participant exposed to a primary, positive >negative
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
  filter(row_number() == 1) #identify any tertiary cases, exposed to secondaries

#remove unnecessary columns
prim_1_data <- primary_1_data[,-c(43:45)]
sec_1_data <- secondary_1_data[,-c(45:47)]
tert_1_data <- tertiary_1_data[,-c(45:47)]
ind_1_data <- index_data1[,-c(45:47)]

#bind exposed participant data together
cluster_1_data <- rbind(prim_1_data, sec_1_data, tert_1_data)
cluster_1_data$PCR_result <- factor(cluster_1_data$PCR_result, levels=c("Positive", "Negative",NA))
cluster_1_data <- cluster_1_data %>% group_by(Participant_ID) %>% arrange(cluster_1_data$PCR_result, cluster_1_data$Date) %>% filter(row_number()==1)
#keep one row per exposed participant, with their outcome arranged by positive -> negative, with earliest first

ind_1_data <- ind_1_data  %>% group_by(Household_ID) %>% mutate(count = n_distinct(Participant_ID)) #count number of index's per household 
single_index_HH_1 <- ind_1_data[ind_1_data$count==1,-21] #select individiual index clusters
cluster_1_data <- anti_join(cluster_1_data, single_index_HH_1, by="Participant_ID")
cluster_1_data <- left_join(single_index_HH_1, cluster_1_data, by="Household_ID")
#results in single index clusters


#separate positives to calculate end dates
positives_cluster1 <- cluster_1_data[cluster_1_data$PCR_result=="Positive",] #keep just positives 
positives_cluster1 <- rbind(positives_cluster1, ind_1_data) #join positives with indexes (in case there were no positives after exposure)

#need to define an end date to the cluster
positives_cluster1 <- positives_cluster1[,c(24,1,4)]
colnames(positives_cluster1) <- c("Household_ID","Participant_ID",  "Date")
test <-  PCR_data[PCR_data$PCR_result=="Positive",] %>% group_by(Participant_ID) #group by ID
test <- test %>% arrange(Date, .by_group = TRUE) #organise dates into order per group
test <- test %>% mutate(num = 1:n()) #add count
test$Household_ID <- as.factor(test$Household_ID)
end_cluster1 <- left_join(positives_cluster1, test, by=c("Participant_ID", "Household_ID"))
#so now num = 1 for first positive, then increasing number count for each subsequent positive PCR

#calcualte difference between subsequent positives and the first (in entire study)
end_cluster1$date_diff <- difftime(end_cluster1$`Date.y`, end_cluster1$`Date.x`, units="days")

#for participants with a repeat positive within 28 days 
end_cluster1 <- end_cluster1[end_cluster1$date_diff<=28,] %>% group_by(Participant_ID) %>% mutate(max_date <- max(`Date.y`))
end_cluster1 <- end_cluster1[,c(1,2,48)]
colnames(end_cluster1)[3] <- "cluster_end_date"
end_cluster1 <- end_cluster1 %>% group_by(Household_ID) %>% summarize(max(cluster_end_date))
colnames(end_cluster1)[2] <- "cluster_end_date"
end_cluster1$Household_ID <- as.factor(end_cluster1$Household_ID)
end_cluster1$cluster_end_date <- as.Date(end_cluster1$cluster_end_date) 



#symptom data
colnames(single_index_HH_1)[4] <- "Index_date"
PCR_index_1 <- left_join(PCR_data, single_index_HH_1, by="Participant_ID")
PCR_index_1$date_diff <- difftime(PCR_index_1$Date, PCR_index_1$Index_date, units="days") #create time between first infection and all visis within 10 days of it 

PCR_index_1 <- PCR_index_1[PCR_index_1$date_diff<10 & PCR_index_1$date_diff>=0,] #include only visits within 10 days of third positive test
PCR_index_1 <- PCR_index_1[!is.na(PCR_index_1$Participant_ID),] #remove NAs
PCR_index_1$symp <- ifelse(is.na(PCR_index_1$`1 - Have you had ILI symptoms in the last week?.x`), "Missing", 
                           ifelse(PCR_index_1$`1 - Have you had ILI symptoms in the last week?.x`=="Yes", 1,0))
PCR_index_1 <- PCR_index_1 %>% arrange(factor(PCR_index_1$symp, levels = c(1, "Missing",0))) #order varibales so that Yes is kept over NA if people had two visits within 10 days 
PCR_index_1 <- PCR_index_1 %>% distinct(PCR_index_1$Participant_ID, .keep_all = TRUE) #keep only one episode per participant ID
PCR_index_1_symp <- PCR_index_1[,c(1,24,50)]

single_index_HH_1 <- left_join(single_index_HH_1, PCR_index_1_symp, by="Participant_ID")



# Cluster 2 ---------------------------------------------------------------

pos_HH_data <- pos_data %>%
  left_join(full_demog_data, by="Participant_ID") %>%
  mutate(Household_ID = as.factor(Household_ID))

PCR_data <- visit_data_fu %>%
  left_join(full_demog_data, by="Participant_ID") %>%
  mutate(Household_ID = as.factor(Household_ID))

# Calculate index data
index_data2 <- pos_data %>%
  left_join(full_demog_data, by="Participant_ID") %>%
  mutate(Household_ID = as.factor(Household_ID)) %>%
  left_join(end_cluster1, by="Household_ID") %>%
  filter(Date > cluster_end_date) %>%
  group_by(Household_ID) %>%
  summarize(first_case_date = min(Date)) %>%
  left_join(pos_HH_data, by="Household_ID") %>%
  filter(Date == first_case_date) %>%
  mutate(
    first_end_date = as.Date(first_case_date) + time_window,
    index_interval = interval(first_case_date, first_end_date)
  ) %>%
  select(Household_ID, Participant_ID, Date, Ct_val, 9, age_cat, Sex, index_interval) %>%
  rename(
     Symp =`1 - Have you had ILI symptoms in the last week?`
  )

# Primary data
primary_2_data <- PCR_data %>%
  left_join(select(index_data2, Household_ID, index_interval), by="Household_ID") %>%
  filter(Date %within% index_interval & !is.na(Participant_ID)) %>%
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>%
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1) %>%
  anti_join(index_data2, by="Participant_ID") %>% 
  select(Participant_ID, Date, PCR_result, Household_ID)

# Positive primary data and primary dates
pos_primary_2_data <- filter(primary_2_data, PCR_result == "Positive")

primary_2_dates <- pos_primary_2_data %>%
  group_by(Household_ID) %>%
  summarize(primary2_max_date = max(Date)) %>%
  mutate(
    primary_2_enddate = as.Date(primary2_max_date) + time_window,
    primary2_interval = interval(primary2_max_date, primary_2_enddate)
  )

# Secondary data
secondary_2_data <- PCR_data %>%
  left_join(primary_2_dates, by="Household_ID") %>%
  filter(Date %within% primary2_interval & !is.na(Participant_ID)) %>%
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>%
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1) %>% 
  select(Participant_ID, Date, PCR_result, Household_ID)

# Positive secondary data and secondary dates
pos_secondary_2_data <- filter(secondary_2_data, PCR_result == "Positive")

secondary_2_dates <- pos_secondary_2_data %>%
  group_by(Household_ID) %>%
  summarize(secondary2_max_date = max(Date)) %>%
  mutate(
    secondary_2_enddate = as.Date(secondary2_max_date) + time_window,
    secondary2_interval = interval(secondary2_max_date, secondary_2_enddate)
  )

# Tertiary data
tertiary_2_data <- PCR_data %>%
  left_join(secondary_2_dates, by="Household_ID") %>%
  filter(Date %within% secondary2_interval & !is.na(Participant_ID)) %>%
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>%
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1) %>% 
  select(Participant_ID, Date, PCR_result, Household_ID)

ind_2_data <- index_data2

cluster_2_data <- rbind(primary_2_data, secondary_2_data, tertiary_2_data)
cluster_2_data$PCR_result <- factor(cluster_2_data$PCR_result, levels=c("Positive", "Negative",NA))
cluster_2_data <- cluster_2_data %>% group_by(Participant_ID) %>% arrange(cluster_2_data$PCR_result, cluster_2_data$Date) %>% filter(row_number()==1)

ind_2_data <- ind_2_data  %>% group_by(Household_ID) %>% mutate(count = n_distinct(Participant_ID))
single_index_HH_2 <- ind_2_data[ind_2_data$count==1,-9] 
#separate positives to calculate end dates
positives_cluster2 <- cluster_2_data[cluster_2_data$PCR_result=="Positive",]
positives_cluster2 <- rbind(positives_cluster2, ind_2_data)
cluster_2_data <- anti_join(cluster_2_data, single_index_HH_2, by="Participant_ID")
cluster_2_data <- left_join(single_index_HH_2, cluster_2_data, by="Household_ID")
#results in single index clusters

#need to define an end date to the cluster
positives_cluster2 <- positives_cluster2 %>%
  select(Household_ID, Participant_ID, Date)
test <-  PCR_data[PCR_data$PCR_result=="Positive",] %>% group_by(Participant_ID) #group by ID
test <- test %>% arrange(Date, .by_group = TRUE) #organise dates into order per group
test <- test %>% mutate(num = 1:n()) #add count
test$Household_ID <- as.factor(test$Household_ID)
end_cluster2 <- left_join(positives_cluster2, test, by=c("Participant_ID", "Household_ID"))
#so now num = 1 for first positive, then increasing number count for each subsequent positive PCR

#calcualte difference between subsequent positives and the first (in entire study)
end_cluster2$date_diff <- difftime(end_cluster2$`Date.y`, end_cluster2$`Date.x`, units="days")

#for participants with a repeat positive within 28 days 
end_cluster2 <- end_cluster2[end_cluster2$date_diff<=28,] %>% group_by(Participant_ID) %>% mutate(max_date <- max(`Date.y`))
end_cluster2 <- end_cluster2[,c(1,2,48)]
colnames(end_cluster2)[3] <- "cluster_end_date"
end_cluster2 <- end_cluster2 %>% group_by(Household_ID) %>% summarize(max(cluster_end_date))
colnames(end_cluster2)[2] <- "cluster_end_date"



#symptom data
colnames(single_index_HH_2)[3] <- "Index_date"
PCR_index_2 <- left_join(PCR_data, single_index_HH_2, by="Participant_ID")
PCR_index_2$date_diff <- difftime(PCR_index_2$Date, PCR_index_2$Index_date, units="days") #create time between first infection and all visis within 10 days of it 

PCR_index_2 <- PCR_index_2[PCR_index_2$date_diff<10 & PCR_index_2$date_diff>=0,] #include only visits within 10 days of third positive test
PCR_index_2 <- PCR_index_2[!is.na(PCR_index_2$Participant_ID),] #remove NAs
PCR_index_2$symp <- ifelse(is.na(PCR_index_2$`1 - Have you had ILI symptoms in the last week?`), "Missing", 
                           ifelse(PCR_index_2$`1 - Have you had ILI symptoms in the last week?`=="Yes", 1,0))
PCR_index_2 <- PCR_index_2 %>% arrange(factor(PCR_index_2$symp, levels = c(1, "Missing",0))) #order varibales so that Yes is kept over NA if people had two visits within 10 days 
PCR_index_2 <- PCR_index_2 %>% distinct(PCR_index_2$Participant_ID, .keep_all = TRUE) #keep only one episode per participant ID
PCR_index_2_symp <- PCR_index_2[,c(1,24,48)]

single_index_HH_2 <- left_join(single_index_HH_2, PCR_index_2_symp, by="Participant_ID")




# Cluster 3 ---------------------------------------------------------------

pos_HH_data <- pos_data %>%
  left_join(full_demog_data, by="Participant_ID") %>%
  mutate(Household_ID = as.factor(Household_ID))

PCR_data <- visit_data_fu %>%
  left_join(full_demog_data, by="Participant_ID") %>%
  mutate(Household_ID = as.factor(Household_ID))

# Calculate index data
index_data3 <- pos_data %>%
  left_join(full_demog_data, by="Participant_ID") %>%
  mutate(Household_ID = as.factor(Household_ID)) %>%
  left_join(end_cluster2, by="Household_ID") %>%
  filter(Date > cluster_end_date) %>%
  group_by(Household_ID) %>%
  summarize(first_case_date = min(Date)) %>%
  left_join(pos_HH_data, by="Household_ID") %>%
  filter(Date == first_case_date) %>%
  mutate(
    first_end_date = as.Date(first_case_date) + time_window,
    index_interval = interval(first_case_date, first_end_date)
  ) %>%
  select(Household_ID, Participant_ID, Date, Ct_val, 9, age_cat, Sex, index_interval) %>%
  rename(
    Symp =`1 - Have you had ILI symptoms in the last week?`
  )

# Primary data
primary_3_data <- PCR_data %>%
  left_join(select(index_data2, Household_ID, index_interval), by="Household_ID") %>%
  filter(Date %within% index_interval & !is.na(Participant_ID)) %>%
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>%
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1) %>%
  anti_join(index_data3, by="Participant_ID") %>% 
  select(Participant_ID, Date, PCR_result, Household_ID)

# Positive primary data and primary dates
pos_primary_3_data <- filter(primary_3_data, PCR_result == "Positive")

primary_3_dates <- pos_primary_3_data %>%
  group_by(Household_ID) %>%
  summarize(primary3_max_date = max(Date)) %>%
  mutate(
    primary_3_enddate = as.Date(primary3_max_date) + time_window,
    primary3_interval = interval(primary3_max_date, primary_3_enddate)
  )

# Secondary data
secondary_3_data <- PCR_data %>%
  left_join(primary_3_dates, by="Household_ID") %>%
  filter(Date %within% primary3_interval & !is.na(Participant_ID)) %>%
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>%
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1) %>% 
  select(Participant_ID, Date, PCR_result, Household_ID)

# Positive secondary data and secondary dates
pos_secondary_3_data <- filter(secondary_3_data, PCR_result == "Positive")

secondary_3_dates <- pos_secondary_3_data %>%
  group_by(Household_ID) %>%
  summarize(secondary3_max_date = max(Date)) %>%
  mutate(
    secondary_3_enddate = as.Date(secondary3_max_date) + time_window,
    secondary3_interval = interval(secondary3_max_date, secondary_3_enddate)
  )

# Tertiary data
tertiary_3_data <- PCR_data %>%
  left_join(secondary_3_dates, by="Household_ID") %>%
  filter(Date %within% secondary3_interval & !is.na(Participant_ID)) %>%
  mutate(PCR_result = factor(PCR_result, levels=c("Positive", "Negative", NA))) %>%
  group_by(Participant_ID) %>%
  arrange(PCR_result, Date) %>%
  filter(row_number() == 1) %>% 
  select(Participant_ID, Date, PCR_result, Household_ID)

ind_3_data <- index_data3

cluster_3_data <- rbind(primary_3_data, secondary_3_data, tertiary_3_data)
cluster_3_data$PCR_result <- factor(cluster_3_data$PCR_result, levels=c("Positive", "Negative",NA))
cluster_3_data <- cluster_3_data %>% group_by(Participant_ID) %>% arrange(cluster_3_data$PCR_result, cluster_3_data$Date) %>% filter(row_number()==1)

ind_3_data <- ind_3_data  %>% group_by(Household_ID) %>% mutate(count = n_distinct(Participant_ID))
single_index_HH_3 <- ind_3_data[ind_3_data$count==1,-9] 
#separate positives to calculate end dates
positives_cluster3 <- cluster_3_data[cluster_3_data$PCR_result=="Positive",]
positives_cluster3 <- rbind(positives_cluster3, ind_3_data)
cluster_3_data <- anti_join(cluster_3_data, single_index_HH_3, by="Participant_ID")
cluster_3_data <- left_join(single_index_HH_3, cluster_3_data, by="Household_ID")
#results in single index clusters

#need to define an end date to the cluster
positives_cluster3 <- positives_cluster3 %>%
  select(Household_ID, Participant_ID, Date)
test <-  PCR_data[PCR_data$PCR_result=="Positive",] %>% group_by(Participant_ID) #group by ID
test <- test %>% arrange(Date, .by_group = TRUE) #organise dates into order per group
test <- test %>% mutate(num = 1:n()) #add count
test$Household_ID <- as.factor(test$Household_ID)
end_cluster3 <- left_join(positives_cluster3, test, by=c("Participant_ID", "Household_ID"))
#so now num = 1 for first positive, then increasing number count for each subsequent positive PCR

#calcualte difference between subsequent positives and the first (in entire study)
end_cluster3$date_diff <- difftime(end_cluster3$`Date.y`, end_cluster3$`Date.x`, units="days")

#for participants with a repeat positive within 28 days 
end_cluster3 <- end_cluster3[end_cluster3$date_diff<=28,] %>% group_by(Participant_ID) %>% mutate(max_date <- max(`Date.y`))
end_cluster3 <- end_cluster3[,c(1,2,48)]
colnames(end_cluster3)[3] <- "cluster_end_date"
end_cluster3 <- end_cluster3 %>% group_by(Household_ID) %>% summarize(max(cluster_end_date))
colnames(end_cluster3)[2] <- "cluster_end_date"



#symptom data
colnames(single_index_HH_3)[3] <- "Index_date"
PCR_index_3 <- left_join(PCR_data, single_index_HH_3, by="Participant_ID")
PCR_index_3$date_diff <- difftime(PCR_index_3$Date, PCR_index_3$Index_date, units="days") #create time between first infection and all visis within 10 days of it 

PCR_index_3 <- PCR_index_3[PCR_index_3$date_diff<10 & PCR_index_3$date_diff>=0,] #include only visits within 10 days of third positive test
PCR_index_3 <- PCR_index_3[!is.na(PCR_index_3$Participant_ID),] #remove NAs
PCR_index_3$symp <- ifelse(is.na(PCR_index_3$`1 - Have you had ILI symptoms in the last week?`), "Missing", 
                           ifelse(PCR_index_3$`1 - Have you had ILI symptoms in the last week?`=="Yes", 1,0))
PCR_index_3 <- PCR_index_3 %>% arrange(factor(PCR_index_3$symp, levels = c(1, "Missing",0))) #order varibales so that Yes is kept over NA if people had two visits within 10 days 
PCR_index_3 <- PCR_index_3 %>% distinct(PCR_index_3$Participant_ID, .keep_all = TRUE) #keep only one episode per participant ID
PCR_index_3_symp <- PCR_index_3[,c(1,24,48)]

single_index_HH_3 <- left_join(single_index_HH_3, PCR_index_3_symp, by="Participant_ID")



