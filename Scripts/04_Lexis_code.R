require(lubridate)
require(tidyverse)
require(Epi)
require(popEpi)

#import the correct dataset - depending on sensitivity analysis to be conducted
#default is base
all_inf_data <- read.csv(here::here("Output", "all_inf_data.csv"))
full_demog_data <- read.csv(here::here("Output", "full_demog_data.csv"))

#Base = all_inf_data_base.csv - median SN and median N vacc
#1 = all_inf_data_1.csv - IQR1 SN and median N vaccinated
#2 = all_inf_data2.csv - IQR1 SN and IQR1 N vaccinated
#3 = all_inf_data3.csv - median S and median N vaccinated
#4 = all_inf_data4.csv - median N and median N vaccinated
#5 = all_inf_data5.csv - threshold S and threshold N vaccinated
#6 = all_inf_data6.csv - threshold N and threshold N vaccinated
#7 = all_inf_data7.csv - median SN and median N vaccinated - 90 day reinfection
#8 = all_inf_data8.csv - PCR+only



#Produce Lexis Object 
lexis_inf <- all_inf_data %>% select(-...1) %>%
  mutate(
    timein = cal.yr(as.Date(start_date, format = "%Y-%m-%d")),
    timeout = cal.yr(as.Date(end_date, format = "%Y-%m-%d")),
    V2_date = cal.yr(as.Date(V2_date, format = "%Y-%m-%d")),
    Baseline_date = cal.yr(as.Date(V1_date, format = "%Y-%m-%d")),
    infection = if_else(!is.na(Inf_date),1,0),
    ID = `Participant_ID`
  ) %>%
  {Lexis(entry = list("calendar_time"=.$timein), 
         exit = list("calendar_time" = .$timeout),
         exit.status = .$infection ==1, id = .$ID, data = ., keep.dropped=TRUE)}


rate(lexis_inf, obs= lex.Xst, pyrs=lex.dur)

full_demog_data <- full_demog_data[,-c(27:35)]

lexis_inf_split <- lexis_inf %>%
  cutLexis(lexis_inf$V2_date, timescale = "calendar_time") %>%
  splitLexis(breaks=c(2021.513, 2021.923)) %>%
  left_join(full_demog_data, by = "Participant_ID") %>%
  mutate(
    household_children = Children_5_17 + Children_2_4 + Children_6m_2 + Children_under_6m,
    crowding = Total_household_no / Total_rooms,
    hh_children = cut(household_children, c(0, 4, 5, Inf), right = FALSE),
    hh_size = cut(Total_household_no, c(0, 5, 7, 10, Inf), right = FALSE),
    sero = if_else(calendar_time < V2_date | is.na(V2_date), V1_spike, V2_spike),
    period = case_when(
      calendar_time < 2021.513 ~ "Pre-Delta",
      calendar_time < 2021.923 ~ "Delta",
      TRUE ~ "Omicron"
    ),
    V1_sero = if_else(V1_spike == "Pos"|is.na(V1_spike), 1, 0),
    V2_sero2 = if_else(V2_spike == "Pos", 1, 0),
    First_positive_date = cal.yr(as.Date(First_positive_date, format = "%Y-%m-%d")),
    Second_positive_date = cal.yr(as.Date(second_positive_date, format = "%Y-%m-%d")),
    third_positive_date = cal.yr(as.Date(third_positive_date, format = "%Y-%m-%d"))
  ) %>%
  mutate(
    prior_infection = case_when(
      V1_sero == 0 & total_inf == 0 ~ 0,
      V1_sero == 1 & total_inf == 0 ~ 1,
      V1_sero == 0 & calendar_time < First_positive_date ~ 0,
      V1_sero == 1 & calendar_time < First_positive_date ~ 1,
      V1_sero == 0 & calendar_time > First_positive_date & calendar_time < Second_positive_date ~ 1,
      V1_sero == 1 & calendar_time > First_positive_date & calendar_time < Second_positive_date ~ 2,
      V1_sero == 0 & calendar_time > First_positive_date & is.na(Second_positive_date) ~ 1,
      V1_sero == 1 & calendar_time > First_positive_date & is.na(Second_positive_date) ~ 2,
      V1_sero == 0 & calendar_time > First_positive_date & calendar_time > Second_positive_date & calendar_time < third_positive_date ~ 2,
      V1_sero == 1 & calendar_time > First_positive_date & calendar_time > Second_positive_date & calendar_time < third_positive_date ~ 3,
      V1_sero == 0 & calendar_time > First_positive_date & calendar_time > Second_positive_date & is.na(third_positive_date) ~ 2,
      V1_sero == 1 & calendar_time > First_positive_date & calendar_time > Second_positive_date & is.na(third_positive_date) ~ 3,
      V1_sero == 0 & calendar_time > First_positive_date & calendar_time > Second_positive_date & calendar_time > third_positive_date ~ 3,
      V1_sero == 1 & calendar_time > First_positive_date & calendar_time > Second_positive_date & calendar_time > third_positive_date ~ 4,
      V1_sero == 0 & total_inf == 0 & calendar_time < V2_date ~ 0,
      V1_sero == 0 & total_inf == 0 & calendar_time > V2_date & V2_sero2 == 0 ~ 0,
      V1_sero == 0 & total_inf == 0 & calendar_time >= V2_date & V2_sero2 == 1 ~ 1,
      V1_sero == 0 & total_inf == 0 & is.na(V2_sero2) ~ 0
    ),
    prior_infection = if_else(!is.na(First_positive_date) & First_positive_date <= Baseline_date & V1_sero == 1, prior_infection - 1, prior_infection)
  ) %>%
  mutate(
    sero = as.factor(sero),
    period = factor(period, levels = c("Pre-Delta","Delta", "Omicron")),
    prior_infection = as.numeric(prior_infection)
  )

lexis_inf_split$vax_date <- cal.yr(as.Date(lexis_inf_split$vax_date, format = "%Y-%m-%d"))

lexis_inf_split <- lexis_inf_split %>% 
  cutLexis(lexis_inf_split$vax_date, timescale = "calendar_time") %>% 
  mutate(
    vaccinated = case_when(
      calendar_time>=vax_date~1,
      calendar_time<vax_date~0, 
      is.na(vax_date) ~0, 
      TRUE ~ NA)
  )


lexis_inf_split <- lexis_inf_split %>% select(lex.id,ID, calendar_time, lex.dur, lex.Cst, lex.Xst, type, inf_num, total_inf,
                                              First_positive_date, second_positive_date, third_positive_date, V1_date, V1_spike,
                                              V1_ncp, V2_date, V2_spike, V2_ncp, V3_date, V3_spike, V3_ncp, start_date, end_date, 
                                              timein, timeout, infection, Sex, Household_ID, Total_household_no, vax_date,
                                              household_children, crowding, hh_children, sero, period, V1_sero, V2_sero2, prior_infection, age_cat, HIV, Steroid,
                                              Cancer, Diabetes, HTN, Smoking, Employed, hh_size, vaccinated, symp, latest_infection)

#Final tidying
lexis_inf_split$prior_infection[lexis_inf_split$prior_infection==3] <- 2
lexis_inf_split$prior_infection[lexis_inf_split$prior_infection==4] <- 2
lexis_inf_split$prior_infection <- factor(lexis_inf_split$prior_infection, levels=c("0", "1", "2"))
lexis_inf_split$latest_infection <- cal.yr(as.Date(lexis_inf_split$latest_infection, format = "%Y-%m-%d"))
lexis_inf_split$inf_interval <- (lexis_inf_split$calendar_time + lexis_inf_split$lex.dur) - lexis_inf_split$latest_infection
lexis_inf_split$inf_interval_cat <- cut(lexis_inf_split$inf_interval, breaks = c(0, 0.25, Inf), labels = c("0-90","90+") )
lexis_inf_split <- lexis_inf_split %>%
  mutate(combined_variable = paste(prior_infection, inf_interval_cat, sep = "+"))


#for now a quick fix to 34-399H as NO serology at baseline, one infection during follow-up
lexis_inf_split$combined_variable[lexis_inf_split$lex.id == "34-399H"] <- "0+NA"

write.csv(lexis_inf_split, "lexis_inf_split.csv")

