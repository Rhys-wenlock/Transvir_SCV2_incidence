length(unique(all_inf_data$Participant_ID[all_inf_data$V1_spike=="Neg"])
       + )
#149
length(unique(all_inf_data$Participant_ID[all_inf_data$V1_spike=="Neg" & all_inf_data$second_positive_date >"2021-07-07" & all_inf_data$second_positive_date <"2021-12-04"]))
#51 delta re-infections
 length(unique(all_inf_data$Participant_ID[all_inf_data$V1_spike=="Neg" & all_inf_data$second_positive_date >="2021-12-04"]))


 
 
 #Thinking about previous infection date
 
 #if 
 #first infection (inf_num=1) & baseline = pos -> prev infection = start date
 #first infection & baseline = neg -> no previous infection
 #second infection = 
 #0 infection (means end of follow-up) -> last infection
 
 date_cols <- c("First_positive_date", "second_positive_date", "third_positive_date")
all_inf_data<- all_inf_data %>%
   mutate(across(all_of(date_cols), as.Date))

gambia_covid_daily <- read.csv(here::here("Data", "gambia_covid_daily.csv"))

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


#create time interval between infection episodes
#amend this in code 04


#time intervals
<180 days = <0.5
180+ days>= 0.2465
lexis_inf_split$latest_infection <- cal.yr(as.Date(lexis_inf_split$latest_infection, format = "%Y-%m-%d"))
lexis_inf_split$inf_interval <- (lexis_inf_split$calendar_time + lexis_inf_split$lex.dur) - lexis_inf_split$latest_infection
lexis_inf_split$inf_interval_cat <- cut(lexis_inf_split$inf_interval, breaks = c(0, 0.25,0.5, Inf), labels = c("0-90","90-180", "180+") )

lexis_inf_split$prior_infection[lexis_inf_split$prior_infection==3] <- 2
lexis_inf_split$prior_infection[lexis_inf_split$prior_infection==4] <- 2

lexis_inf_split <- lexis_inf_split %>%
mutate(combined_variable = paste(prior_infection, inf_interval_cat, sep = "+"))


#for now a quick fix to 34-399H
lexis_inf_split$combined_variable[lexis_inf_split$lex.id == "34-399H"] <- "0+NA"




 