###Table 1 


##need to load full_demog_data

# If not already installed, install the necessary packages
# install.packages("dplyr")
# install.packages("tidyverse")

# Load the packages
library(dplyr)
library(tidyverse)

# Create a function to compute the count and percentage
count_pct <- function(var) {
  df <- data.frame(value = var)
  df %>%
    group_by(value) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup() %>%
    transmute(Level = value, Count_Percentage = paste0(n, " (", round(freq * 100, 1), "%)"))
}

# Compute the summary for each categorical variable
age_cat_summary <- count_pct(full_demog_data$age_cat)
sex_summary <- count_pct(full_demog_data$Sex)
V1_spike_summary <- count_pct(full_demog_data$V1_spike)
V2_spike_summary <- count_pct(full_demog_data$V2_spike)
V3_spike_summary <- count_pct(full_demog_data$V3_spike)
V1_ncp_summary <- count_pct(full_demog_data$V1_ncp)
V2_ncp_summary <- count_pct(full_demog_data$V2_ncp)
V3_ncp_summary <- count_pct(full_demog_data$V3_ncp)
vac_6m_summary <- count_pct(full_demog_data$vaccinated_6m)
vac_12m_summary <- count_pct(full_demog_data$vaccinated_12m)

# Compute the median and IQR for the continuous variable
number_of_participants_summary <- tibble(
  Level = c("Median"),
  Count_Percentage = paste0(
    formatC(median(full_demog_data$Total_household_no, na.rm = TRUE), format = "f", digits = 1),
    " (",
    formatC(quantile(full_demog_data$Total_household_no, 0.25, na.rm = TRUE), format = "f", digits = 1),
    " - ",
    formatC(quantile(full_demog_data$Total_household_no, 0.75, na.rm = TRUE), format = "f", digits = 1),
    ")"
  )
)

number_of_rooms_summary <- tibble(
  Level = c("Median"),
  Count_Percentage = paste0(
    formatC(median(full_demog_data$Total_rooms, na.rm = TRUE), format = "f", digits = 1),
    " (",
    formatC(quantile(full_demog_data$Total_rooms, 0.25, na.rm = TRUE), format = "f", digits = 1),
    " - ",
    formatC(quantile(full_demog_data$Total_rooms, 0.75, na.rm = TRUE), format = "f", digits = 1),
    ")"
  )
)

# Combine all summaries into a single data frame
summary_table <- bind_rows(
  mutate(age_cat_summary, Variable = "Age Category"),
  mutate(sex_summary, Variable = "Sex"),
  mutate(V1_spike_summary, Variable="V1_spike"), 
  mutate(V2_spike_summary, Variable="V2_spike"), 
  mutate(V3_spike_summary, Variable="V3_spike"), 
  mutate(V1_ncp_summary, Variable="V1_ncp"),
  mutate(V2_ncp_summary, Variable="V2_ncp"), 
  mutate(V3_ncp_summary, Variable="V3_ncp"), 
  mutate(vac_6m_summary, Variable = "Vaccination at 6months"),
  mutate(vac_12m_summary, Variable = "Vaccination at 12months"),
  mutate(number_of_participants_summary, Variable = "Number in household"), 
  mutate(number_of_rooms_summary, Variable = "Number of rooms")
)

# Print the summary table
View(summary_table)
write.csv(summary_table, "Table1.csv")
