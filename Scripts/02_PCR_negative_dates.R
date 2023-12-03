require(lubridate)
require(dplyr)
require(dplyr)
require(purrr)
require(lubridate)

gambia_covid_daily <- read.csv(here::here("Data", "gambia_covid_daily.csv"))

# Assuming your original dataframe is 'df_original' with 'date' and 'infections'
df_original <- gambia_covid_daily
df_original$date <- as.Date(df_original$date)
colnames(df_original)[4] <- "infections"

# Assuming your original dataframe is 'df_original' with 'date' and 'infections'
df_original$date <- as.Date(df_original$date)

# Create a week_start variable
df_original$week_start <- floor_date(df_original$date, "week")

# Sum infections by week_start
infections_by_week_start <- aggregate(infections ~ week_start, data = df_original, FUN = sum)

# Create a week variable that represents the middle of the week
infections_by_week_start$week <- infections_by_week_start$week_start + days(3)
# Now infections_by_week_start has the middle date of each week (in the 'week' column) 
# and the sum of infections for that week (in the 'infections' column)

infections_by_week_start


# Function to sample weeks for a single participant
sample_weeks <- function(start_date, end_date, size = 1) {
  infections_by_week_start %>%
    filter(week >= start_date & week <= end_date) %>%
    mutate(prob = infections / sum(infections)) %>%
    {sample(.$week, size, replace = TRUE, prob = .$prob)}
}
