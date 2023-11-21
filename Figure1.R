# Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)


visit_data$PCR_result[visit_data$PCR_result=="N/A"] <- NA
df <- visit_data
full_demog_data
# Function to assign results based on conditions
assign_results <- function(data, column) {
  data %>%
    mutate(Result = case_when(
      !!rlang::sym(column) == "Pos" ~ "Seropositive",
      !is.na(!!rlang::sym(column)) ~ "Seronegative",
      TRUE ~ NA_character_
    )) %>%
    rename(Date = 2)  # Rename the second column to 'Date'
}
full_demog_data <- full_demog_data %>%
  left_join(All_visit_serol, by="Participant_ID")

# Apply function to each data frame
V1_seroldate <- assign_results(full_demog_data[,c(1,25,26)], "V1_spike")
V2_seroldate <- assign_results(full_demog_data[,c(1,28,29)], "V2_spike")
V3_seroldate <- assign_results(full_demog_data[,c(1,31,32)], "V3_spike")

# Combine data frames and replace NA in Result with 'Missing'
df <- bind_rows(visit_data %>% mutate(Result = ifelse(PCR_result=="Positive", "RT-PCR positive", ifelse(!is.na(PCR_result), "RT-PCR negative", NA))),
                V1_seroldate, V2_seroldate, V3_seroldate)

df$Date <- as.Date(df$Date)

# Add year_week column
df <- df %>%
  mutate(Date = as.Date(Date),
         year = year(Date),
         week = week(Date),
         year_week = paste(year, sprintf("%02d", week), sep="-"))

# Define priority function
prioritize_results <- function(results) {
  priority_order = c("Seropositive", "Seronegative", "RT-PCR positive", "RT-PCR negative", "Missing")
  result = intersect(priority_order, results)
  if(length(result) > 0) result[1] else NA_character_
}

# Group and summarise
df_summary <- df %>%
  group_by(Participant_ID, year, week) %>%
  summarize(result = prioritize_results(Result), .groups = "drop")  %>%
  mutate(year_week = paste(year, sprintf("%02d", week), sep="-"))  # Recombine into year_week

# Filter the data for the desired week range
df_summary_filtered <- df_summary %>%
  filter((year == 2021 & week >= 10) | (year == 2022 & week <= 52))

df_summary_filtered$result[is.na(df_summary_filtered$result)] <- "Missing"

# Manually specify the breaks
breaks_year_week <- c("2021-10", "2021-23", "2021-36", "2021-49", "2022-10", "2022-23") # Assuming that these weeks correspond to the middle of the months you mentioned

# Corresponding labels
breaks_labels <- c("Mar-21", "Jun-21", "Sep-21", "Dec-21", "Mar-22", "Jun-22")

#Reorder results
df_summary_filtered$result <- factor(df_summary_filtered$result, levels = c("RT-PCR positive", "RT-PCR negative", "Seropositive", "Seronegative", "Missing"), 
                                     labels=c("RT-PCR positive", "RT-PCR negative", "Seropositive", "Seronegative", "Missing"))

# Plot your data

Figure_1 <- ggplot(df_summary_filtered, aes(x = year_week, y = Participant_ID, fill = result)) +
  geom_tile(colour = "white", width = 0.9, height = 1.8) +
  scale_fill_manual(values = c("RT-PCR positive" = "red", "RT-PCR negative" = "skyblue1", "Seropositive" = "purple", "Seronegative" = "orange", "Missing" = "white")) +
  scale_x_discrete(breaks = breaks_year_week, labels = breaks_labels) +
  scale_y_discrete(breaks = NULL) + scale_x_discrete(
    breaks = breaks_year_week, 
    labels = breaks_labels,
    sec.axis = sec_axis(~., breaks = breaks_year_week, labels = breaks_labels)  # Secondary axis
  )+
  labs(x = "Date", y = "Study participants", fill = "") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size=15),
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 15),  # Adjust x-axis title text size
    axis.title.y = element_text(size = 15),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 15),  # Adjust legend text size
    legend.key.size = unit(1.5, "lines"),  # Increase legend key size
    legend.margin = margin(10, 0, 10, 0)  # Adjust legend margin
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(shape = 15, size = 1, colour = "black"),
      keywidth = unit(1, "lines"),  # Increase legend key width
      keyheight = unit(1, "lines")  # Adjust legend key height
    )
  )

ggsave("figure_1.png", Figure_1, width = 10, height = 10, dpi=600)
