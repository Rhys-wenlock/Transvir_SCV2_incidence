
require(readr)
require(vcdExtra)
require(readxl)
require(ggplot2)
require(tidyverse)

#data needed
#All inf data -> my cleaned/reformatted infection episode data
#owid_covid_data -> gambia case data
#long serol visits -> long format list of serology dates
long_serol_visits <- read_csv(here::here("Data", "long_serol_visits.csv"))


owid_covid_data <- read_csv(here::here("Data","owid-covid-data.csv"))
gambia_covid_daily <- owid_covid_data[owid_covid_data$location == "Gambia", c(3,4,6)]
gambia_case_linelist<- expand.dft(gambia_covid_daily, freq="new_cases")
gambia_case_linelist$date <- as.Date(gambia_case_linelist$date)
gambia_case_linelist$location <- "Gambia"

Tranvsir_curve <- as.data.frame(all_inf_data[all_inf_data$type=="PCR_pos",3])###Consider including PCR negatives here too
Tranvsir_curve$location <- "Transvir"
colnames(Tranvsir_curve)[1] <- "date"
Tranvsir_curve$date <- as.Date(Tranvsir_curve$date)

colnames(long_serol_visits)[2] <- "date"
long_serol_visits <- long_serol_visits[,-1]
complete_epi_curve <- rbind(gambia_case_linelist, Tranvsir_curve, long_serol_visits)

ancestral_period <- interval("2020-01-01", "2021-07-07")
delta_period <- interval("2021-07-07", "2021-12-04")
complete_epi_curve$wave <- 
  ifelse(complete_epi_curve$date%within%ancestral_period, "Pre-Delta", 
         ifelse(complete_epi_curve$date%within%delta_period, "Delta", "Omicron"))


weekly_breaks <- seq.Date(from = as.Date("2021-03-02"),
                          to = as.Date("2022-05-31"),
                          by = "week")

complete_epi_curve$location <- factor(complete_epi_curve$location, levels = 
                                        c("Gambia", "Transvir", "Serology"))

complete_epi_curve$type <- ifelse(complete_epi_curve$location=="Gambia"|complete_epi_curve$location=="Transvir", complete_epi_curve$wave, "Serology")
complete_epi_curve$type <- factor(complete_epi_curve$type, levels=c("Pre-Delta", "Delta", "Omicron", "Serology"))
complete_epi_curve$type[complete_epi_curve$location=="Serology" & complete_epi_curve$date <"2021-07-01"] <- "Bleed 1"
complete_epi_curve$type[complete_epi_curve$location=="Serology" & complete_epi_curve$date >="2021-07-01" & complete_epi_curve$date<"2022-01-01"] <- "Bleed 2"

complete_epi_curve <- complete_epi_curve %>%
  mutate(type=case_when(
    location == "Serology" & date < "2021-07-01"~"Bleed 1", 
    location == "Serology" & date >= "2021-07-01" & date <"2022-01-02"~"Bleed 2", 
    location == "Serology" & date >= "2021-01-01"~"Bleed 3",
    TRUE ~ type
  )) %>% mutate(location = case_when(
    location == "Serology"~"Bleed dates", 
    TRUE ~ location
  ))

cols <- c("Pre-Delta"="#36013F","Delta"="#1974D2","Omicron"= "#EC9706","Bleed 1"= "black", "Bleed 2" = "darkgrey", "Bleed 3" = "lightgrey")
complete_epi_curve$type <- factor(complete_epi_curve$type, levels=c(
  "Pre-Delta","Bleed 1", "Delta",   "Bleed 2","Omicron", "Bleed 3"
))
complete_epi_curve$location <- factor(complete_epi_curve$location, levels = c(
  "Gambia", "Transvir", "Bleed dates"
))

figure_2 <-  ggplot(complete_epi_curve) + 
  
  geom_histogram(
    mapping = aes(
      x = date, color = type, fill= type),    # arguments inside aes() apply by group
    stat = "count",
    # arguments outside aes() apply to all data
    
    # histogram breaks
    breaks = weekly_breaks, # pre-defined date vector (see earlier in this page)
    closed = "left" # count cases from start of breakpoint
  ) + scale_fill_manual(values=cols) + scale_color_manual(values=cols)+
  theme_bw() +
  
  facet_wrap(
    location~., ncol=1, scales = "free_y") + ylab("") + xlab("Week") + 
  theme(
    # Remove panel border
    panel.border = element_blank(),
    axis.line = element_line(colour = "grey"), legend.position = "bottom", legend.title = element_blank(),
    axis.title.x = element_text(vjust=-0.5), plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"), ,
    panel.grid.major = element_blank()
  ) + scale_x_date(breaks = as.Date(c("2020-03-01", "2020-06-01", "2020-09-01", "2020-12-01",
                                      "2021-03-01", "2021-06-01", "2021-09-01", "2021-12-01",
                                      "2022-03-01", "2022-06-01", "2022-09-01")),
                   labels = c("Mar-20", "Jun-20", "Sept-20", "Dec-20",
                              "Mar-21", "Jun-21", "Sept-21", "Dec-21",
                              "Mar-22", "Jun-22", "Sept-22"))

ggsave_workaround <- function(g){survminer:::.build_ggsurvplot(x = g,
                                                               surv.plot.height = NULL,
                                                               risk.table.height = NULL,
                                                               ncensor.plot.height = NULL)}

g_to_save <- ggsave_workaround(figure_2)
ggsave("figure_2.pdf", figure_2, width = 10, height = 6, dpi=600)


