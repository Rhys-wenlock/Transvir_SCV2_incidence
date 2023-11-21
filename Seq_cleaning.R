##Sequence Results Code

Transvir_sequencing_summary_clean_02072023 <- read.csv("~/Documents/LSHTM/MSc Project/Manuscript 1/Transvir1/Transvir_sequencing_summary_clean_02072023.csv")

Transvir_sequencing_summary_clean_02072023 <- Transvir_sequencing_summary_clean_02072023 %>%
  group_by(Participant.ID, Date) %>% filter(row_number()==1)

colnames(Transvir_sequencing_summary_clean_02072023)[1] <- "Participant_ID"
Transvir_sequencing_summary_clean_02072023$Date <- as.Date(Transvir_sequencing_summary_clean_02072023$Date, format="%d/%m/%Y")

Transvir_seq <- left_join(all_inf_data[all_inf_data$type=="PCR_pos",], 
                          Transvir_sequencing_summary_clean_02072023, by=c("Participant_ID", "Date"))

#97 sequenced results 
#64 delta, 26 omicron, 7 delta

Transvir_seq <- Transvir_seq %>%  mutate(period = case_when(
  Date < "2021-07-07" ~ "Pre-Delta",
  Date < "2021-12-04" ~ "Delta",
  TRUE ~ "Omicron"
))

table(Transvir_seq$period, Transvir_seq$Variant_type)

write.csv(Transvir_seq, "Transvir_seq.csv")
