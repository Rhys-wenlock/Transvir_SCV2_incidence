lexis_inf_split$period <- factor(lexis_inf_split$period, levels=c("Delta", "Omicron"))

test <- lexis_inf_split[lexis_inf_split$prior_infection!=0,]

test$prior_infection <- factor(test$prior_infection, levels = c("1", "2"))
newmodel <- coxph(Surv(lex.dur,lex.Xst)~period*inf_interval_cat+age_cat+vaccinated+cluster(ID,Household_ID), data = test[test$prior_infection!="0",])
summary(newmodel)




length(unique(lexis_inf_split$lex.id[lexis_inf_split$period=="Omicron" & lexis_inf_split$prior_infection!=0]))
#279 previously infected + at risk during omicron
length(unique(lexis_inf_split$lex.id[lexis_inf_split$period=="Omicron" & lexis_inf_split$prior_infection!=0 & lexis_inf_split$lex.Xst==TRUE]))
#157 people re-infected during Omicron
#56.3 (50.3 - 62.2)

length(unique(lexis_inf_split$lex.id[lexis_inf_split$period=="Delta" & lexis_inf_split$prior_infection!=0]))
#281 previously infected + at risk during delta
length(unique(lexis_inf_split$lex.id[lexis_inf_split$period=="Delta" & lexis_inf_split$prior_infection!=0 & lexis_inf_split$lex.Xst==TRUE]))
#75 people re-infected during delta
#26.7 (21.6 - 32.3)




length(unique(lexis_inf_split$lex.id[lexis_inf_split$period=="Omicron" & lexis_inf_split$prior_infection==0]))

