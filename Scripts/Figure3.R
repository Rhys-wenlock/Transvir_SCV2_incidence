##Figure 3 Plotting

require(survminer)
require(survival)

#data requirements
#lexis_inf_split

ggsave_workaround <- function(g){survminer:::.build_ggsurvplot(x = g,
                                                               surv.plot.height = NULL,
                                                               risk.table.height = NULL,
                                                               ncensor.plot.height = NULL)}

lexis_inf_split$prior_infection[lexis_inf_split$prior_infection==3] <- 2
lexis_inf_split$prior_infection[lexis_inf_split$prior_infection==4] <- 2
lexis_inf_split$prior_infection <- factor(lexis_inf_split$prior_infection, levels=c("0", "1", "2"))
lexis_inf_split$period <- factor(lexis_inf_split$period, levels=c("Delta", "Omicron"))


plot_data <- lexis_inf_split

plot_data1 <- plot_data %>% group_by(ID, prior_infection) %>% summarize(lex.dur2=sum(lex.dur), keep.all=TRUE)
plot_data$lex.Xst <- factor(plot_data$lex.Xst, levels=c(TRUE, FALSE))
plot_data2 <- plot_data %>% group_by(ID, prior_infection) %>% arrange(lex.Xst) %>% filter(row_number()==1)
plot_data_use <- left_join(plot_data1, plot_data2, by=c("ID","prior_infection"))
plot_data_use$infection <- ifelse(plot_data_use$lex.Xst==TRUE,1,0)
plot_data_use$prior_infection <- as.factor(plot_data_use$prior_infection)
single_plot <- survfit(Surv(lex.dur2, infection)~prior_infection,data = plot_data_use)


plot_data_use$prior_infection <- factor(plot_data_use$prior_infection, levels = c("0", "1", "2"), labels=c("No prior infections", 
                                                                                                           "1 prior infection", "2+ prior infections"))
# Create the initial ggsurvplot
figure_3a <- ggsurvplot(
  single_plot, 
  palette = c("#36013F", "#1974D2", "#EC9706"),
  pval = TRUE, 
  risk.table = TRUE, 
  risk.table.col = "strata",
  risk.table.height = 0.2, 
  fun = "event", 
  conf.int = TRUE,
  tables.theme = clean_theme(), 
  xlim = c(0, 1), 
  legend.title = "", 
  censor = FALSE, 
  legend.labs = levels(plot_data_use$prior_infection),
  xlab = "Follow-up time (years)",
  ylab = "Cumulative incidence (risk)")

# Modify the theme of the entire plot to add gridlines
figure_3a$plot <- figure_3a$plot + 
  theme(legend.text = element_text(size = 15),  # Increase the font size of the legend text
        axis.title.y = element_text(margin = margin(r = 10)), 
         panel.grid.major = element_line(color = "grey95", size = 0.5)) + scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1.0))+
  xlab("Follow-up time (years)") + 
  ylab("Cumulative incidence (risk)")

# Adjust the risk table's theme to remove gridlines
figure_3a$table <- figure_3a$table + 
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(), 
    text =  element_text(size=15)
  )

# Print the final plot
print(figure_3a)

g_to_save <- ggsave_workaround(figure_3a)
ggsave("figure_3a.png", g_to_save, width = 10, height = 6, dpi=600)


#Figure 3b


figure_3b_plot <- survfit(Surv(lex.dur, lex.Xst)~prior_infection,data = lexis_inf_split)

figure_3b <- ggsurvplot_facet(figure_3b_plot, pval = FALSE, data=lexis_inf_split,
                                           risk.table = TRUE,
                                           facet.by = "period",
                                           risk.table.col = "strata",
                                           risk.table.height = 0.2, , fun="event", conf.int = TRUE,
                                           tables.theme = clean_theme(), xlim=c(0,0.4), 
                                           legend.title="", censor=FALSE,  palette = c("#36013F", "#1974D2",  "#EC9706")) + xlab("Follow-up time (years)") + 
  ylab("Cumulative incidence (risk)") +  theme(
    panel.background = element_rect(fill = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_line(color = "grey95", size = 0.5),
    legend.text = element_text(size=10), strip.text = element_text(size=15)
  ) 



print(figure_3b)

ggsave("figure_3b.png", figure_3b, width = 10, height = 6, dpi=600)

