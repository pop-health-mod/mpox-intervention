# author: Fanyu Xiu
library(MpoxModelPack)
library(tidyverse)
cal_cities <- c("mtl", "trt", "van")
# plot_analysis <- c("main", "contact_15", "contact_10", "VE_lb", "VE_ub", "RR_in")
plot_analysis <- "main"

# reported case data
init.pop.fn(cal_cities, 1)
load.params.fn()
data_incid <- do.call(rbind.data.frame, c(case_city, make.row.names = FALSE))
data_incid <- data_incid %>% 
  mutate(city_name = case_when(prov == "QC" ~ "Panel A): Montréal",
                               prov == "ON" ~ "Panel B): Toronto",
                               prov == "BC" ~ "Panel C): Vancouver"))
data_incid <- subset(data_incid, time_conti < 150)

for(analysis in plot_analysis){
  data_mod_fit <- readRDS(sprintf("./out/%s/data_mod_fit.rds", analysis)) %>% 
    mutate(vaccine_date = ifelse(city_name == "Montréal", as.Date("2022-06-03"), ifelse(city_name == "Toronto", as.Date("2022-06-12"), as.Date("2022-06-20")))) %>% 
    mutate(city_name = case_when(city_name == "Montréal" ~ "Panel A): Montréal",
                          city_name == "Toronto" ~ "Panel B): Toronto",
                          city_name == "Vancouver" ~ "Panel C): Vancouver")) 
  
  data_nothing_fit <- readRDS(sprintf("./out/%s/data_nothing_fit.rds", analysis)) %>% 
    mutate(city_name = case_when(city_name == "Montréal" ~ "Panel A): Montréal",
                                 city_name == "Toronto" ~ "Panel B): Toronto",
                                 city_name == "Vancouver" ~ "Panel C): Vancouver"))
  data_behavourial_fit <- readRDS(sprintf("./out/%s/data_behavourial_fit.rds", analysis)) %>% 
    mutate(city_name = case_when(city_name == "Montréal" ~ "Panel A): Montréal",
                                 city_name == "Toronto" ~ "Panel B): Toronto",
                                 city_name == "Vancouver" ~ "Panel C): Vancouver"))
  data_vaccine_fit <- readRDS(sprintf("./out/%s/data_vaccine_fit.rds", analysis)) %>% 
    mutate(city_name = case_when(city_name == "Montréal" ~ "Panel A): Montréal",
                                 city_name == "Toronto" ~ "Panel B): Toronto",
                                 city_name == "Vancouver" ~ "Panel C): Vancouver"))
  data_tracing_fit <- readRDS(sprintf("./out/%s/data_tracing_fit.rds", analysis)) %>% 
    mutate(city_name = case_when(city_name == "Montréal" ~ "Panel A): Montréal",
                                 city_name == "Toronto" ~ "Panel B): Toronto",
                                 city_name == "Vancouver" ~ "Panel C): Vancouver"))

  p_AF <- ggplot(data_mod_fit, aes(x = date, y = cases)) +
  geom_line(aes(linetype = "Observed (all three interventions)", col = "Observed (all three interventions)"), linewidth = 0.5) +
  geom_ribbon(aes(ymin = cases_lci, ymax = cases_uci, fill = "Observed (all three interventions)"), alpha = 0.4) +
  geom_line(data_nothing_fit, 
            mapping = aes(x = date, y = cases, 
                          linetype = "Unmitigated epidemic \n (without any of the three interventions)", 
                          col = "Unmitigated epidemic \n (without any of the three interventions)",
                          fill = "Unmitigated epidemic \n (without any of the three interventions)"), linewidth = 0.5 ) +
  # geom_ribbon(data_nothing_fit, mapping = aes(ymin = cases_lci, ymax = cases_uci, col = NULL, fill = "Unmitigated epidemic \n (without any of the three interventions)"), alpha = 0.3) +

  geom_line(data_behavourial_fit, mapping = aes(x = date, y = cases, linetype = "With only partner number change", col = "With only partner number change", fill = "With only partner number change"), linewidth = 0.5 ) +
  # geom_ribbon(data_behavourial_fit, mapping = aes(ymin = cases_lci, ymax = cases_uci, col = NULL, fill = "With only partner number change"), alpha = 0.3) +

  geom_line(data_vaccine_fit, mapping = aes(x = date, y = cases, linetype = "With only first-dose vaccination", col = "With only first-dose vaccination", fill = "With only first-dose vaccination"), linewidth = 0.5) +
  # geom_ribbon(data_vaccine_fit, mapping = aes(ymin = cases_lci, ymax = cases_uci, col = NULL, fill = "With only first-dose vaccination"), alpha = 0.3) +

  geom_line(data_tracing_fit, mapping = aes(x = date, y = cases, linetype = "With only contact tracing", col = "With only contact tracing", fill =  "With only contact tracing"), linewidth = 0.5) +
  # geom_ribbon(data_tracing_fit, mapping = aes(ymin = cases_lci, ymax = cases_uci, col = NULL, fill = "With only contact tracing"), alpha = 0.3) +

  # plot scale and facetting
  facet_wrap(~ city_name, ncol = 1, scales = "free_y",
             strip.position = "bottom") +
  # coord_cartesian(xlim = range(data_incid$date)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("Observed (all three interventions)" = "green4", "With only partner number change" = "#6899c9", "With only contact tracing" = "purple2", "With only first-dose vaccination" = "firebrick4", "Unmitigated epidemic \n (without any of the three interventions)" = "ivory3"),
                     breaks = c("Observed (all three interventions)", "With only partner number change", "With only contact tracing", "With only first-dose vaccination", "Unmitigated epidemic \n (without any of the three interventions)")) +
  scale_fill_manual(values = c("green4", "white", "white", "white", "white"),
                    breaks = c("Observed (all three interventions)", "With only partner number change", "With only contact tracing", "With only first-dose vaccination", "Unmitigated epidemic \n (without any of the three interventions)")) +
  scale_linetype_manual(values = c(1, 1, 1, 1, 1),
                        breaks = c("Observed (all three interventions)", "With only partner number change", "With only contact tracing", "With only first-dose vaccination", "Unmitigated epidemic \n (without any of the three interventions)")) +
  geom_vline(aes(xintercept = vaccine_date), col = "firebrick4", linetype = 2, linewidth = 0.3) +
  labs(x = "Date", y = "Reported mpox cases (day)", linetype = NULL, shape = NULL, color  = NULL, fill = NULL, alpha = NULL) +
  theme_bw() +
  theme(
    legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
    legend.background = element_rect(colour = "transparent", fill = alpha("white", 0.4)),
    legend.position = c(.77, .92),
    legend.key.size = unit(0.4, 'cm'),
    legend.direction = "vertical",
    legend.text = element_text(size = 5.3),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(size = 10)
  )
  # p_AF
  ggsave(p_AF,
  file = sprintf("./result-fig/AF_ci/%s.png", analysis),
  width = 10, height = 15, units = "cm", dpi = 600)

}
