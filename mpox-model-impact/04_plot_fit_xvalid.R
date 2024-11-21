# author: Fanyu Xiu, Jorge Luis Flores Anato
library(MpoxModelPack)
library(tidyverse)
library(gridExtra)
cal_cities <- c("mtl", "trt", "van")
# reported case data
init.pop.fn(cal_cities, 1)
load.params.fn()

data_incid <- do.call(rbind.data.frame, c(case_city, make.row.names = FALSE))
data_incid <- data_incid %>% 
  mutate(city_name = case_when(prov == "QC" ~ "Panel A): Montréal",
                               prov == "ON" ~ "Panel B): Toronto",
                               prov == "BC" ~ "Panel C): Vancouver"))
data_incid <- subset(data_incid, time_conti < 150)

data_xval_obs <- read.csv("./observed_case_features.csv")
data_xval_prop_case_age <- readRDS("./out/main/data_xval_prop_case_age.rds") %>% 
  mutate(metric = "Panel A): Proportion of cumulative cases by age groups")
data_xval_prop_vaccine_age <- readRDS("./out/main/data_xval_prop_vaccine_age.rds") %>% 
  mutate(metric = "Panel B): Proportion of cumulative vaccination by age groups")
data_xval_prop_case_hiv <- readRDS("./out/main/data_xval_prop_case_hiv.rds") %>% 
  mutate(grp = ifelse(grp == "0", "Seronegative/Unknown", "Seropositive"),
         metric = "Panel C): Proportion of cumulative cases by HIV-serostatus")
data_xval_cum_vaccine <- readRDS("./out/main/data_xval_cum_vaccine.rds") %>% 
  mutate(metric = "Cumulative vaccination") 

data_xval <- rbind(data_xval_prop_case_age,
                   data_xval_prop_vaccine_age, 
                   data_xval_prop_case_hiv, 
                   data_xval_cum_vaccine) %>% 
  rename(city = city_name)



cal_analysis <- c("main",
                  "RR_in", "RR_1",
                  "contact_15", "contact_10", 
                  "VE_lb", "VE_ub", "standardized_vaccine_date",
                  "VE_1", "prioritize_vaccine")[c(3,9,10)]
for(analysis in cal_analysis){
data_mod_fit <- readRDS(sprintf("./out/%s/data_mod_fit.rds", analysis)) %>% 
  mutate(city_name = case_when(city_name == "Montréal" ~ "Panel A): Montréal",
                               city_name == "Toronto" ~ "Panel B): Toronto",
                               city_name == "Vancouver" ~ "Panel C): Vancouver"))
p_fit <- ggplot(data_mod_fit,
                      aes(x = date,
                          y = cases,
                          col = city_name)) +
    # model fit and CrI's
    geom_ribbon(aes(ymin = cases_lci,
                    ymax = cases_uci,
                    fill = city_name,
                    col = NULL),
                alpha = 0.4) +
    geom_line(aes(linetype = "Model fit"), linewidth = 0.7) +
    # observed incidence data
    geom_point(
      data = data_incid,
      aes(y = incidence, shape = "Observed"),
      size = 1.4,
      alpha = .4
    ) +
    # plot scale and facetting
    facet_wrap(~ city_name, ncol = 1, scales = "free_y", strip.position = "bottom") +
    coord_cartesian(xlim = range(data_incid$date)) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    # colour scheme
    scale_colour_viridis_d(option = "C", end = 0.8) +
    scale_fill_viridis_d(option = "C", end = 0.8) +
    # legends
    labs(x = "Date", y = "Reported mpox cases (day)", linetype = NULL, shape = NULL) +
    guides(
      shape = guide_legend(order = 1, override.aes = list(size = 2.8)),
      colour = "none",
      fill = "none"
    ) +
    theme(axis.title = element_text(size = 9),
          axis.text = element_text(size = 7),
          strip.text = element_text(size = 7),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 7),
          plot.title = element_text(size = 10),
          plot.title.position = "plot",
          plot.caption.position = "plot") +
    theme_bw() +
    theme(
      legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
      legend.background = element_rect(colour = "transparent", fill = alpha("white", 0.4)),
      legend.position = c(.82, .91),
      legend.spacing.y = unit(-2, "pt"),
      legend.direction = "vertical",
      strip.background = element_blank(),
      strip.placement = "outside"
    )

p_fit
ggsave(p_fit,
       file = sprintf("./result-fig/model_fit/%s.jpeg", analysis),
       width = 90, height = 135, units = "mm", dpi = 500)}

p_xvalid <- ggplot(filter(data_xval_obs, grp != "vaccine"),
                   aes(x = grp,
                       y = estimate,
                       col = city)) +
  geom_point(size = 0.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(filter(data_xval, grp != "vaccine"),
                mapping = aes(ymin = lci, 
                              ymax = uci), 
                width = 0.3, 
                linewidth = 0.5,
                position = position_dodge(width = 0.5)) + 
  facet_wrap(~ metric, ncol = 1, scales = "free_x", strip.position = "bottom") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(labels = scales::percent) + 
  # colour scheme
  scale_colour_viridis_d(option = "C", end = 0.8) +
  scale_fill_viridis_d(option = "C", end = 0.8) +
  # legends
  labs(x = "Group", y = "Proportion", linetype = NULL, shape = NULL) +
  guides(
    shape = guide_legend(order = 1, override.aes = list(size = 2.8))
  ) +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 7),
        strip.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        plot.title = element_text(size = 10),
        plot.title.position = "plot",
        plot.caption.position = "plot") +
  theme_bw() +
  theme(
    legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
    legend.background = element_rect(colour = "transparent", fill = alpha("white", 0.4)),
    legend.position = c(.89, .92),
    legend.spacing.y = unit(1, "pt"),
    legend.direction = "vertical",
    legend.key.size = unit(0.7, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    strip.background = element_blank(),
    strip.placement = "outside"
  ) 
p_xvalid

ggsave(p_xvalid,
       file = "./result-fig/p_fit_xvalid.jpeg",
       width = 12, height = 18, units = "cm", dpi = 600)
