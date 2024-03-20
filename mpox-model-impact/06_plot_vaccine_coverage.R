# author: Fanyu Xiu
library(MpoxModelPack)
library(tidyverse)
cal_cities <- c("mtl", "trt", "van")

# reported case and vaccine data
init.pop.fn(cal_cities, 1)
load.params.fn()
data_incid <- do.call(rbind.data.frame, c(case_city, make.row.names = FALSE))
data_incid <- data_incid %>% 
  mutate(city_name = case_when(prov == "QC" ~ "Montréal",
                               prov == "ON" ~ "Toronto",
                               prov == "BC" ~ "Vancouver")) %>% 
  group_by(prov) %>% 
  mutate(first_doses_coverage = case_when(prov == "QC" ~ cumsum(first_doses) / N[["mtl"]],
                                          prov == "ON" ~ cumsum(first_doses) / N[["trt"]],
                                          prov == "BC" ~ cumsum(first_doses) / N[["van"]])) %>% 
  ungroup() %>% 
  rename(city = city_name) %>% 
  filter(date <= as.Date("2022-10-15"))

p_vc <- ggplot(data_incid, aes(x = date, y = first_doses_coverage)) +
  geom_line(aes(col = city), linewidth = 0.5) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Date", y = "Vaccine coverage") +
  scale_colour_viridis_d(option = "C", end = 0.8) +
  theme_bw() +
  theme(
    legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
    legend.background = element_rect(colour = "transparent", fill = alpha("white", 0.4)),
    legend.position = c(.15, .85),
    legend.key.size = unit(0.7, 'cm'),
    legend.direction = "vertical",
    legend.text = element_text(size = 10)
  )
p_vc
ggsave(p_vc,
       file = "./result-fig/p_vaccine_coverage.png",
       width = 15, height = 15, units = "cm", dpi = 600)
