# author: Fanyu Xiu
library(tidyverse)
cal_cities <- c("mtl", "trt", "van")
cal_analysis <- c("main", "RR_in", "contact_15", "contact_10", "VE_lb", "VE_ub", "standardized_vaccine_date")
data_results_list <- list()

for(analysis in cal_analysis){
  data_results_list[[analysis]] <- readRDS(sprintf("./out/%s/data_results.rds", analysis)) %>% 
    mutate(analysis = analysis)
}

all_results <- do.call(rbind.data.frame, c(data_results_list, make.row.names = FALSE)) %>% 
  # subset(name != "all three combined") %>%
  mutate(analysis = case_when(analysis == "contact_10" ~ "contract traced 10%",
                              analysis == "contact_15" ~ "contract traced 15%",
                              analysis == "RR_in" ~ "informative prior for RR",
                              analysis == "standardized_vaccine_date" ~ "same vaccination start relative to cases",
                              analysis == "VE_lb" ~ "use VE lower 95% CI",
                              analysis == "VE_ub" ~ "use VE upper 95% CI",
                              analysis == "main" ~ "main"),
         name = factor(name, 
                       levels = c("imported cases (tau)",
                         "transmission parameter (beta)",
                         "assortativity (omega)",
                         "rate ratio in high risk group (RR_H)",
                         "rate ratio in low risk group (RR_L)",
                         "duration infectiousness (1/gamma1)",
                         "behavourial change", 
                         "contact tracing", 
                         "first-dose vaccination", 
                         "behavourial change and contact tracing",
                         "behavourial change and first-dose vaccination",
                         "contact tracing and first-dose vaccination", 
                         "all three combined")), 
         analysis = factor(analysis, levels = c("main", "informative prior for RR", 
                                                "contract traced 10%", "contract traced 15%",
                                                "use VE lower 95% CI", "use VE upper 95% CI",
                                                "same vaccination start relative to cases"))) %>% 
  mutate(AF = as.logical(name %in% c("behavourial change", "contact tracing", "first-dose vaccination", 
                                     "behavourial change and contact tracing",
                                     "behavourial change and first-dose vaccination",
                                     "contact tracing and first-dose vaccination", 
                                     "all three combined")),
         estimate = ifelse(AF, paste0(round(estimate * 100, 0), "%"), round(estimate, 2)),
         lci = ifelse(AF, paste0(round(lci * 100, 0), "%"), round(lci, 2)),
         uci = ifelse(AF, paste0(round(uci * 100, 0), "%"), round(uci, 2)),
         `estimate (95% CrI)` = ifelse(AF, paste0(estimate, " (", lci, "-", uci, ")"), paste0(estimate, " (", lci, ",", uci, ")")),
         city_name = ifelse((!AF & !name %in% c("imported cases (tau)", "assortativity (omega)")), "all", city_name)) %>% 
  dplyr::select(analysis, city_name, name, `estimate (95% CrI)`, AF)  %>%
  distinct(analysis, city_name, name, .keep_all = TRUE) 
 
AF_results <- all_results %>% filter(AF) %>% arrange(name, analysis, city_name) %>%  dplyr::select(-AF) %>% 
  pivot_wider(names_from = "name", values_from = "estimate (95% CrI)")
# view(AF_results)
para_results <- all_results %>% filter(!AF) %>% arrange(name, analysis, city_name) %>%  dplyr::select(-AF)
write.csv( AF_results,
          "./result-tbl/AF_results_tbl.csv")
view(para_results)
write.csv( para_results,
           "./result-tbl/para_results_tbl.csv")
