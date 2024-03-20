# author: Fanyu Xiu, Jorge Luis Flores Anato
# library and data --------
library(dplyr)
library(tidyr)
library(readr)
source("./src/helper_dates.R")
source("./src/helper_ipcw.R")
data_full <- read_csv("./data-3cities-feb-2023/engage_visits_3cities.csv")
select <- dplyr::select

## Create datasets for every time period ----
## pre-pandemic time period
data_fu_pre_covid <- data_full %>% 
  group_by(part_id) %>%
  subset(visit_num == 1) %>% 
  mutate(time_pt = "Pre-COVID")

data_fu_pre_covid %>% 
  group_by(city) %>%
  summarize(sum_rds_wt = sum(wt_rds_norm), 
            n = n())
length(unique(data_fu_pre_covid$part_id))

# for contact rate: data for participants who have visits in 2022 before mpox
data_fu_pre_mpox <- data_full %>% 
  # keep only visits during 2022 and before mpox started
  filter(year(date_intv) == 2022 & date_intv < DATE_MPOX_START) %>%
  # keep latest visit 
  group_by(part_id) %>% 
  filter(date_intv == max(date_intv) & visit_num == max(visit_num)) %>% 
  mutate(time_pt = "2022-Pre-Mpox")

data_visit_key_date <- bind_rows(data_fu_pre_covid, data_fu_pre_mpox)

# retention rate
# compute proportion retained at every time point
tbl_retain <- data_visit_key_date %>% 
  mutate(time_pt = factor(time_pt, levels = c("Pre-COVID", "2022-Pre-Mpox"))) %>% 
  group_by(city, time_pt) %>% 
  summarize(nb = n(), .groups = "drop_last") %>% 
  mutate(prop = round(nb / max(nb) * 100)) %>% 
  ungroup()

# format proportion (apostrophe used to prevent excel from reading as negatives)
tbl_retain <- tbl_retain %>% 
  mutate(prop = sprintf("'(%s%%)", prop))

# format to one column per city
tbl_retain <- tbl_retain %>% 
  pivot_wider(names_from = "city", values_from = c("nb", "prop")) %>% 
  select(time_pt, ends_with("_mtl"), ends_with("_trt"), ends_with("_van"))

## Compute IPCWs ----
# lists to store each city's dataset with IPCWs
ipcw_pre_mpox <- create_city_list()

covariates <- c("apps_partn_m",
                "apps_partn_d",
                "rel_status",
                "age_cats",
                "hiv_cats",
                "education_level_cat",
                "income_level_cat",
                "ethnicity_cat",
                "sex_work_m",
                "sex_work_d",
                "nb_part_ttl") ### variables associated with ltfu

## check for LTFU and create indicator variable
# check if participant has a defined visit
data_ipcw_pre_mpox <- data_fu_pre_covid %>% 
  mutate(
    ltfu = (part_id %in% unique(data_fu_pre_mpox$part_id))
  )

# transform variables

# apps_partn, hiv_stat, sex_work are all already binary
# education_level as well but need to turn from text
data_ipcw_pre_mpox[, covariates]
for(cur_var in covariates){ print(table(data_ipcw_pre_mpox[, cur_var], useNA = "ifany")) }

data_ipcw_pre_mpox$edu <- as.integer(data_ipcw_pre_mpox$education_level_cat == "bachelor or higher")

## rearrange all other variables
# relationship status (reference is single)
data_ipcw_pre_mpox$rel_stat_exclusive <- as.integer(data_ipcw_pre_mpox$rel_status == "exclusive")
data_ipcw_pre_mpox$rel_stat_open <- as.integer(data_ipcw_pre_mpox$rel_status == "open")
data_ipcw_pre_mpox$rel_stat_unclear <- as.integer(data_ipcw_pre_mpox$rel_status == "unclear")

data_ipcw_pre_mpox %>% group_by(across(starts_with("rel_stat"))) %>% count()

# age group (reference is 18-29)
data_ipcw_pre_mpox$age_30_39 <- as.integer(data_ipcw_pre_mpox$age_cats == "30-39")
data_ipcw_pre_mpox$age_40_49 <- as.integer(data_ipcw_pre_mpox$age_cats == "40-49")
data_ipcw_pre_mpox$age_50_59 <- as.integer(data_ipcw_pre_mpox$age_cats == "50-59")
data_ipcw_pre_mpox$age_60plus <- as.integer(data_ipcw_pre_mpox$age_cats == "60+")

data_ipcw_pre_mpox %>% select(-age) %>% group_by(across(starts_with("age"))) %>% count()

# income (reference is <20k)
data_ipcw_pre_mpox$income_20_40k <- as.integer(data_ipcw_pre_mpox$income_level_cat == "20,000-40,000")
data_ipcw_pre_mpox$income_40kplus <- as.integer(data_ipcw_pre_mpox$income_level_cat == "40,000+")

data_ipcw_pre_mpox %>% 
  select(-income_level) %>% 
  group_by(across(starts_with("income"))) %>% 
  count()

# ethnicity (reference is White Canadian or European)
data_ipcw_pre_mpox$ethnic_vismin <- as.integer(!data_ipcw_pre_mpox$ethnicity_cat %in% c("White Canadian", "European"))

table(data_ipcw_pre_mpox$ethnicity_cat, data_ipcw_pre_mpox$ethnic_vismin)

# nb_partners (do <=5 vs >5)
data_ipcw_pre_mpox$nb_part_over5 <- as.integer(data_ipcw_pre_mpox$nb_part_ttl > 5)

table(data_ipcw_pre_mpox$nb_part_over5)
table(`>5` = data_ipcw_pre_mpox$nb_part_over5, true_nb = data_ipcw_pre_mpox$nb_part_ttl)

## vector of formatted variables
covariates_formatted <- c("apps_partn_m",
                          "rel_stat_exclusive", "rel_stat_open", "rel_stat_unclear",
                          "age_30_39", "age_40_49", "age_50_59", "age_60plus", 
                          "hiv_cats",
                          "edu",
                          "income_20_40k",
                          "income_40kplus",
                          "ethnic_vismin",
                          "sex_work_m",
                          "sex_work_d",
                          "nb_part_over5")

## generate IPCW weights for all cities
for (cur_city in CITIES){
  ipcw_pre_mpox <- compute_ipcw(cur_city, covariates_formatted,
                            data_timept = data_ipcw_pre_mpox,
                            data_list = ipcw_pre_mpox)
}

# turn into a single dataset
ipcw_pre_mpox <- do.call(bind_rows, ipcw_pre_mpox)

ipcw_pre_covid <- data_visit_key_date %>%
  subset(time_pt == "Pre-COVID") %>%
  mutate(ipw_rds = wt_rds_norm)

ipcw_pre_mpox <- merge.data.frame(data_visit_key_date, ipcw_pre_mpox, by = "part_id") %>%
  subset(time_pt == "2022-Pre-Mpox")

ipcw_pre_mpox %>%
  group_by(city) %>%
  summarize(n = n(), 
            sum(wt_rds_norm), 
            sum(ipw_rds))

## save datasets
write.csv(ipcw_pre_covid,
          "./data-3cities-feb-2023/survey_baseline_3cities.csv",
          row.names = F)
write.csv(ipcw_pre_mpox,
          "./data-3cities-feb-2023/data_contact_rate_3cities.csv",
          row.names = F)
