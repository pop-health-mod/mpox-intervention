# author: Fanyu Xiu
# library and data --------
library(dplyr)
library(readr)
library(survey)
source("./src/helper_dates.R")
source("./src/helper_ipcw.R")

select <- dplyr::select
data_full <- read_csv("./data-3cities-feb-2023/engage_visits_3cities.csv")

# for regression analysis: prepare a data including all visits in 2022 -------
data_22 <- data_full %>% 
  filter(year(date_intv) %in% 2022)

recall_period <- as.numeric(difftime(DATE_MPOX_START %m+% months(6), DATE_MPOX_START, units = "days"))

data_22 <- data_22 %>% 
  mutate(date_intv = as.Date(date_intv),
         # period_mpox indicate a visit's coverage of peak of mpox outbreak over the 6 months recall period
         period_mpox = ifelse(
                              date_intv >= DATE_MPOX_START & date_intv < DATE_MPOX_VA2,
                              (date_intv - DATE_MPOX_START + 1) / recall_period,
                              0))

# for sensitivity analysis (6.14 endpoint)
data_22_sen6 <- data_22 %>% 
  mutate(date_intv = as.Date(date_intv),
         period_mpox = ifelse(
           date_intv >= DATE_MPOX_START & date_intv < DATE_MPOX_VA,
           (date_intv - DATE_MPOX_START + 1) / recall_period,
           0))

# for sensitivity analysis (7.14 endpoint)
data_22_sen7 <- data_22 %>% 
  mutate(date_intv = as.Date(date_intv),
         period_mpox = ifelse(
           date_intv >= DATE_MPOX_START & date_intv < DATE_MPOX_VA %m+% months(1),
           (date_intv - DATE_MPOX_START + 1) / recall_period,
           0))

# ## check retention from survey baseline
length(unique(data_22$part_id)) / length(unique(filter(data_full, visit_num == 1)$part_id)) # 0.59
sum(data_22$period_mpox != 0) # number of mpox visits

# "baseline" characteristics for participants who had a visit in 2022 
# "baseline" in this analysis is their last visit before 2022
# to prepare for "baseline" covariates in regression model
data_last_visit_before_22 <- data_full %>% 
  subset(part_id %in% unique(data_22$part_id) & date_intv < as.Date("2022-01-01")) %>%
  group_by(part_id) %>%
  mutate(last_visit_date = max(date_intv),
         last_visit_num = max(visit_num)) %>%
  subset(date_intv == last_visit_date & visit_num == last_visit_num) %>%
  select(- c("last_visit_date", 
             "last_visit_num"))
  
nrow(data_last_visit_before_22) == length(unique(data_22$part_id))

# categorize sexual activity category based on "baseline" reported # partners
data_sa_cat <- data_last_visit_before_22 %>% 
  mutate(sa_cats_ttl = ifelse(nb_part_ttl > 7, ">7", "0-7"),
         sa_cats_anal = ifelse(nb_part_anal > 7, ">7", "0-7")) %>% 
  mutate(sa_cats_ttl = factor(sa_cats_ttl, levels = c("0-7", ">7")),
         sa_cats_anal = factor(sa_cats_anal, levels = c("0-7", ">7"))) %>%
  ungroup() %>%
  dplyr::select(part_id, 
                sa_cats_ttl,
                sa_cats_anal)

table(data_sa_cat$sa_cats_ttl)

### at analysis baseline, small group constituted of 72% of the participants
sum(data_last_visit_before_22$nb_part_ttl <= 5)/nrow(data_last_visit_before_22)

data_last_visit_before_22 <- left_join(data_last_visit_before_22, data_sa_cat, by = "part_id")

# save "baseline" data to make table 1 
write.csv(data_last_visit_before_22, "./data-3cities-feb-2023/data_last_visit_before_22.csv")

# combine "baseline" characteristics with 2022 data for regression analysis
df_list <- list(data_22, data_22_sen6, data_22_sen7)
for(i in 1:3){
  df_list[[i]] <- full_join(df_list[[i]],
                            data_last_visit_before_22,
                            suffix = c("_visit", "_last"), 
                            by = c("part_id", "city", "wt_rds", "wt_rds_norm", "network_sz_rds"))
  df_list[[i]] <- df_list[[i]] %>% 
    mutate(month_intv = month(date_intv_visit),
         month_intv_centre = scale(month_intv, scale = F),
         period_mpox_centre = scale(period_mpox, scale = F), ### centre the continuous variable to avoid colinearity issue
         mpox_visit = ifelse(period_mpox == 0, 0, 1)) # if a visit is a visit that took place during the peak of the mpox outbreak
}
data_22 <- df_list[[1]]
data_22_sen6 <- df_list[[2]]
data_22_sen7 <- df_list[[3]]
# "visit" means the value of a variable at the interview date
# "last" means the value reported at the "baseline"
table(data_22$sa_cats_ttl, data_22$nb_part_ttl_visit) 
### most people who had higher number of partners at the analysis baseline 
### reported high number at the following visit in 2022 as well

## Compute IPCWs ----
## survey baseline
data_fu_pre_covid <- data_full %>% 
  group_by(part_id) %>%
  subset(visit_num == 1) %>% 
  mutate(time_pt = "Pre-COVID")
data_fu_pre_covid %>% 
  group_by(city) %>%
  summarize(sum_rds_wt = sum(wt_rds_norm), 
            n = n())
length(unique(data_fu_pre_covid$part_id))

# during 2022 
data_fu_2022 <- data_full %>% 
  # keep only visits during 2022 and before mpox started
  filter(year(date_intv) == 2022) %>%
  # use latest visit 
  group_by(part_id) %>% 
  filter(date_intv == max(date_intv) & visit_num == max(visit_num)) %>% 
  mutate(time_pt = "2022")
data_visit_key_date <- bind_rows(data_fu_pre_covid, data_fu_2022)

# retention rate
# compute proportion retained at every time point
tbl_retain <- data_visit_key_date %>% 
  mutate(time_pt = factor(time_pt, levels = c("Pre-COVID", "2022"))) %>% 
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

# lists to store each city's dataset with IPCWs
ipcw_2022 <- create_city_list()

covariates <- c("apps_partn_m",
                "apps_partn_d",
                "rel_status",
                "age", 
                "hiv_cats",
                "education_level_cat",
                "income_level_cat",
                "ethnicity_cat",
                "sex_work_m",
                "sex_work_d",
                "nb_part_ttl") # variables associated with ltfu

## check for LTFU and create indicator variable
# check if participant has a 2022 visit
data_ipcw_2022 <- data_fu_pre_covid %>% 
  mutate(
    ltfu = (part_id %in% unique(data_fu_2022$part_id))
  )

# transform variables
data_ipcw_2022 <- data_ipcw_2022 %>% 
  mutate(across(all_of(covariates), as.factor)) %>%
  mutate(across(all_of(covariates), as.integer))

## generate IPCW weights for all cities
for (cur_city in CITIES){
  ipcw_2022 <- compute_ipcw(cur_city, 
                            covariates,
                            data_timept = data_ipcw_2022,
                            data_list = ipcw_2022)
}

# turn into a single dataset
ipcw_2022 <- do.call(bind_rows, ipcw_2022)

ipcw_pre_covid <- data_visit_key_date %>%
  subset(time_pt == "Pre-COVID") %>%
  mutate(ipw_rds = wt_rds_norm)

ipcw_2022 <- merge.data.frame(data_visit_key_date, ipcw_2022, by = "part_id") %>%
  subset(time_pt == "2022")

ipcw_2022 %>% 
  group_by(city) %>%
  summarize(n = n(), sum(wt_rds_norm), sum(ipw_rds)) 

data_22 <- left_join(data_22,
                     dplyr::select(ipcw_2022, c("part_id", "ipw_rds")),
                     by = "part_id")

## save datasets
write.csv(data_22, 
          "./data-3cities-feb-2023/behavourial_regression_3cities.csv", 
          row.names = F)
write.csv(data_22_sen6, 
          "./data-3cities-feb-2023/behavourial_regression_3cities_sen_date6.csv", 
          row.names = F)
write.csv(data_22_sen7, 
          "./data-3cities-feb-2023/behavourial_regression_3cities_sen_date7.csv", 
          row.names = F)
# ESS
data_22 <- data_22 %>% 
  group_by(part_id) %>%
  mutate(at_least1_mpox_visit = ifelse(sum(mpox_visit) == 0, 0, 1))

ipcw_2022 <- ipcw_2022 %>% 
  left_join(select(data_22, c("part_id", "at_least1_mpox_visit")))

df_ess <- vector("list", 2)
names(df_ess) <- c(0, 1)

for (mpox_vis in names(df_ess)){
  print(sprintf("%s ===================", mpox_vis))
  data_tmp <- ipcw_2022 %>% 
    filter(at_least1_mpox_visit == mpox_vis)
  
  print( nrow(data_tmp) )
  print( sum(data_tmp$wt_rds) )
  print( sum(data_tmp$wt_rds_norm) )
  
  # get deff
  df_svy <- svydesign(ids = ~0, data = data_tmp, weights = ~wt_rds)
  avg <- svymean(~ nb_part_ttl, design = df_svy, deff = TRUE)
  
  df_ess[[mpox_vis]] <- data.frame(mpox_visit = mpox_vis,
                              n = nrow(data_tmp),
                              deff = data.frame(avg)$deff,
                              ess = nrow(data_tmp) / data.frame(avg)$deff)
  
  print( round(df_ess[[mpox_vis]]$ess, 2) )
}
df_ess <- bind_rows(df_ess)