# authors: Fanyu Xiu, Jorge Luis Flores Anato
# library and data --------
library(dplyr)
source("./src/helper_dates.R")
data_22 <- read_csv("./data-3cities-feb-2023/behavourial_regression_3cities.csv")
select <- dplyr::select
table(data_22$city, useNA = "ifany")

# Data cleaning and formating ------------------
## create dummy variables

# age; reference is 16-29
table(data_22$age_cats_visit, useNA = "ifany")

# partnership status
table(data_22$rel_status_last, useNA = "ifany")
# 
# # check SPVs and groupsex
# table(data_22$bath_m_last, useNA = "ifany")
# table(data_22$groupsex_m_last, useNA = "ifany")
# 
# data_22$bath_m_last[data_22$bath_d == 1] <- "Missing"
# data_22$groupsex_m_last[data_22$groupsex_d == 1] <- "Missing"
# 
# # format sex work variable 
# table(data_22$sex_work_m, useNA = "ifany")
# data_22$sex_work_m_last[data_22$sex_work_d == 1] <- "Missing"

# hiv status 
table(data_22$hiv_cats_visit, useNA = "ifany")

# health infor too many missing values!
# table(data_22$hlth_infor_m, useNA = "ifany")
# table(data_22$hlth_infor_d, useNA = "ifany")
# data_22$hlth_infor_m[data_22$hlth_infor_d == 1] <- "Missing"

# turn into factors to ensure correct ordering and coding
data_22 <- data_22 %>%
  mutate(
    period_mpox = period_mpox * 100, # show precentage on table 
    nb_part_ttl = nb_part_ttl_visit,
    nb_part_anal = nb_part_anal_visit,
    mpox_visit = factor(mpox_visit, levels = c("0", "1")),
    sa_cats_ttl = factor(sa_cats_ttl, levels = c("0-7", ">7")),
    sa_cats_anal = factor(sa_cats_anal, levels = c("0-7", ">7")),
    age_cats = age_cats_visit,
    rel_status = factor(rel_status_last,
                        levels = c("single", "open", "exclusive", "unclear"),
                        labels = c("Single", "Open", "Exclusive", "Unclear")),
    # education = factor(education_level_cat, levels = c("lower than bachelor",
    #                                                    "bachelor or higher")),
    # sex_work = factor(sex_work_m, levels = c(1, 0, "Missing"),
    #                   labels = c("Yes", "No", "Missing")),
    hiv_cats = factor(hiv_cats_visit,
                      levels = c(1, 0),
                      labels = c("Seropositive", "Seronegative/Unknown"))
    # groupsex = factor(groupsex_m,
    #                   levels = c(1, 0, "Missing"),
    #                   labels = c("Yes", "No", "Missing")),
    # bath = factor(bath_m,
    #               levels = c(1, 0, "Missing"),
    #               labels = c("Yes", "No", "Missing")),
    # hlth_infor = factor(hlth_infor_m,
    #                          levels = c(1, 0, "Missing"),
    #                          labels = c("Yes", "No", "Missing"))
  )

# Table 1, unadjusted & RDS (the rest of 2022 vs during mpox) ----
## Table for proportions ----
# duplicate data to have outputs for 'overall' (both mpox visit and not mpox visit)
data_dbl <- bind_rows(
  data_22, # pool analysis does not include city
  mutate(data_22, 
         mpox_visit = "all")
)
nrow(data_dbl)
# generate variable of total nb of visits
data_dbl <- data_dbl %>% 
  group_by(mpox_visit) %>%
  mutate(mpox_visit_n = n(),
         mpox_visit_n_rds = sum(wt_rds_norm)) %>%
  ungroup()

# note the numbers are for visits, not number of people
tbl_1 <- data_dbl %>% 
  group_by(mpox_visit) %>%
  summarize(nb = n(), 
            .groups = "drop")

ls_var <- c("age_cats", 
            "rel_status", 
            "hiv_cats",
            "sa_cats_ttl",
            "sa_cats_anal")

# for each categorical variable compute the proportions
for(var in ls_var){
  tbl_1 <- tbl_1 %>% 
    bind_rows(
      data_dbl %>% 
        group_by(mpox_visit_n, mpox_visit_n_rds, mpox_visit, response = get(var)) %>% 
        summarize(nb = n(), 
                  nb_rds = sum(wt_rds_norm),
                  .groups = "drop_last") %>% 
        mutate(
          # unadjusted proportion
          prop = nb / mpox_visit_n,
          # RDS-adjusted
          prop_rds = nb_rds / mpox_visit_n_rds,
          # RDS 95% CI
          me = sqrt(prop_rds * (1 - prop_rds)) / sqrt(mpox_visit_n_rds) * 1.96,
          ci.lb = prop_rds - me,                                   
          ci.ub = prop_rds + me,
          # variable label
          char = var
        )
    )
}

## Table for mean number of partners ----
ls_vars_partn <- c("period_mpox", "month_intv", "nb_part_ttl", "nb_part_anal")
tbl_1_partn <- tibble()

for(var in ls_vars_partn){
  tbl_1_partn <- tbl_1_partn %>% 
    bind_rows(
      data_dbl %>%
        mutate(nb_part_rds = get(var) * wt_rds_norm) %>% 
        group_by(mpox_visit) %>%
        summarize(
          # unadjusted mean
          mean_nb = mean(get(var)),
          sd_nb = sprintf("(SD=%s)", round(sd(get(var)), 1)),
          # RDS-adjusted
          mpox_visit_n = n(),
          mpox_visit_n_rds = sum(wt_rds_norm),
          mean_rds = sum(nb_part_rds) / mpox_visit_n_rds,
          # RDS 95% CI
          me = sqrt( sum((nb_part_rds - mean_rds)^2) ) / mpox_visit_n_rds,
          ci.lb = mean_rds - me * 1.96,
          ci.ub = mean_rds + me * 1.96,
          .groups = "drop"
        ) %>%
        mutate(char = var)
    )
}

## Pivot and join tables ----
## proportions
# round proportions (need to turn into character to join with SD later)
tbl_1 <- tbl_1 %>% 
  mutate(across(c(prop, prop_rds, ci.lb, ci.ub), \(x) round_prop(x))) %>%
  mutate(prop = paste0("(", prop, "%)"),
         prop_rds = paste0(prop_rds, "%"),
         nb = as.integer(round(nb)))

# format for table
tbl_1 <- tbl_1 %>% 
  mutate(rds_ci = sprintf("(%s\u2013%s%s)", ci.lb, ci.ub, "%"))

# switch to one column per city
tbl_1 <- tbl_1 %>% 
  select(mpox_visit, nb, response, prop, prop_rds, char, rds_ci) %>% 
  pivot_wider(names_from = "mpox_visit", 
              values_from = c("nb", "prop", "prop_rds", "rds_ci"))

## number of partners
# round
tbl_1_partn <- tbl_1_partn %>% 
  mutate(across(c(mean_nb, mean_rds, ci.lb, ci.ub), \(x) round(x, 1)))

# turn RDS mean into character to join later
tbl_1_partn <- tbl_1_partn %>% 
  mutate(mean_rds = as.character(mean_rds))

# format for table
tbl_1_partn <- tbl_1_partn %>% 
  mutate(rds_ci = sprintf("(%s\u2013%s)", 
                          ci.lb, 
                          ci.ub))

# switch to column per mpox visit or not
tbl_1_partn <- tbl_1_partn %>% 
  select(mpox_visit, mean_nb, sd_nb, mean_rds, rds_ci, char) %>%
  pivot_wider(names_from = "mpox_visit",
              values_from = c("mean_nb", "sd_nb", "mean_rds", "rds_ci"))

# rename columns (match mean with counts and SD with proportions)
names(tbl_1_partn) <- gsub("mean_nb", "nb", names(tbl_1_partn))
names(tbl_1_partn) <- gsub("mean_rds", "prop_rds", names(tbl_1_partn))
names(tbl_1_partn) <- gsub("sd_nb_", "prop_", names(tbl_1_partn))

## join together and reorder columns
tbl_1 <- bind_rows(tbl_1, tbl_1_partn)

## Output table ----
tbl_1 <- tbl_1 %>% 
  select(char, 
         response,
         ends_with("_0"), 
         ends_with("_1"), 
         ends_with("all"))

write.csv(tbl_1, "./behavourial-tables/table_S3_unadj_and_rds.csv", 
          row.names = FALSE)