## author: Fanyu Xiu, Jorge Luis Flores Anato
## library and data --------
library(lubridate)

## cities
CITIES <- c("mtl", "trt", "van")
PROV <- c("QC", "ON", "BC")
create_city_list <- function(list_names = CITIES){
  ls_new <- vector("list", length = length(list_names))
  names(ls_new) <- list_names
  return(ls_new)
}

make_ind_age <- function(data){
  data %>% 
    mutate(age_30_39 = as.integer(age_cats == "30-39"),
           age_40_49 = as.integer(age_cats == "40-49"),
           age_50_59 = as.integer(age_cats == "50-59"),
           age_60_   = as.integer(age_cats == "60+"))
}

# partnership status; reference is not partnered
make_ind_rel <- function(data){
  data %>% 
    mutate(rel_status = ifelse(is.na(rel_status), "unclear", rel_status),
                  rel_status = factor(rel_status, levels = c("single", "open", "exclusive", "unclear"))) %>%
    mutate(rel_y_excl = as.integer(rel_status == "exclusive"),
           rel_y_open = as.integer(rel_status == "open"),
           rel_y_uncl = as.integer(rel_status == "unclear"))
}

# summarize effective sample size from Stan coefficients in a tibble format
summarize_ess <- function(df, beta_only = FALSE){
  # display ESS only for the beta coefficients
  if(beta_only){
    df_coeff_names <- df$coeff
    df <- df[grep("beta", df_coeff_names), ]
  }
  
  # compute the summary stats for the ESS
  df <- df %>% 
    summarize(mean = mean(n_eff), median = median(n_eff), min = min(n_eff), max = max(n_eff),
              .groups = "drop")
  
  return(df)
}

## definition of key dates
#### use first paper's start point of post-restrictions period instead
DATE_RES_END_TEXT <- "Ease of COVID-19-related travel restrictions"
DATE_RES_END <- as.Date("2021-09-07")  # new measures for fully vaccinated international travelers to Canada came into force

DATE_MPOX_START_TEXT <- "First case of Mpox in Canada"
DATE_MPOX_START <- as.Date("2022-05-19") 

### mpox peak period "end" when PrEP vaccination was scaled up in Montreal
DATE_MPOX_VA_TEXT <- "Mpox PrEP vaccination expanded in Montreal"
DATE_MPOX_VA <- as.Date("2022-06-14") 

#### gives 304 participants for restriction
DATE_MPOX_VA2_TEXT <- "Aug 14th - 2 months after first doses were expanded"
DATE_MPOX_VA2 <-  as.Date("2022-08-14")

# for table 1
round_prop <- function(x){
  as.character(round(x * 100, 0))
}