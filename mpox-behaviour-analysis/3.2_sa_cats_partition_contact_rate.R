# author: Fanyu Xiu, Jorge Luis Flores Anato
## library and data --------
library(dplyr)
library(readr)
library(tidyr)
source("./src/helper_dates.R")
load("./data-misc/age_matrix.RData")
# age_prop_census <- read.csv("./data-misc/age_prop_census.csv")
# age_prop_census_fine <- read.csv("./data-misc/age_prop_census_fine.csv")
data_pre_mpox <- read.csv("./data-3cities-feb-2023/data_contact_rate_3cities.csv")
outcome_var <- "nb_part_ttl"

ifelse(outcome_var == "nb_part_ttl",
       data_pmf <- read_csv("./parametrization-outputs/fitted-distr-ttl/pmf_weighted_all.csv"),
       data_pmf <- read_csv("./parametrization-outputs/fitted-distr-anal/pmf_weighted_anal.csv"))

type <- ifelse(outcome_var == "nb_part_ttl",
               "ttl",
               "anal")

# contact rate-------
data_pmf <- data_pmf %>% 
  filter(grepl("2022-Pre-Mpox", data_pt))

AGES <- sort(unique(data_pmf$age_cats))
HIV <- sort(as.character(unique(data_pmf$hiv_cats)))

# categorize based on groups
df_results <- create_city_list()

for(cty in c("mtl", "trt", "van")){
  
  cdf_cutoff <- c(.30, .40, .50, .60, .70, 
                  .80, .90, .95, .96, .97, 
                  .98, .99, .995, .998, 1) 
  
  df_results[[cty]] <- vector("list", length(AGES))
  names(df_results[[cty]]) <- AGES
  
  for(a in AGES){
    
    df_results[[cty]][[a]] <- vector("list", length(HIV))
    names(df_results[[cty]][[a]]) <- HIV
    
    for(h in HIV){
    data_city_ah <- data_pmf %>% 
      filter(
        grepl(case_when(cty == "mtl" ~ "Montreal",
                        cty == "trt" ~ "Toronto",
                        cty == "van" ~ "Vancouver"),
              data_pt)
      ) %>%
      filter(age_cats == a & hiv_cats == h)
    
    # fitted distributions of sexual partnerships
    pmf_size <- data_city_ah$mean
    pmf_rate <- data_city_ah$y_pred
    
    # sexual activity groups to create
    k_size <- vector("double", 0)
    k_rate <- vector("double", 0)
    
    # which points of the PMF are included
    k_min <- vector("double", 0)
    k_max <- vector("double", 0)
    
    # correct the PMF by reweighing to 1
    if( sum(pmf_size) < 1 ){
      print(sprintf("%s %s %s pmf was reweighted, previously added up to %s", cty, a, h, sum(pmf_size)))
      pmf_size <- pmf_size / sum(pmf_size)
    }
    
    for(i in 1:length(cdf_cutoff)){
      # print(i)
      # verify before running loop that the PMF adds up to 1
      if( i == 1 & abs(sum(pmf_size) - 1) > 0.0001) { stop("error; sum(PMF) < 1") }
      
      # take x% from the current group
      x_perc <- ifelse(i == 1,
                       cdf_cutoff[i],
                       cdf_cutoff[i] - cdf_cutoff[i - 1])
      
      if(pmf_size[1] >= x_perc){
        ## if we have enough people we take the required density, assign rate and go to next group
        k_size[i] <- min(x_perc, pmf_size[1])
        pmf_size[1] <- pmf_size[1] - k_size[i]
        
        k_rate[i] <- pmf_rate[1]
        
        # store min and max
        k_min[i] <- k_max[i] <- pmf_rate[1]
      } else {
        ## otherwise pull from additional groups
        ## and get a rate that is the weighted average of the partnerships in all groups
        vec_rate <- vector("double", 0)
        vec_size <- vector("double", 0)
        
        # tracks how much we have left to assign
        x_perc_left <- x_perc
        j <- 1
        
        while(round(sum(vec_size), 10) < round(x_perc, 10)){ ## seems to hang here
          # first note the rate and the contribution of the group
          vec_size[j] <- min(x_perc_left, pmf_size[1])
          pmf_size[1] <- pmf_size[1] - vec_size[j]
          
          vec_rate[j] <- pmf_rate[1]
          
          # check progress
          # if(pmf_rate[1] %% 50 == 0) { print(pmf_rate[1]) }
          
          # then remove PMF point from vector
          if(pmf_size[1] <= 0 & length(pmf_size) > 1){
            pmf_size <- pmf_size[2:length(pmf_size)]
            pmf_rate <- pmf_rate[2:length(pmf_rate)]
          }
          
          # note how much we have added to vec_size
          x_perc_left <- x_perc_left - vec_size[j]
          
          # increase counter
          j <- j + 1
        }
        
        # now we add the size and weighted rate to the group compartment
        k_size[i] <- sum(vec_size)
        k_rate[i] <- sum(vec_size * vec_rate) / sum(vec_size)
        
        # verify that the weighted average makes sense
        if( !(vec_rate[1] <= k_rate[i] & k_rate[i] <= vec_rate[j - 1]) ) { print("error; weighted rate calculation") }
        
        # store min and max
        k_min[i] <- vec_rate[1]
        k_max[i] <- vec_rate[j-1]
      }
      
      # remove PMF point from vectors if we have taken all the density
      if(pmf_size[1] <= 0 & length(pmf_size) > 1){
        pmf_size <- pmf_size[2:length(pmf_size)]
        pmf_rate <- pmf_rate[2:length(pmf_rate)]
      }
    }
    
    ## verify results
    # verify group sizes
    
    if( abs(sum(k_size) - 1) > 0.0001 ) { print(sprintf("%s %s final group sizes do not add up to 1, added up to %s", cty, ah, sum(k_size))) }
    
    # verify mean partnership rate
    sum(k_rate * k_size) / sum(k_size)
    sum(data_city_ah$y_pred * data_city_ah$mean) / sum(data_city_ah$mean)
    
    # format results
    df_out <- data.frame(city = cty,
                         age_cats = a,
                         hiv_cats = h,
                         sa_cats = paste("group_", 1:length(cdf_cutoff), sep = ""),
                         mean_rate = k_rate,
                         prop = k_size,
                         min_grp = k_min,
                         max_grp = k_max)
    
    df_results[[cty]][[a]][[h]] <- df_out
  } # for hiv
    df_results[[cty]][[a]] <- bind_rows(df_results[[cty]][[a]])
  } # for age
    df_results[[cty]] <- bind_rows(df_results[[cty]])
}# for city

df_results <- bind_rows(df_results) 
df_results$prop; round(df_results$prop, 4)

df_results$prop <- round(df_results$prop, 5)

df_results %>% 
  group_by(age_cats, hiv_cats, city) %>% 
  summarize(sum_prop = sum(prop))

### multiply by proportion of each age-hiv combination in each city
### method 1: use Engage weighted age-hiv proportion
df_ah_freq <- data_pre_mpox %>%
  mutate(hiv_cats = as.character(hiv_cats)) %>%
  group_by(city, age_cats, hiv_cats) %>%
  summarize(n_cah = sum(ipw_rds)) %>%
  ungroup() %>%
  complete(city, age_cats, hiv_cats, fill = list(n_cah = 0)) %>%
  group_by(city) %>%
  mutate(freq = n_cah / sum(n_cah)) %>%
  ungroup() %>%
  dplyr::select( - n_cah)

### method 2: use census age proportion, times the Engage weighted hiv proportion for each age group
# df_ah_freq <- data_pre_mpox %>% 
#   mutate(hiv_cats = as.character(hiv_cats)) %>% 
#   group_by(city, age_cats, hiv_cats) %>% 
#   summarize(n_cah = sum(ipw_rds)) %>%
#   ungroup() %>% 
#   complete(city, age_cats, hiv_cats, fill = list(n_cah = 0)) %>% 
#   group_by(city, age_cats) %>% 
#   mutate(hiv_prop_age = n_cah / sum(n_cah)) %>% 
#   ungroup() %>%
#   dplyr::select( - n_cah)
# 
# df_ah_freq <- df_ah_freq %>% 
#   left_join(filter(age_prop_census, age_cats != "ttl"),
#                         by = c("city", "age_cats")) %>%
#   dplyr::select( - pop_size) %>% 
#   mutate(freq = hiv_prop_age * age_prop)
# end method 2 

### check they sum to 1 
df_ah_freq %>% 
  group_by(city) %>% 
  summarize(sum_freq = sum(freq))

### multiplying the sexual activity strata proportion for each age-hiv combination
df_results <- df_results %>% 
  left_join(df_ah_freq, by = c("city", "age_cats", "hiv_cats")) %>%
  mutate(prop = prop * freq ) %>%
  dplyr::select(- freq)

### check they sum to 1 
df_results %>% 
  group_by(city) %>% 
  summarize(sum_prop = sum(prop))

write_csv(df_results, "./parametrization-outputs/fitted_distr_partn_p6m_standard_size.csv")
  
## store data in correct format to data folder ----
contact_rate_prop <- df_results %>%
  mutate(type = type,
         c_ash = mean_rate / 180,
         sa_cats = ifelse(parse_number(sa_cats) >= 10, sa_cats, paste0("group_0", parse_number(sa_cats)))) 
save(contact_rate_prop, file = "./parametrization-outputs/contact_rate_prop.rda")

# city-specific age matrix -------
age_cat_name <- sort(unique(as.character(data_pre_mpox$age_cats)))

## method 1: use age frequency from the Engage
df_a_freq <- data_pre_mpox %>%
  group_by(city, age_cats) %>%
  summarize(n_ca = sum(ipw_rds)) %>%
  ungroup() %>%
  complete(city, age_cats, fill = list(n_ca = 0)) %>%
  group_by(city) %>%
  mutate(freq = n_ca / sum(n_ca)) %>%
  ungroup() %>%
  dplyr::select( - n_ca)
## freq of each age category in each city according to the original A matrix
df_ori_a_freq <- data_pre_mpox %>%
  mutate(ori_age_cats = case_when(
    age %in% 16:19 ~ "16-19",
    age %in% 20:24 ~ "20-24",
    age %in% 25:29 ~ "25-29",
    age %in% 30:34 ~ "30-34",
    age %in% 35:39 ~ "35-39",
    age %in% 40:44 ~ "40-44",
    age %in% 45:49 ~ "45-49",
    age %in% 50:54 ~ "50-54",
    age %in% 55:59 ~ "55-59",
    age %in% 60:64 ~ "60-64",
    age %in% 65:69 ~ "65-69",
    age > 69 ~ "70+",
  )) %>% group_by(city, ori_age_cats) %>%
  summarize(n_ori_ca = sum(ipw_rds)) %>%
  ungroup() %>%
  complete(city, ori_age_cats, fill = list(n_ori_ca = 0)) %>%
  group_by(city) %>%
  mutate(freq = n_ori_ca / sum(n_ori_ca)) %>%
  ungroup() %>%
  dplyr::select( - n_ori_ca)
df_a_freq %>% group_by(city) %>% summarize(sum(freq))

### note that there is no "16-19" age groups in all cities so manually  add those
df_add <- data.frame(city = CITIES, ori_age_cats = rep("16-19", 3), freq = rep(0, 3))
df_ori_a_freq <- rbind(df_ori_a_freq, df_add) %>%
  arrange(ori_age_cats)

# ## method 2: use census age freq
# df_a_freq <- filter(age_prop_census, age_cats != "ttl") %>% 
#   rename(freq = age_prop)
# df_ori_a_freq <- filter(age_prop_census_fine, age_cats != "ttl") %>% 
#   rename(freq = age_prop,
#          ori_age_cats = age_cats)
# 
# df_ori_a_freq %>% group_by(city) %>% summarize(sum(freq))

### original matrix from the ref (12 age group)
A_matrix_origin <- age_matrix
ori_age_cat_name <- sort(unique(df_ori_a_freq$ori_age_cats))
dimnames(A_matrix_origin) <- rep(list(ori_age_cat_name), 2)

### create new age matrix
A_matrix <- list()

for(cur_city in CITIES){
  
  ### create the new 5*5 age matrix for each city
  A_matrix[[cur_city]] <- matrix(nrow = 5, ncol = 5)
  dimnames(A_matrix[[cur_city]]) <- rep(list(age_cat_name), 2)
  
  ### derive the data set for each city
  city_a_freq <- df_a_freq %>% 
    filter(city == cur_city)
  city_ori_a_freq <- df_ori_a_freq %>% 
    filter(city == cur_city)

  for(a_new_i in age_cat_name){
    
    denom = city_a_freq[city_a_freq$age_cats == a_new_i, ]$freq # denominator is the same across one row
    
    ### a vector indicating which of the original groups are encompassed by
    ### the ith new age group
    if(a_new_i == "60+"){
    bd_i <- 60:70
    }else{
    lower_bd_i <- as.numeric(unlist(strsplit(a_new_i,'-'))[1])
    upper_bd_i <- as.numeric(unlist(strsplit(a_new_i,'-'))[2])
    bd_i <- lower_bd_i:upper_bd_i
    }
    
    ori_age_cat_included_i <- ori_age_cat_name[parse_number(ori_age_cat_name) %in% bd_i]
      
    for(a_new_j in age_cat_name){
      
      numer = 0 # numerator is specific to each cell (i.e. varies across both rows and cols
      
      if(a_new_j == "60+"){
        bd_j <- 60:70
      }else{
        lower_bd_j <- as.numeric(unlist(strsplit(a_new_j,'-'))[1])
        upper_bd_j <- as.numeric(unlist(strsplit(a_new_j,'-'))[2])
        bd_j <- lower_bd_j:upper_bd_j
      }
      
      ori_age_cat_included_j <- ori_age_cat_name[parse_number(ori_age_cat_name) %in% bd_j]
      
      for(a_ori_i in ori_age_cat_included_i){
      for(a_ori_j in ori_age_cat_included_j){
          numer = numer + sum(A_matrix_origin[a_ori_i, a_ori_j] * city_ori_a_freq[city_ori_a_freq$ori_age_cats == a_ori_i, ]$freq)
      }}
      
    A_matrix[[cur_city]][a_new_i, a_new_j] <- numer / denom
  } ### for new matrix column
  } ### for new matrix row
  print(rowSums(A_matrix[[cur_city]]))
  } ### for city

save(A_matrix, file = "./parametrization-outputs/A_matrix.rda")

# case data -----
case_data <- read.csv("./data-misc/mpox_case_data_phac.csv")
save(case_data, file = "./parametrization-outputs/case_data.rda")