#' @title init.pop.fn
#' @author: Fanyu Xiu, Carla Doyle
#' @description Sets up the initial model population of interest (according to city and type of sexual partners).
#' @param city "mtl" = Montreal, "trt" = Toronto, "van" = Vancouver (can input a vector), Default: c("mtl", "trt", "van")
#' @param type 1 = all-type, 2 = anal, Default: 1
#' @return pop_dat
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  init.pop.fn(city = c("trt", "van"), type = 1)
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{select}}
#' @rdname init.pop.fn
#' @export 

init.pop.fn <- function(city = c("mtl", "trt", "van"), type = 1) { 
  
  TYPE <- c("ttl", "anal")
  N.GBM <- c(54000, 78000, 26100)
  names(N.GBM) <- c("mtl", "trt", "van")
  
  pop_dat <- list(mtl = NULL,
                  trt = NULL,
                  van = NULL)
  
  name_age_cats <- sort(as.character(unique(contact_rate_prop$age_cats)))
  name_hiv_cats <- sort(as.character(unique(contact_rate_prop$hiv_cats)))
  name_sa_cats <- sort(as.character(unique(contact_rate_prop$sa_cats)))
  
  n_age_cats <- length(name_age_cats) 
  n_hiv_cats <- length(name_hiv_cats)
  n_sa_cats <- length(name_sa_cats)
  
  for(cur_city in city){
    
    type_name <- TYPE[type]
    
    city_df <-  dplyr::filter(contact_rate_prop, type == type_name & city == cur_city)
    city_df <-  dplyr::mutate(city_df, n_pop = prop * N.GBM[cur_city])
    city_df <- dplyr::arrange(city_df, factor(sa_cats, levels = name_sa_cats)) 
    city_df <- dplyr::select(city_df, - c("city", "prop"))
    
    for(par in c("n_pop", "c_ash")){
      par_index <- which(colnames(city_df) == par)
      
      pop_dat[[cur_city]][[par]] <- array(data = NA, 
                                           dim = c(n_age_cats, n_sa_cats, n_hiv_cats),
                                           dimnames = list(name_age_cats, name_sa_cats, name_hiv_cats))
      
      
      for(a in name_age_cats){
        for(s in name_sa_cats){
         for(h in name_hiv_cats){
        pop_dat[[cur_city]][[par]][a, s, h] <- unlist(unname(city_df[city_df$age_cats == a & city_df$hiv_cats == as.numeric(h) & city_df$sa_cats == s, 
                                                                      par_index]))
        }}}
    }}

  assign("pop_dat", pop_dat, envir = .GlobalEnv)

  return(pop_dat)
}