#' Aggregate-level socio-demographic and behavioural data derived from the Engage Cohort Study ( 2017 - 2023 )
#'
#' @description Engage Cohort Study is a prospective, 
#' population-based cohort study of GBM in Montréal, Toronto, and Vancouver. 
#' Eligible participants were self-identified cis or trans men living in one of the three cities,
#' aged ≥16 years, who reported sex with another man in the past 6 months (P6M), 
#' understood English or French, and provided written consent. From February 2017 to August 2019, 
#' participants were recruited using respondent-driven sampling (RDS), 
#' a method used to sample hard-to-reach populations to estimate more representative population characteristics.
#' Initial participants were purposively selected to represent diverse characteristics of the GBM community, 
#' and all participants were invited to recruit up to six peers in their social networks.
#'
#' @format The data frame has 450 rows and the following 10 columns:
#' \describe{
#'   \item{city}{name of the city, Montreal, Toronto, or Vancouver}
#'   \item{age_cats}{age groups}
#'   \item{hiv_cats}{HIV-serostatus, 0: seronegative/unknown, 1: seropositive}
#'   \item{sa_cats}{sexual acitivty groups acertained by on a type (ttl/anal) of sexual partners}
#'   \item{mean_rate}{contact rate per past 6 months before the mpox outbreak for age group a, sexual activity group s, and HIV-serostatus h}
#'   \item{prop}{proportion of the size of an age-sexual activity-HIV-serostatus group relative to the MSM population size in a city for a type (ttl or anal) of sexual partners}
#'   \item{min_grp}{minimum number of sexual partner numbers in the past 6 months of a sexual activity group}
#'   \item{max_grp}{maximum number of sexual partner numbers in the past 6 months of a sexual activity group}
#'   \item{type}{ttl: all-type sexual partners, anal: anal sexual partners}
#'   \item{c_ash}{contact rate per day before the mpox outbreak for age group a, sexual activity group s, and HIV-serostatus h}
#' }
"contact_rate_prop"


#' Age mixing matrix based on Milwid et al, 2022 and the Engage Cohort Study age distribution (2022)
#'
#' @description Probability of sexual partnership formation between gay, bisexual, and other men who have sex with men belonging to different age groups.
#' @format The list has 3 sublists (for Montreal, Toronto, and Vancouver), each contains a 5*5 matrix.
"A_matrix"


#' Daily confirmed mpox cases in Montréal, Toronto, and Vancouver (May-Oct 2022.
#'
#' @description case data derived from the Public Health Agency of Canada and city surveillance reports
#' @format The data frame has 513 rows and the following 9 columns:
#' \describe{
#'   \item{prov}{name of the provinces: Quebec, Ontario, and British Columbia, where Montreal, Toronto, and Vancouver are located, respectively}
#'   \item{date}{date}
#'   \item{time_conti}{day since first reported case in each city}
#'   \item{prov_incidence}{provincial daily mpox incidence}
#'   \item{prov_cumul_cases}{provincial cumulative mpox incidence}
#'   \item{city_prov}{city fraction of provincial mpox incidence, varying weekly}
#'   \item{confirmed_total}{confirmed fraction of confirmed + probable mpox incidence}
#'   \item{first_dose}{mpox vaccine doses administered, varying weekly and allocated uniformly to each day}
#' }
"case_data"