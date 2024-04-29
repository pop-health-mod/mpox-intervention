# author: Jesse Knight
# load data
source("./src/helper_mixing.R")
load("./parametrization-outputs/A_matrix.rda")
load("./parametrization-outputs/contact_rate_prop_smooth.rda")
X0 = contact_rate_prop; rm(contact_rate_prop)
A  = A_matrix;          rm(A_matrix)    # age mixing (probability) from Milwid2022
X0$city     = as.character(X0$city)     # fix factor issues
X0$age_cats = as.character(X0$age_cats) # fix factor issues
# define labels
lab = list( # strata labels
  y = unique(X0$city),
  a = unique(X0$age_cats),
  h = unique(X0$hiv_cats),
  s = unique(X0$sa_cats))
lab$ah = c(outer(lab$a, lab$h, paste)) # age:hiv
# data pre-processing & constants
X_ah = aggregate(cbind(x = prop * c_ash) ~ age_cats + hiv_cats + city,
                 X0,
                 sum) # x = total contacts
H_mtl = matrix(c(.924, .660, .076, .340), 2, 2) # hiv mixing (probability) from Milwid2022
H = list(mtl = H_mtl,
         trt = NA * H_mtl,
         van = NA * H_mtl)
N = list(mtl = 54000, 
         trt = 78000, 
         van = 26100) # population size (not actually needed)
# run fitting
mix_odds = list() # for each city
fix = 0 # fix odds[1] (hiv) after estimating from mtl
for (y in lab$y){ # cities
  x = N[[y]] * setNames(subset(X_ah, city == y)$x,
                        lab$ah)      # total contacts
  mix_odds[[y]] = fit.mix.ah(x = x,
                             A = A[[y]],
                             H = H[[y]],
                             fix = fix) # odds vector
  if (y == 'mtl'){ fix = mix_odds[[y]][1] }
  else         { mix_odds[[y]][1] = fix }
}
# save output
save(mix_odds,
     file = "./parametrization-outputs/mix_odds.rda")
