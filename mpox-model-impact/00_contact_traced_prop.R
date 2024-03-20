# authors:  Fanyu Xiu, Mathieu Maheu-Giroux
library(MpoxModelPack)

cal_cities <- c("mtl", "trt", "van")
traced_prop <- 0.20
tracing_dur <- 2 # days it takes for contact tracing 

for(cty in cal_cities){
 # loading reporting fraction in each city
init.pop.fn(cty, 1)
load.params.fn() 
# mtl  trt  van 
# 0.82 0.86 0.77 

x = seq(0, 15, 0.05)
y0 = dexp(x, 1/5.1) # distribution of how long people spend in the exposed compartment = latent period = 5.1 days (Xiu et al, 2024)
plot(y0 ~ x, type = "l") # pdf of the distribution 
y1 = pexp(x, 1/5.1, lower.tail = T) # cdf of the distribution: at a given x, y1 is the density of people stay for <=x days in the exposed compartment
y = 1 - y1 # ccdf of the distribution: at a given x, y is the density of people stay for >=x days in the exposed compartment
# or ccdf can also be computed by: y = pexp(x, 1/5.1, lower.tail = F)
plot(y ~ x, type = "l") # plot of the ccdf

# y[which(x == 4)] is the density of people who remains in the exposed compartment on the 4th days 
# i.e. the % of exposed people who hasn't become infectious on the 4th days in the exposed state
# 4 comes from 2 days for reporting delay plus contact tracing takes 2 days
contact_prop <- traced_prop * report_frac[[cty]] * y[which(x == 2 + tracing_dur)] # is the fraction of isolated among exposed 
# OR contact_prop <- traced_prop * report_frac[[cty]] * pexp(q = 2 + tracing_dur, rate = 1/5.1, lower.tail = F)
print(cty)
print(contact_prop)
}

