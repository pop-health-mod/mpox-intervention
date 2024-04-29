
city = "van"
import_cases_city = 4
bbeta_city = 0.9
omega_city = 6
RR_H_city = 0.94
RR_L_city = 0.95
gamma1_city = 1/7
period_city = 160
sd_city = 2
TRACING = 1
VACCINATING = 1
cpp = T
init.pop.fn()
load.params.fn()
resultc <- fn_model(city,
                    import_cases_city,
                    bbeta_city,
                    omega_city,
                    RR_H_city,
                    RR_L_city,
                    gamma1_city,
                    period_city,
                    sd_city,
                    TRACING,
                    VACCINATING,
                    cpp
)

cases <- resultc$cases[resultc$time > days_imported[["van"]]]
time <- seq(0, (length(cases) - 1) * 0.25 , 0.25)
plot(cases ~ time,
     main = city,
     ylim = c(0, 10))
points(case_data[case_data$prov == "BC", ]$city_cases)

city = "trt"
import_cases_city = 5.40
bbeta_city = 0.91
omega_city = 5
RR_H_city = 0.94
RR_L_city = 0.95
gamma1_city = 1/5.76
period_city = 170
sd_city = 2
TRACING = 1
VACCINATING = 1
cpp = T
init.pop.fn()
load.params.fn()
resultc <- fn_model(city,
                    import_cases_city,
                    bbeta_city,
                    omega_city, 
                    sd_city,
                    RR_H_city,
                    RR_L_city,
                    gamma1_city,
                    period_city,
                  
                    TRACING,
                    VACCINATING,
                    cpp
)

cases <- resultc$cases[resultc$time > days_imported[["trt"]]]
time <- seq(0, (length(cases) - 1) * 0.25 , 0.25)
plot(cases ~ time,
     main = city,
     ylim = c(0, 20))
points(case_data[case_data$prov == "ON", ]$city_cases)





city = "mtl"
import_cases_city = 4
bbeta_city = 0.9
omega_city = 6
RR_H_city = 0.94
RR_L_city = 0.95
gamma1_city = 1/6
period_city = 150
sd_city = 2
TRACING = 1
VACCINATING = 1
cpp = T
init.pop.fn()
load.params.fn()
resultc <- fn_model(city,
                    import_cases_city,
                    bbeta_city,
                    omega_city,
                    sd_city,
                    RR_H_city,
                    RR_L_city,
                    gamma1_city,
                    period_city,
                    TRACING,
                    VACCINATING,
                    cpp
)

cases <- resultc$cases[resultc$time > days_imported[["mtl"]]]
time <- seq(0, (length(cases) - 1) * 0.25 , 0.25)
plot(cases ~ time,
     main = city,
     ylim = c(0, 10))
points(case_data[case_data$prov == "QC", ]$city_cases)
