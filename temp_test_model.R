city = "van"
import_cases_city = 5
bbeta_city = 0.9
omega_city = 6
RR_H_city = 0.94
RR_L_city = 0.95
gamma1_city = 1/8
period_city = 150
sd = 2
TRACING = 1
VACCINATING = 1
cpp = F
init.pop.fn()
load.params.fn()
resultr <- fn_model(city,
                    import_cases_city,
                    bbeta_city,
                    omega_city,
                    RR_H_city,
                    RR_L_city,
                    gamma1_city,
                    period_city,
                    sd,
                    TRACING,
                    VACCINATING,
                    cpp
)

plot(resultr$cases ~ resultr$time,
     main = city,
     ylim = c(0, 10))
points(case_data[case_data$prov == "BC", ]$city_cases)


city = "trt"
import_cases_city = 5
bbeta_city = 0.9
omega_city = 5
RR_H_city = 0.94
RR_L_city = 0.95
gamma1_city = 1/8
period_city = 150
sd = 2
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
                    sd,
                    TRACING,
                    VACCINATING,
                    cpp
)

plot(resultc$cases ~ resultc$time,
     main = city,
     ylim = c(0,20))
points(case_data[case_data$prov == "ON", ]$city_cases)





city = "mtl"
import_cases_city = 3
bbeta_city = 0.90
omega_city = 3.4
RR_H_city = 0.94
RR_L_city = 0.95
gamma1_city = 1/9
period_city = 150
sd = 2.7
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
                    sd,
                    TRACING,
                    VACCINATING,
                    cpp
)

plot(resultc$cases ~ resultc$time,
     main = city,
     ylim = c(0,10))
points(case_data[case_data$prov == "QC", ]$city_cases)
