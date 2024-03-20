# mpox-model-impact
 Code to calibrate the model in package MpoxModelPack with daily reported mpox cases and compute averted fraction of new infections (AF) from interventions. 


## data source (for calibration)
I used surveillance reports by the Direction régionale de santé publique de Montréal(https://www.acpjournals.org/doi/10.7326/M22-2699), Public Health Ontario (https://www.publichealthontario.ca/-/media/Documents/M/2022/monkeypox-episummary.pdf?rev=ccdc118970104a4c9e634eb46e52839c&sc_lang=en), and British Columbia Centre for Disease Control (http://www.bccdc.ca/Health-Info-Site/Documents/Monkeypox/Epidemiological_Summary/Mpox_Surveillance_20230109.pdf) to infer the fraction of city over provincial cases that varied weekly. I obtained the city-level cases by multiplying the provincial cases (https://health-infobase.canada.ca/mpox/) and the city-to-province fraction.

## code structure
- [`00_contact_traced_prop.R`](00_contact_traced_prop.R) estimation of proportion of contact traced and isolated among exposed for each city.
- [`01_run_calibration.R`](01_run_calibration.R) model calibration using Bayesian sampling importance resampling.
- [`02_run_calibration_sen_RR.R`](02_run_calibration_sen_RR.R) model calibration for sensitivity analysis of fixing RR to estimates from emipirical analysis.
- [`03_tabulate_result.R`](03_tabulate_result.R) tabulate calibrated output and AF with their 95% credible intervals.
- [`04_plot_fit_xvalid.R`](04_plot_fit_xvalid.R) plot model fit and proportion of cumulative cases among each age group and HIV-serostatus using calibrated parameters as compared to observed cases from surveillance data.
- [`05_plot_AF.R`](05_plot_AF.R) plot counterfactual senarios as compared to the senario with all interventions present.

