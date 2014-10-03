### Code to export a CSV file with the disease cases and R_t plus credibility intervals.
### Rt is the raw estimate Y(t+1)/Y(t) and Rt2 is the mean of the inverse beta distribution
### See the tech report for details: https://github.com/maxbiostat/CODE/blob/master/Alerta/confidence_Rt.pdf

source("Rt_aux.R")
source("elicit_Rt_prior.R")
library(lubridate)

dat <- read.csv("../raw_data_Rt_dengue_complete.csv", header = TRUE)
pars <- elicit.Rt.hyperpars(m0 = 1/2, v0 = 2)  
Rtb <- Rt.beta(dat$cases, a0 = pars[1], b0 = pars[2])
dadoexp <- data.frame(dat, Rtb)
write.csv(dadoexp, file = "data_Rt_dengue_complete.csv", row.names = FALSE)
