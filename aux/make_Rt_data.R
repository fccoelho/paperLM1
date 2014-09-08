### Code to export a CSV file with the disease cases and R_t plus credibility intervals.
### Rt is the raw estimate Y(t+1)/Y(t) and Rt2 is the mean of the inverse beta distribution
### Se the tech report for details: https://github.com/maxbiostat/CODE/blob/master/Alerta/confidence_Rt.pdf

source("Rt_aux.R")
source("elicit_Rt_prior.R")

dat <- read.csv("dadojuntos_201433.csv", header=TRUE)

N <- length(unique(dat$SE))
cases <- rowSums(matrix(dat$casos, nrow = N))
pars <- elicit.Rt.hyperpars(m0 = 1/2, v0 = 2)  
Rtb <- Rt.beta(cases, a0 = pars[1], b0 = pars[2])
dadoexp <- data.frame(cases, Rtb)
write.csv(dadoexp, file = "data_Rt_dengue.csv", row.names = FALSE)
