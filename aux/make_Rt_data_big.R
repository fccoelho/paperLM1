### Code to export a CSV file with the disease cases and R_t plus credibility intervals.
### Rt is the raw estimate Y(t+1)/Y(t) and Rt2 is the mean of the inverse beta distribution
### See the tech report for details: https://github.com/maxbiostat/CODE/blob/master/Alerta/confidence_Rt.pdf

source("Rt_aux.R")
source("elicit_Rt_prior.R")
library(lubridate)

dat <- read.csv("../dengue_cases_Rio_perEW_1996-2009.csv", sep = "\t", header = TRUE)
# please check these everytime you run the script with new data 
# use http://www.cve.saude.sp.gov.br/htm/Cve_se03.htm to get 'start' and 'end' right
start <- ymd("1996-01-08")
end <- ymd("2009-12-20")
starts <- format(seq(start, end, by = "week"), "%Y-%m-%d")
ends <- format(seq(as.Date(start) + 6, as.Date(end) + 6, by = "week"), "%Y-%m-%d")
seDates <- data.frame( start = starts, end = ends)
#  
cases <- dat$CASES
pars <- elicit.Rt.hyperpars(m0 = 1/2, v0 = 2)  
Rtb <- Rt.beta(cases, a0 = pars[1], b0 = pars[2])
dadoexp <- data.frame(SE = rep(1:52, 14), seDates, cases, Rtb)
write.csv(dadoexp, file = "data_Rt_dengue_big.csv", row.names = FALSE)
