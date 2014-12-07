### Code to calculate and plot R_t plus credibility intervals.
### Rt is the raw estimate Y(t+1)/Y(t) and Rt2 is the mean of the inverse beta distribution
### See the tech report for details: https://github.com/maxbiostat/CODE/blob/master/Alerta/confidence_Rt.pdf

source("Rt_aux.R")
source("elicit_Rt_prior.R")
library(lubridate)
#  
dadosexp <- data.frame(read.table("../DATA/data_Rt_dengue_complete.csv", sep = ",", header = TRUE))
dadoexp <- na.omit(dadosexp)
library(ggplot2)
library(gridExtra)
dadoexp$start <- as.Date(dadoexp$start)
svg("../plots/cases_and_Rt.svg")
par(mar = c(0, 0, 0, 0))
bottom.panel <- ggplot(dadoexp, aes(x = start, y = Rt2, ymin = lwr, ymax = upr)) +
    geom_line() +
  geom_ribbon(alpha = 0.25) + 
  geom_line(data = dadoexp, size = 1.0, colour = "blue" ) +
  labs(x = "Time (weeks)",
       y = expression(R[t]),
#        y = "Effective reproductive number",
       title = "")
top.panel <- ggplot(dadoexp, aes(x = start, y = cases)) +
  geom_line() +
  geom_line(data = dadoexp, size = 1.0, colour = "red" ) +
  labs(x = "Time (weeks)",
              y = "Dengue Cases",
       title = "")
grid.arrange(top.panel, bottom.panel, heights = c(1/2, 1/2)) 
dev.off()