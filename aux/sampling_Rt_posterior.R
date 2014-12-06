### Example code showing how to sample from posterior of R_t
### A raw estimate of R_t is Y(t+1)/Y(t) but we can also take the  mean of the (inverse beta) posterior distribution
### WARNING: please note that we know the posteriors exactly, so none of this random variate generation is actually needed.
### We are just illustrating how one can sample from the prior/posterior if one so desires.
source("Rt_aux.R")
source("elicit_Rt_prior.R")
set.seed(8632)
## getting the data
cases <- data.frame(read.table("../DATA/data_Rt_dengue.csv", sep = ",", header = TRUE))$cases

# now lets choose some values of cases
sorted <- sort(cases)
low <- sample(head(sorted, 10), 1)
medium <- sample(sorted[round(length(cases)/2):(round(length(cases)/2)+10)], 1)
high <- sample(tail(sorted, 10), 1)
values <- c(low, medium, high)
post <- sapply(values, function(x) match(x, cases)) # Y_{t}
post1 <- post + 1 # Y_{t+1}

(hatRt <- cases[post1]/cases[post])

## creating the priors

(pars.cons <- elicit.Rt.hyperpars(m0 = .5, c = .5)) # conservative prior
(pars.lax <- elicit.Rt.hyperpars(m0 = 1, c = 1)) # more lax prior
(pars.crazy <- elicit.Rt.hyperpars(m0 = 2, c = .2)) # crazy prior

# the actual priors
N <- 1000
cons.prior <- rRt(n = N, a = as.numeric(pars.cons["a_0"]), b = as.numeric(pars.cons["b_0"]))
lax.prior <- rRt(n = N, a = as.numeric(pars.lax["a_0"]), b = as.numeric(pars.lax["b_0"]))
crazy.prior <- rRt(n = N, a = as.numeric(pars.crazy["a_0"]), b = as.numeric(pars.crazy["b_0"]))

# now lets get the posteriors, shall we?
## with few cases
cons.post.low <- rRt(n = N, a = cases[post1][1] + as.numeric(pars.cons["a_0"]), b = cases[post][1] + as.numeric(pars.cons["b_0"]))
lax.post.low <- rRt(n = N, a = cases[post1][1] + as.numeric(pars.lax["a_0"]), b = cases[post][1] + as.numeric(pars.lax["b_0"]))
crazy.post.low <- rRt(n = N, a = cases[post1][1] +as.numeric(pars.crazy["a_0"]), b = cases[post][1] +as.numeric(pars.crazy["b_0"]))
## somewhat in the middle
cons.post.med <- rRt(n = N, a = cases[post1][2] + as.numeric(pars.cons["a_0"]), b = cases[post][2] + as.numeric(pars.cons["b_0"]))
lax.post.med <- rRt(n = N, a = cases[post1][2] + as.numeric(pars.lax["a_0"]), b = cases[post][2] + as.numeric(pars.lax["b_0"]))
crazy.post.med <- rRt(n = N, a = cases[post1][2] +as.numeric(pars.crazy["a_0"]), b = cases[post][2] +as.numeric(pars.crazy["b_0"]))
## big counts
cons.post.high <- rRt(n = N, a = cases[post1][3] + as.numeric(pars.cons["a_0"]), b = cases[post][3] + as.numeric(pars.cons["b_0"]))
lax.post.high <- rRt(n = N, a = cases[post1][3] + as.numeric(pars.lax["a_0"]), b = cases[post][3] + as.numeric(pars.lax["b_0"]))
crazy.post.high <- rRt(n = N, a = cases[post1][3] +as.numeric(pars.crazy["a_0"]), b = cases[post][3] +as.numeric(pars.crazy["b_0"]))

## Plotting
svg("../plots/prior_sensitivity_Rt.svg")
plot(density(cons.prior), xlim = c(0, 5), lwd = 2, col = 1, main = "Priors  & posteriors versus case counts", 
     xlab = expression(R[t]), ylab = "Density", lty = 1)
lines(density(lax.prior), lwd = 2, col = 2, lty = 1)
lines(density(crazy.prior), lwd = 2, col = 3, lty = 1)
#
lines(density(cons.post.low), col = 1, lwd = 2, lty = 2)
lines(density(lax.post.low), col = 2, lwd = 2, lty = 2)
lines(density(crazy.post.low), col = 3, lwd = 2, lty = 2)
#
lines(density(cons.post.med), col = 1, lwd = 2, lty = 3)
lines(density(lax.post.med), col = 2, lwd = 2, lty = 3)
lines(density(crazy.post.med), col = 3, lwd = 2, lty = 3)
#
lines(density(cons.post.high), col = 1, lwd = 2, lty = 4)
lines(density(lax.post.high), col = 2, lwd = 2, lty = 4)
lines(density(crazy.post.high), col = 3, lwd = 2, lty = 4)
legend(x = "topright", col = 1:3, legend = c("Conservative", "Diffuse", "Crazy") , lwd = 2, bty = "n", title = "Prior")
legend(x = "right", lty = 2:4, legend = c("Low", "Medium", "High") , lwd = 2, bty = "n", title = "Case counts")
dev.off()