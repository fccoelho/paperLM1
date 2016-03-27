#############################
## obtain 100\% alpha confidence/credibility intervals for the success probability \theta
betaconf <- function(alpha = .95, x, n, a = 1, b = 1, CP = "FALSE"){
  if(CP=="TRUE"){  
    lower <- 1 - qbeta((1-alpha)/2, n + x - 1, x)
    upper <- 1 - qbeta((1+alpha)/2, n - x, x + 1)
  }else{
    lower <- qbeta((1-alpha)/2, a + x, b + n - x)
    upper <- qbeta((1+alpha)/2, a + x, b + n - x)  
  } 
  return(c(lower, upper))
  # CP stands for Clopper-Pearson
  # Default is 'Bayesian' with an uniform prior over p (a=b=1)
}
#############################
# R = theta/(1-theta)
odds <- function(x) x/(1-x)
#############################
Rt.beta <- function(d, gt = 3, alpha = .95, a0 = 2 , b0 = 3){
  ## Default a0 and b0 chosen so that E[R_t] = 1 and Var(R_t) = 2 
  ## gt = generation time, default 3 weeks
  ## alpha = *100% level of confidence
  N <- length(d)
  ac <- d[gt:N]
  for(i in 1:(gt-1)){
    ac <- ac + d[(gt-i):(N-i)] 
  } 
  K <- length(ac)
  yk1 <- ac[2:K] # Y[k+1] 
  yk  <- ac[1:(K-1)]# Y[k]  
  Rt <- Rt2 <- NA
  Rt[(gt + 1):N] <- yk1/yk
  Rt2[(gt + 1):N] <- (yk1 + a0)/(yk + b0 - 1)
  CIs <- data.frame (p1 = rep(NA, N), lwr = rep(NA, N), upr = rep(NA, N) )
  ## p1 = Pr(R>1) 
  for( k in 1: (N-gt)){
    CIs[k + gt, 1] <- 1 - pbeta(.5, shape1 = yk1[k] + a0 , shape2 = yk[k] + b0)
    CIs[k + gt, 2:3] <- odds(betaconf(alpha = alpha, x = yk1[k], 
                                  n = yk1[k] + yk[k], a = a0, b = b0))
  } 
  return(data.frame(Rt, Rt2, CIs))
}
###################################################
rRt <- function(n, a, b){ # let's sample from the distribution of R_t
  x <- rbeta(n, shape1 = a, shape2 = b)
  y <- x/(1-x) # yup! simple as that!
return(y)
}
###################################################
Rt.gamma <- function(d, gt = 3, alpha = .95){
  ## It's all quite simple: we have the conditional distribution of R_{t}:
  ## p(R_{t}| Y_{t}, Y_{t + 1}) = (R_{t}*Y_{t})^Y_{t + 1} * exp(- R_{t}*Y_{t})
  ## Nishiura et al (2010) -- JRSInterface
  ## Lets just sample analytically from it, since its a(n unnormalised) 
  ## Gamma distribution  
  ## gt = generation time, default 3 weeks
  ## alpha = *100% level of confidence
  N <- length(d)
  ac <- d[gt:N]
  for(i in 1:(gt-1)){
    ac <- ac + d[(gt-i):(N-i)] 
  } 
  K <- length(ac)
  yk1 <- ac[2:K] # Y[k+1] 
  yk  <- ac[1:(K-1)]# Y[k]  
  Rt <- NA
  Rt[(gt + 1):N] <- yk1/yk
    CIs <- data.frame (p1 = rep(NA, N), lwr = rep(NA, N), upr = rep(NA, N) )
     ## p1 = Pr(R>1) 
  for( k in 1: (N-gt)){
    CIs[k + gt, 1] <- 1-pgamma(1, shape = yk1[k] + 1, rate = yk[k] )
    CIs[k + gt, 2:3] <- qgamma(p = c((1-alpha)/2, (1+alpha)/2 ),
                             shape = yk1[k] + 1, rate = yk[k] )
  } 
  return(data.frame(Rt, CIs))
}
###
estBeta <- function(y1, y2, alpha = .95, a0 = 2, b0 = 3){
  mean <- (y2 + a0)/(y1 + b0 - 1)
  CI <- odds(betaconf(alpha = alpha, x = y2, 
                    n = y1 + y2, a = a0, b = b0 ))
  prob <- 1 - pbeta(.5, shape1 = y2 + a0, shape2 = y1 + b0)
  return(list(post.mean = mean, CI = CI, pmt1 = prob))
}
estGamma <- function(y1, y2, alpha = .95){
  mean <- y2/y1
  prob <- 1-pgamma(1, shape = y2 + 1, rate = y1)
  CI <- qgamma(p = c((1-alpha)/2, (1+alpha)/2 ),
               shape = y2 + 1, rate = y1)
  return(list(mean = mean, CI = CI, pmt1 = prob))
}
compareCI <- function(c1, c2){
  ## CIs are expected in the format c(lower, upper)
  if(length(c1) != 2 || length(c2) != 2) stop("Please supply 2-dimensional vectors")
  res <- NA
  test <- 0
  if(c1[1] < c2[1]) test <- test + 1
  if(c1[2] > c2[1]) test <- test + 2
  # test = 0 : c1 is contained in c2 
  # test = 1 or 2 : there is overlap
  # test = 3 : c1 contains c2
  return(test)
}
