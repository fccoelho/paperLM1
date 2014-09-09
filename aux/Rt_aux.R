#############################
## obtain 100\lapha confidence/credibility intervals for the success probability \theta
betaconf <- function(alpha = .95, x, n, a = 1, b = 1, CP = "FALSE"){
  if(CP=="TRUE"){  
    lower <- 1 - qbeta((1-alpha)/2, n + x - 1, x)
    upper <- 1 - qbeta((1+alpha)/2, n - x, x + 1)
  }else{
    lower <- qbeta( (1-alpha)/2, a + x, b + n - x)
    upper <- qbeta(( 1+alpha)/2, a + x, b + n - x)  
  } 
  return(c(lower, upper))
  #CP stands for Clopper-Pearson
  #Default is 'Bayesian' with an uniform prior over p (a=b=1)
}
#############################
# R = theta/(1-theta)
ll <- function(x) x/(1-x)
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
  jk1 <- ac[2:K] # J[k+1] 
  jk  <- ac[1:(K-1)]# J[k]  
  Rt <- NA
  Rt[(gt + 1):N] <- jk1/jk
  CIs <- data.frame (p1 = rep(NA, N), lwr = rep(NA, N), upr = rep(NA, N) )
  ## p1 = Pr(R>1) 
  for( k in 1: (N-gt)){
    CIs[k+gt, 1] <- 1 - pbeta(.5, shape1 = jk1[k], shape2 = jk[k])
    CIs[k+gt, 2:3] <- ll(betaconf(alpha = alpha, x = jk1[k], 
                                  n = jk1[k] + jk[k], a = a0, b = b0 ))
  } 
  return(data.frame(Rt, CIs))
}
###################################################
Rt.gamma <- function(d, gt = 3, alpha = .95){
  ## It's all quite simple: we have the conditional distribution of R_{t}:
  ## p(R_{t}| J_{t}, J_{t + 1}) = (R_{t}*J_{t})^J_{t + 1} * exp(- R_{t}*J_{t})
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
  jk1 <- ac[2:K] # J[k+1] 
  jk  <- ac[1:(K-1)]# J[k]  
  Rt <- NA
  Rt[(gt + 1):N] <- jk1/jk
    CIs <- data.frame (p1 = rep(NA, N), lwr = rep(NA, N), upr = rep(NA, N) )
     ## p1 = Pr(R>1) 
  for( k in 1: (N-gt)){
    CIs[k+gt, 1] <- 1-pgamma(1, shape = jk1[k] + 1, rate = jk[k] )
    CIs[k+gt, 2:3] <- qgamma(p = c((1-alpha)/2, (1+alpha)/2 ),
                             shape = jk1[k] + 1, rate = jk[k] )
  } 
  return(data.frame(Rt, CIs))
}
