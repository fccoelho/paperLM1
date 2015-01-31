## Code to obtain the hyperparamenters for the prior for R_{t}
## Copyleft (or the one to blame): Carvalho, LMF (2014)
## You can supply mean and variance (m0 and v0) or
## mean and coefficient of variation (m0 and c)
elicit.Rt.hyperpars <- function(m0, v0 = NULL, c = .5){
  if(is.null(v0)){
    a0 <- (m0^3*c^2 + m0^3 + m0^2)/(m0^2*c^2) 
    b0 <- (2*m0^2*c^2 + m0^2 + m0)/(m0^2*c^2)
  } else{
    a0 <- (m0*v0 + m0^3 + m0^2)/v0
    b0 <- (2*v0 + m0^2 + m0)/v0
  }  
  pars <- c(a0, b0)
  names(pars) <- c("a_0", "b_0")
return(pars)
}
