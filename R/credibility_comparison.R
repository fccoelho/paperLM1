## This script will compare the credibility intervals obtained using the methods of
## Nishiura et al. (2010) [slightly modified] and Coelho & Carvalho (2015)
# generate the data
set.seed(666)
lambdas <- expand.grid(l1 = seq(2,  1000, length.out = 10),
                       l2 = seq(2, 1000, length.out = 10))
counts <- apply(lambdas, 1, function(x){
  data.frame(Y1 = rpois(10000, lambda = x[1]),
             Y2 = rpois(10000, lambda = x[2]))
})

# calculate empirical "Rt"
getempRt <- function(dt){
  dtt <- subset(subset(dt, Y1 > 0), Y2 > 0)
  return(dtt$Y2/dtt$Y1)
}
empRts <- lapply(counts, getempRt)
# lapply(empRts, hist)

# calculating posterior means and credible intervals
source("Rt_aux.R")
nzcounts <- lapply(counts, function(d) subset(subset(d, Y1 > 0), Y2 > 0))

betas <- lapply(nzcounts,
                function(d){
                  apply(d, 1, 
                        function(x) estBeta(y1 = x[1], y2 = x[2]))})

gammas <- lapply(nzcounts,
                function(d){apply(d, 1,
                                  function(x) estGamma(y1 = x[1], y2 = x[2]))})

CIcomparisons <- vector(length(betas), mode = "list")
for (i in 1:length(betas)){
  indexes <- 1:length(betas[[i]])
  tempcompz <- sapply(indexes, function(j) compareCI(c1 = betas[[i]][[j]]$CI,
                                                           c2 = gammas[[i]][[j]]$CI))
  CIcomparisons[[i]] <- table(tempcompz)/sum(table(tempcompz))
}
CIcomparisons
## integrating over the grid
betas.total <- unlist(betas, recursive = FALSE)
gammas.total <- unlist(gammas, recursive = FALSE)
inds <- 1:length(betas.total)
compz.total <- sapply(inds, function(j) compareCI(c1 = betas.total[[j]]$CI,
                                            c2 = gammas.total[[j]]$CI))
prob.compz.total <- sapply(inds,
                           function(j) betas.total[[j]]$pmt1 < gammas.total[[j]]$pmt1)

table(prob.compz.total)/sum(table(prob.compz.total))
table(compz.total)/sum(table(compz.total))
