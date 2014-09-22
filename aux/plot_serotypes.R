dt <- read.table("..//denv_serotypes_Rio.csv", sep = "\t", header = TRUE)
matplot(dt, type = "l", lwd = 2, col = 1:3, lty =1)
legend(x="topleft", legend = c("DENV-1", "DENV-2", "DENV-3"), col = 1:3, lty = 1, bty = "n")
