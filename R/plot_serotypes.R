dt <- read.table("..//denv_serotypes_Rio.csv", sep = "\t", header = TRUE)
matplot(1993:2010, dt, type = "l", lwd = 2, col = 1:3,
        lty = 1, xlab = "Time (years)", ylab = "Viral isolation")
legend(x="topleft", legend = c("DENV-1", "DENV-2", "DENV-3"), col = 1:3, lty = 1, bty = "n")
