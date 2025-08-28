## Forest plot for MDD

library(metafor)

#bart
dat <- data.frame(tau_hat = c(0.852, 3.476, 3.520, 0.687),
                  var = c(9.028, 6.484, 6.125, 5.187))
mod <- rma(tau_hat, var, data=dat, method="REML")
forest(mod, addpred = T)

#lm
dat <- data.frame(tau_hat = c(-0.218, 2.823, 2.470, 0.649),
                  var = c(1.723, 1.507, 1.223, 1.707))
mod <- rma(tau_hat, var, data=dat, method="REML")
forest(mod, addpred = T)

#lm - id3
dat <- data.frame(tau_hat = c(1.777, 3.411, 4.785, 0.914),
                  var = c(1.554, 1.358, 1.601, 2.033))
mod <- rma(tau_hat, var, data=dat, method="REML")
forest(mod, addpred = T)
