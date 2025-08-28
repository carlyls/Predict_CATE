## Confirm Prediction Interval Code

library(tidyverse)
library(multicate)
library(metafor)

#read in simulated data
setwd("/Users/cb644/Library/CloudStorage/Box-Box/Home Folder cb644/Documents/Moderation/R Package")
dat <- readRDS("sample_dat_estimateCATE.RDS")
K <- 10

#estimate cate
feat <- dat %>%
  dplyr::select(c(S, X1, X2, X3, X4, X5)) %>%
  fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)
tau_forest <- grf::causal_forest(X = feat,
                                 W = dat$W,
                                 Y = dat$Y)

# cate <- estimate_cate(trial_tbl = dat,
#                       estimation_method = "causalforest",
#                       aggregation_method = "studyindicator",
#                       study_col = "S",
#                       treatment_col = "W",
#                       outcome_col = "Y",
#                       covariate_col = c("X1", "X2", "X3", "X4", "X5"))
# res <- cate$model

#simulate target data
target_dat <- data.frame(X1 = rnorm(100),
                         X2 = rnorm(100),
                         X3 = rnorm(100),
                         X4 = rnorm(100),
                         X5 = rnorm(100))

#predict cate - manual
#assign study
new_dat <- target_dat %>%
  slice(rep(1:n(), each=K)) %>%
  mutate(S = rep(1:K, nrow(target_dat)))
new_feat <- new_dat %>%
  dplyr::select(c(S, X1, X2, X3, X4, X5)) %>%
  fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)

#predict CATE
target_pred <- predict(tau_forest, newdata = new_feat, estimate.variance = T)
new_dat <- new_dat %>%
  mutate(tau_hat = target_pred$predictions,
         var = target_pred$variance.estimates)

#focus on profile 1
new_dat_1 <- new_dat[1:10,]

#try metafor package
mod <- rma(tau_hat, var, data = new_dat_1, method = "DL")
summary(mod)
predict(mod)
predict(mod, pi.type = "Riley")

#compare with my approach - raw
avg <- mean(new_dat_1$tau_hat)
within_var <- mean(new_dat_1$var)
between_var <- var(new_dat_1$tau_hat)
(lower_raw <- avg - qt(.975, K-2)*sqrt(within_var + between_var))
(upper_raw <- avg + qt(.975, K-2)*sqrt(within_var + between_var))

#compare with my approach - inv var weight
wts <- 1/(new_dat_1$var + between_var)
wt_avg <- weighted.mean(new_dat_1$tau_hat, wts)
(lower_wt <- wt_avg - qt(.975, K-2)*sqrt(within_var + between_var))
(upper_wt <- wt_avg + qt(.975, K-2)*sqrt(within_var + between_var))

#compare with manual approach
# Extract treatment effects and variances
theta_i <- data[[theta_col]]
var_i <- data[[var_col]]

# Confidence interval z-value
z_crit <- qnorm((1 + conf_level) / 2)

# Fixed-Effects Meta-Analysis
theta_i <- new_dat_1$tau_hat
var_i <- new_dat_1$var
w_i <- 1 / var_i
theta_FE <- sum(w_i * theta_i) / sum(w_i)
var_FE <- 1 / sum(w_i)
se_FE <- sqrt(var_FE)

# Random-Effects Meta-Analysis
# Step 1: Cochran's Q
theta_FE_all <- theta_FE
Q <- sum(w_i * (theta_i - theta_FE_all)^2)

# Step 2: Tau^2 (heterogeneity)
K <- length(theta_i)
tau2 <- max(0, (Q - (K - 1)) / (sum(w_i) - (sum(w_i^2) / sum(w_i))))

# Step 3: Random-effects weights
w_star <- 1 / (var_i + tau2)

# Step 4: Random-effects pooled estimate
theta_RE <- sum(w_star * theta_i) / sum(w_star)
var_RE <- 1 / sum(w_star)
se_RE <- sqrt(var_RE)
ci_RE <- c(
  theta_RE - z_crit * se_RE,
  theta_RE + z_crit * se_RE
)

# Prediction Interval (Random-Effects)
var_new <- tau2 + var_RE
PI <- c(
  theta_RE - qt(0.975, K-2) * sqrt(var_new),
  theta_RE + qt(0.975, K-2) * sqrt(var_new)
)


#for each X profile, calculate avg tau_hat, avg within-study var, variance of tau_hats
cis <- new_dat %>%
  group_by(W, sex, smstat, weight, age, age2, madrs, Y, tau) %>%
  summarise(mean = mean(tau_hat),
            var_within = mean(var),
            var_btwn = var(tau_hat),
            n_K = n()) %>%
  mutate(var_tot = var_within + var_btwn,
         sd = sqrt(var_tot),
         lower = mean - qt(.975, K-2)*sd,
         upper = mean + qt(.975, K-2)*sd)

#reorder to match target data
cis_ord <- target_dat %>%
  left_join(cis, by = c("W","sex","smstat","weight","age",
                        "age2","madrs","Y","tau"))


