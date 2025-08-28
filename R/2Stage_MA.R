### Trying Methods for OOS Estimation in Random Effects Meta-Analysis and Causal Forests ###

library(tidyverse)
library(lme4)
library(rsample)
library(multcomp)
library(MASS)
library(grf)
library(nnet)
library(metafor)


#check results for all methods ####
assess_interval <- function(target_dat) {
  
  ## within-iteration results
  #calculate mean absolute bias
  target_bias <- mean(abs(target_dat$pred - target_dat$tau))
  
  #calculate mse
  target_mse <- mean((target_dat$pred - target_dat$tau)^2)
  
  #calculate coverage
  target_coverage <- sum(target_dat$tau >= target_dat$pi.lb & 
                           target_dat$tau <= target_dat$pi.ub)/nrow(target_dat)
  
  #calculate length
  target_length <- mean(target_dat$pi.ub - target_dat$pi.lb)
  
  #calculate signficance
  target_significance <- sum(sign(target_dat$pi.lb) == sign(target_dat$pi.ub))/nrow(target_dat)
  
  return(c(target_bias = target_bias,
           target_mse = target_mse,
           target_coverage = target_coverage,
           target_length = target_length,
           target_significance = target_significance))
  
}

#add diagnostics for each row of target sample according to method
target_metrics <- function(target_update, method) {
  
  df <- target_update %>%
    mutate(method = method,
           coverage = as.numeric(tau >= pi.lb & tau <= pi.ub),
           bias = pred - tau,
           abs_bias = abs(pred - tau),
           length = pi.ub - pi.lb,
           significant = as.numeric(sign(pi.lb) == sign(pi.ub)))
  
  return(df)
}

#form prediction intervals using 2-stage meta-analysis (REML, Riley methods)
form_pi_t <- function(df) {
  
  #fit 2-stage model
  mod_stage2 <- rma(tau_hat, var, data = df, method = "REML")
  
  #form prediction intervals
  res <- predict(mod_stage2, pi.type = "Riley") %>% data.frame()
  
  return(res)
  
}

#form prediction intervals using 2-stage meta-analysis (REML, Riley methods)
form_pi_z <- function(df) {
  
  #fit 2-stage model
  mod_stage2 <- rma(tau_hat, var, data = df, method = "REML")
  
  #form prediction intervals
  res <- predict(mod_stage2) %>% data.frame()
  
  return(res)
  
}

#apply 2-stage meta-analysis and return target data
ma_stage2 <- function(df, K, target_dat, stat="t") {
  
  #split by id and then form pis within each id
  if (stat == "t") {
    
    pis <- df %>%
      group_split(id) %>%
      map(form_pi_t) %>%
      bind_rows()  
    
  } else if (stat == "z") {
    
    pis <- df %>%
      group_split(id) %>%
      map(form_pi_z) %>%
      bind_rows()
    
  }
  
  #combine with target data
  target_dat_final <- target_dat %>%
    bind_cols(pis)
  
  return(target_dat_final)
    
}


#overall function ####
compare_2stage <- function(K=10, n_mean=500, n_sd=0, covars_fix="age", covars_rand="age",
                        covars=c("age","weight","madrs","sex","smstat"),
                        lin=T, eps_study_m=0.05, eps_study_tau=0.05, eps_study_inter=0.05, 
                        distribution="same", target_sample) {
  
  
  ## Simulate training and target (OOS) data
  sim_dat <- gen_mdd(K, n_mean, n_sd, covars_fix, covars_rand, lin,
                     eps_study_m, eps_study_tau, eps_study_inter, 
                     distribution, target_sample)
  
  #training
  train_dat <- sim_dat[["train_dat"]]
  study_dat <- train_dat %>%
    group_split(S)
  
  #target
  target_dat <- sim_dat[["target_dat"]]
  
  
  ## Linear model
  target_dat_study_lm <- lm_stage1(study_dat, target_dat, K, covars)
  target_dat_lm <- ma_stage2(target_dat_study_lm, K, target_dat, stat="t")
  lm_res <- assess_interval(target_dat_lm)
  
  ## Linear model - Z-statistic
  # target_dat_lmz <- ma_stage2(target_dat_study_lm, K, target_dat, stat="z")
  # lmz_res <- assess_interval(target_dat_lmz)
  
  
  ## Causal forest
  target_dat_study_cf <- cf_stage1(study_dat, target_dat, K, covars)
  target_dat_cf <- ma_stage2(target_dat_study_cf, K, target_dat, stat="t")
  cf_res <- assess_interval(target_dat_cf)
  
  ## Causal forest - Z-statistic
  # target_dat_cfz <- ma_stage2(target_dat_study_cf, K, target_dat, stat="z")
  # cfz_res <- assess_interval(target_dat_cfz)
  
  ## Causal forest - adaptive
  target_dat_study_cfa <- cf_stage1(study_dat, target_dat, K, covars, honesty = F)
  target_dat_cfa <- ma_stage2(target_dat_study_cfa, K, target_dat, stat="t")
  cfa_res <- assess_interval(target_dat_cfa)
  
  
  ## BART
  target_dat_study_sbart <- sbart_stage1(study_dat, target_dat, K, covars)
  target_dat_sbart <- ma_stage2(target_dat_study_sbart, K, target_dat, stat="t")
  sbart_res <- assess_interval(target_dat_sbart)
  
  ## BART - Z-statistic
  # target_dat_sbartz <- ma_stage2(target_dat_study_sbart, K, target_dat, stat="z")
  # sbartz_res <- assess_interval(target_dat_sbartz)
  
  
  ## Save results
  #data frame of parameters
  params <- data.frame(K=K, n_mean=n_mean, n_sd=n_sd,
                       covars_fix=paste(covars_fix, collapse = ", "),
                       covars_rand=paste(covars_rand, collapse = ", "), lin=lin,
                       eps_study_m=eps_study_m, eps_study_tau=eps_study_tau, 
                       eps_study_inter=paste(eps_study_inter, collapse = ", "),
                       distribution=distribution)
  
  #data frame of results
  all_res <- cbind(lm_res, cf_res, cfa_res, sbart_res) %>%
    data.frame() %>%
    rownames_to_column("Metric") %>%
    cbind(params)
  
  #individual prediction errors from each method
  ipe <- target_dat %>%
    mutate(lm_ipe = target_dat_lm$tau - target_dat_lm$pred,
           cf_ipe = target_dat_cf$tau - target_dat_cf$pred,
           cfa_ipe = target_dat_cfa$tau - target_dat_cfa$pred,
           sbart_ipe = target_dat_sbart$tau - target_dat_sbart$pred)
  
  #save target dataset individual results
  target_res <- target_metrics(target_dat_lm, "LM") %>%
    bind_rows(target_metrics(target_dat_cf, "HCF")) %>%
    bind_rows(target_metrics(target_dat_cfa, "ACF")) %>%
    bind_rows(target_metrics(target_dat_sbart, "SBART")) %>%
    cbind(params)
  
  
  return(list(all_res=all_res, ipe=ipe, target_res=target_res))
}

