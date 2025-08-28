## Plotting PI vs CIs ##

library(tidyverse)
library(lme4)
library(rsample)
library(multcomp)
library(MASS)
library(grf)
library(dbarts)
library(fastDummies)

source("R/MDD_Generation.R")
source("R/MA_Stage1.R")
source("R/CF_Stage1.R")
source("R/BART_Stage1.R")
source("R/2Stage_MA.R")



### NEW
K <- 10
n_mean <- 500
n_sd <- 0
covars_fix <- "age"
covars_rand <- "age"
eps_study_m <- 1
eps_study_tau <- 1
eps_study_inter <- 0.5
lin <- T
distribution <- "variable"
seed <- 100
target_sample_hetintercept <- readRDS("Data/target_sample_hetintercept.RDS")
target_sample_hetinteraction <- readRDS("Data/target_sample_hetinteraction.RDS")

#low het
target_sample <- target_sample_hetinteraction
set.seed(seed)
sim_dat <- gen_mdd(K, n_mean, n_sd, covars_fix, covars_rand, lin,
                                               eps_study_m, eps_study_tau, eps_study_inter, 
                                               distribution, target_sample)
#training
train_dat <- sim_dat[["train_dat"]]
study_dat <- train_dat %>% group_split(S)
#target
target_dat <- sim_dat[["target_dat"]]

## Causal forest
target_dat_study_cf <- cf_stage1(study_dat, target_dat, K, covars)
#View(target_dat_study_cf)
target_dat_study_cf <- target_dat_study_cf %>% 
  mutate(lower = tau_hat - 1.96*sqrt(var), 
         upper = tau_hat + 1.96*sqrt(var),
         S = as.character(S))
target_dat_cf <- ma_stage2(target_dat_study_cf, K, target_dat, stat="t")
#View(target_dat_cf)

target_dat_cf <- target_dat_cf %>% 
  rename(lower = pi.lb, 
         upper = pi.ub,
         tau_hat = pred) %>%
  mutate(S = "Target")

target_dat_cf_all <- bind_rows(target_dat_study_cf, target_dat_cf)
rand_id <-sample(unique(target_dat_cf_all$id), 1)

target_dat_cf_all %>%
  filter(id == rand_id) %>%
  mutate(Study = paste("Study", S),
         Study = ifelse(Study == "Study Target", "Target", Study),
         color = ifelse(Study == "Target", 1, 0)) %>%
  ggplot(aes(x = Study, y = tau_hat, color = color)) +
  geom_point(show.legend = F,
             position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), 
                position=position_dodge(width=0.5),
                show.legend = F) +
  geom_hline(aes(yintercept = tau), linetype = "dashed", color="grey") +
  coord_flip()



### OLD

## Set up parameters
K <- 10
n_mean <- 500
n_sd <- 0
covars_fix <- "age"
covars_rand <- "age"
eps_study_m <- 1
eps_study_tau <- 0.25
eps_study_inter <- 0.25
lin <- T
distribution <- "variable"
seed <- 100


## Main function
plot_intervals <- function(K=10, n_mean=500, n_sd=0,
                           covars_fix="age", covars_rand="age",
                           eps_study_m=1, eps_study_tau=0.25,
                           eps_study_inter=0.25, lin=T,
                           distribution="same", seed=100) {
  
  target_sample_hetintercept <- readRDS("Data/target_sample_hetintercept.RDS")
  target_sample_hetinteraction <- readRDS("Data/target_sample_hetinteraction.RDS")
  
  
  ## Simulate training and target (OOS) data
  set.seed(seed)
  sim_dat <- gen_mdd(K, n_mean, n_sd, covars_fix, covars_rand, lin,
                     eps_study_m, eps_study_tau, eps_study_inter, 
                     distribution, target_sample_hetintercept)
  train_dat <- sim_dat[["train_dat"]]
  target_dat <- sim_dat[["target_dat"]]
  
  
  ## Fit causal forest
  covars <- c("sex", "smstat", "weight", "age", "madrs")
  
  feat <- dplyr::select(train_dat, c(S, all_of(covars))) %>%
    fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)
  
  tau_forest <- grf::causal_forest(X=feat, Y=train_dat$Y, W=train_dat$W, 
                                   num.threads=3, honesty=T, num.trees=1000)
  tau_hat <- predict(tau_forest, estimate.variance=T)
  
  
  ## Get all intervals in target
  #pi in target sample
  cf_target <- cf_pi_target(K, target_dat, tau_forest, covars) %>%
    dplyr::select(-c(var_within, var_btwn, var_tot, n_K)) %>%
    mutate(interval_type = "Target",
           interval_cat = "PI")
  cf_all <- cf_target
  
  #pretend like target sample came from each study
  target_feat <- target_dat %>%
    dplyr::select(sex, smstat, weight, age, madrs) %>%
    mutate(S_1=0, S_2=0, S_3=0, S_4=0, S_5=0, S_6=0, S_7=0, S_8=0, S_9=0, S_10=0)
  #predict according to each study
  for (i in 1:10) {
    #set features to be in study i
    target_feati <- target_feat
    target_feati[,i+5] <- 1
    
    #predict for study i
    cf_targeti <- predict(tau_forest, target_feati, estimate.variance = T)
    
    #create data frame with predictions
    target_dati <- target_dat %>%
      mutate(mean = cf_targeti$predictions,
             sd = sqrt(cf_targeti$variance.estimates),
             lower = mean - 1.96*sd,
             upper = mean + 1.96*sd,
             interval_type = paste0("Study ", i),
             interval_cat = "Pooled CI")
    
    cf_all <- bind_rows(cf_all, target_dati)
  }
  
  
  ## Add in study-specific cfs
  for (i in 1:10) {
    #subset to study i
    train_dat_i <- train_dat %>%
      filter(S == i)
    train_feat <- train_dat_i %>%
      dplyr::select(all_of(covars))
    
    #fit causal forest to study i
    tau_forest_study <- grf::causal_forest(X=train_feat, Y=train_dat_i$Y, W=train_dat_i$W, 
                                           num.threads=3, honesty=T, num.trees=1000)
    
    #predict on target setting
    target_feat_nostudy <- target_dat %>%
      dplyr::select(all_of(covars))
    cf_target_study <- predict(tau_forest_study, target_feat_nostudy, estimate.variance = T)
    
    #create data frame with predictions
    target_dat_studyi <- target_dat %>%
      mutate(mean = cf_target_study$predictions,
             sd = sqrt(cf_target_study$variance.estimates),
             lower = mean - 1.96*sd,
             upper = mean + 1.96*sd,
             interval_type = paste0("Study ", i),
             interval_cat = "Study-Specific CI")
    
    #add to main data frame
    cf_all <- bind_rows(cf_all, target_dat_studyi)
  }
  
  cf_all <- cf_all %>%
    mutate(setting = ifelse(eps_study_m == 0.05, "Low Heterogeneity",
                            ifelse(eps_study_tau == 0.05, "Heterogeneous Intercept",
                                   ifelse(eps_study_tau == 0.5, "Heterogeneous Main",
                                          "High Heterogeneity"))))
  
  return(cf_all)
}


## Save one for each setting
cf_all_lowhet <- plot_intervals(eps_study_m=0.05, eps_study_tau=0.05,
                                eps_study_inter=0.05)
cf_all_hetint <- plot_intervals(eps_study_m=1, eps_study_tau=0.05,
                                eps_study_inter=0.05)
cf_all_hetmain <- plot_intervals(eps_study_m=1, eps_study_tau=0.5,
                                 eps_study_inter=0.05)
cf_all_highhet <- plot_intervals(eps_study_m=1, eps_study_tau=1,
                                 eps_study_inter=0.5)

#combine all
cf_all_settings <- bind_rows(cf_all_lowhet, cf_all_hetint,
                             cf_all_hetmain, cf_all_highhet)

## Plot intervals
#just one id, low het
# cf_all_lowhet %>%
#   filter(id == 1) %>%
#   mutate(interval_type = factor(interval_type,
#                                 levels = c("S1","S2","S3","S4","S5","S6",
#                                            "S7","S8","S9","S10","Target")),
#          interval_cat = factor(interval_cat,
#                                levels = c("Study-Specific CI", "Pooled CI", "PI"))) %>%
#   ggplot(aes(x=interval_type, y=mean, color=interval_cat)) +
#   geom_point(position=position_dodge(width=0.5)) +
#   coord_flip() +
#   geom_errorbar(aes(ymin=lower, ymax=upper), position=position_dodge(width=0.5)) +
#   geom_hline(aes(yintercept = tau), linetype = "dashed", color="grey") +
#   xlab("CATE Prediction Setting") +
#   ylab("CATE for Covariate Profile (95% Interval)") +
#   labs(color="Interval Type") +
#   scale_color_manual(values=c("#00BFC4","#F8766D","black"))
cf_all_lowhet %>%
  mutate(interval_type = factor(interval_type,
                                levels = c("Study 1","Study 2","Study 3",
                                           "Study 4","Study 5","Study 6",
                                           "Study 7","Study 8","Study 9",
                                           "Study 10","Target")),
         interval_cat = factor(interval_cat,
                               levels = c("Study-Specific CI", "Pooled CI", "PI"))) %>%
  filter(id == 1,
         interval_cat != "Study-Specific CI") %>%
  ggplot(aes(x=interval_type, y=mean, color=interval_cat)) +
  geom_point(position=position_dodge(width=0.5), 
             show.legend = F) +
  coord_flip() +
  geom_errorbar(aes(ymin=lower, ymax=upper), 
                position=position_dodge(width=0.5),
                show.legend = F) +
  geom_hline(aes(yintercept = tau), linetype = "dashed", color="grey") +
  xlab("Estimation/Prediction Setting") +
  ylab("CATE for Covariate Profile (95% Interval)") +
  scale_color_manual(values=c("#00BFC4","black")) +
  theme(text = element_text(size = 13))
ggsave(paste("Results/EHR/allintervals_oneid", seed,
             0.05, 0.05, 0.05, lin,
             "24May2024.jpeg", sep="_"),
       width = 7, height = 5)

cf_all_hetmain %>%
  mutate(interval_type = factor(interval_type,
                                levels = c("Study 1","Study 2","Study 3",
                                           "Study 4","Study 5","Study 6",
                                           "Study 7","Study 8","Study 9",
                                           "Study 10","Target")),
         interval_cat = factor(interval_cat,
                               levels = c("Study-Specific CI", "Pooled CI", "PI"))) %>%
  filter(id == 1,
         interval_cat != "Study-Specific CI") %>%
  ggplot(aes(x=interval_type, y=mean, color=interval_cat)) +
  geom_point(position=position_dodge(width=0.5), 
             show.legend = F) +
  coord_flip() +
  geom_errorbar(aes(ymin=lower, ymax=upper), 
                position=position_dodge(width=0.5),
                show.legend = F) +
  geom_hline(aes(yintercept = tau), linetype = "dashed", color="grey") +
  xlab("Estimation/Prediction Setting") +
  ylab("CATE for Covariate Profile (95% Interval)") +
  scale_color_manual(values=c("#00BFC4","black")) +
  theme(text = element_text(size = 13))
ggsave(paste("Results/EHR/allintervals_oneid", seed,
             1, 0.5, 0.05, lin,
             "24May2024.jpeg", sep="_"),
       width = 7, height = 5)

#all settings
cf_all_settings %>%
  filter(id == 1) %>%
  mutate(interval_type = factor(interval_type,
                                levels = c("S1","S2","S3","S4","S5","S6",
                                           "S7","S8","S9","S10","Target")),
         interval_cat = factor(interval_cat,
                               levels = c("Study-Specific CI", "Pooled CI", "PI"))) %>%
  ggplot(aes(x=interval_type, y=mean, color=interval_cat)) +
  geom_point(position=position_dodge(width=0.5)) +
  coord_flip() +
  geom_errorbar(aes(ymin=lower, ymax=upper), position=position_dodge(width=0.5)) +
  geom_hline(aes(yintercept = tau), linetype = "dashed", color="grey") +
  facet_wrap(~setting) +
  xlab("Model Setting") +
  ylab("CATE for Covariate Profile (95% Interval)") +
  labs(color="Interval Type") +
  scale_color_manual(values=c("#00BFC4","#F8766D","black"))

#two settings
cf_all_settings %>%
  filter(id == 1,
         setting %in% c("Low Heterogeneity", "Heterogeneous Main")) %>%
  mutate(interval_type = factor(interval_type,
                                levels = c("S1","S2","S3","S4","S5","S6",
                                           "S7","S8","S9","S10","Target")),
         interval_cat = factor(interval_cat,
                               levels = c("Study-Specific CI", "Pooled CI", "PI")),
         setting = factor(setting,
                          levels = c("Low Heterogeneity", "Heterogeneous Main"))) %>%
  ggplot(aes(x=interval_type, y=mean, color=interval_cat)) +
  geom_point(position=position_dodge(width=0.5)) +
  coord_flip() +
  geom_errorbar(aes(ymin=lower, ymax=upper), position=position_dodge(width=0.5)) +
  geom_hline(aes(yintercept = tau), linetype = "dashed", color="grey") +
  facet_wrap(~setting) +
  xlab("Model Setting") +
  ylab("CATE for Covariate Profile (95% Interval)") +
  labs(color="Interval Type") +
  scale_color_manual(values=c("#00BFC4","#F8766D","black")) +
  theme(text = element_text(size = 13))
ggsave(paste("Results/EHR/allintervals_twosettings_oneid", seed,
             "24May2024.jpeg", sep="_"),
       width = 8, height = 5)


## Calculate ratio of lengths
#pooled / study-specific ci
cf_all_settings %>%
  mutate(interval_length = upper - lower) %>%
  filter(interval_type != "Target") %>%
  dplyr::select(id, interval_length, interval_type, interval_cat, setting) %>%
  pivot_wider(names_from = c(interval_cat),
              values_from = c(interval_length)) %>%
  mutate(ratio = `Pooled CI`/`Study-Specific CI`) %>% 
  group_by(setting) %>% 
  summarise(pooled_over_studyspecific = mean(ratio))

#pooled ci / target pi
#first get avg pooled ci
id_setting_avg <- cf_all_settings %>%
  mutate(interval_length = upper - lower) %>%
  filter(interval_cat == "Pooled CI") %>%
  group_by(id, setting) %>%
  summarise(avg_pooled = mean(interval_length))

#combine with target
cf_all_settings %>%
  filter(interval_type == "Target") %>%
  left_join(id_setting_avg, by=c("id","setting")) %>%
  mutate(ratio = avg_pooled/(upper - lower)) %>%
  group_by(setting) %>%
  summarise(pooled_over_pi = mean(ratio),
            pi_over_pooled = mean((upper-lower)/avg_pooled))



