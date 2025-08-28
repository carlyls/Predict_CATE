### Running Simulations for Methods to Estimate on OOS Individuals ###

library(tidyverse)
library(lme4)
library(rsample)
library(multcomp)
library(MASS)
library(grf)
library(dbarts)
library(fastDummies)
library(here)

source("R/MDD_Generation.R")
source("R/MA_Stage1.R")
source("R/CF_Stage1.R")
source("R/BART_Stage1.R")
source("R/2Stage_MA.R")

# set up parameters
n_mean <- 500
n_sd <- 0
distribution <- "variable"

#sims
mods <- list(list(covars_fix="age", covars_rand="age", lin=T,
                  eps_study_m=1, eps_study_tau=0.25, eps_study_inter=0.25,
                  target_lab = "hetintercept"),
             list(covars_fix="age", covars_rand="age", lin=T,
                  eps_study_m=1, eps_study_tau=0.5, eps_study_inter=0.25,
                  target_lab = "hetmain"),
             list(covars_fix="age", covars_rand="age", lin=T,
                  eps_study_m=1, eps_study_tau=1, eps_study_inter=0.5,
                  target_lab = "hetinteraction"),
             list(covars_fix="age", covars_rand="age", lin=F,
                  eps_study_m=1, eps_study_tau=0.25, eps_study_inter=0.25,
                  target_lab = "hetintercept"),
             list(covars_fix="age", covars_rand="age", lin=F,
                  eps_study_m=1, eps_study_tau=0.5, eps_study_inter=0.25,
                  target_lab = "hetmain"),
             list(covars_fix="age", covars_rand="age", lin=F,
                  eps_study_m=1, eps_study_tau=1, eps_study_inter=0.5,
                  target_lab = "hetinteraction"))

settings <- expand.grid(iteration = c(1:500),
                        K = c(4, 30),
                        moderators = c(1:length(mods)))

#set row of settings and define parameters
#i=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

#run main function
for (i in 1:nrow(settings)) {
  
  K <- settings$K[i]
  moderators <- settings$moderators[i]
  covars_fix <- mods[[moderators]]$covars_fix
  covars_rand <- mods[[moderators]]$covars_rand
  eps_study_m <- mods[[moderators]]$eps_study_m
  eps_study_tau <- mods[[moderators]]$eps_study_tau
  eps_study_inter <- mods[[moderators]]$eps_study_inter
  lin <- mods[[moderators]]$lin
  target_lab <- mods[[moderators]]$target_lab
  
  iteration <- settings$iteration[i]
  target_sample <- readRDS(here("Data", paste0("target_sample_", target_lab, ".RDS")))
  seed <- i
  
  set.seed(seed)
  results <- try({compare_2stage(K=K, n_mean=n_mean, n_sd=n_sd, covars_fix=covars_fix,
                                 covars_rand=covars_rand, lin=lin, eps_study_m=eps_study_m, 
                                 eps_study_tau=eps_study_tau, eps_study_inter=eps_study_inter, 
                                 distribution=distribution, target_sample=target_sample)},
                 silent = T)
  if (inherits(results, "try-error")) {
    print(paste("Converging error on iteration", i))
  } else{
    save(results, file=here("Data", "Different K", paste(paste("results",seed,iteration,K,n_mean,n_sd,target_lab,moderators,
                                                             lin,eps_study_m,eps_study_tau,mean(eps_study_inter),distribution,sep = "_"),
                                                       ".Rdata",sep="")))
  }
  
}

