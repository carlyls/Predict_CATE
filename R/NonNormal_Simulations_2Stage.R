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

source("R/MDD_Generation_NonNormal.R")
source("R/MA_Stage1.R")
source("R/CF_Stage1.R")
source("R/BART_Stage1.R")
source("R/2Stage_MA.R")

# set up parameters
K <- 10
n_mean <- 500
n_sd <- 0
covars_fix <- "age"
covars_rand <- "age"
distribution <- "variable"
target_sample <- readRDS(here("Data", "target_sample_nonnormal.RDS"))

settings <- expand.grid(lin = c(T, F),
                        iteration = c(1:500))

#set row of settings and define parameters
#i=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

#run main function
for (i in 1:nrow(settings)) {
  
  lin <- settings$lin[i]
  iteration <- settings$iteration[i]
  seed <- i
  
  set.seed(seed)
  results <- try({compare_2stage(K=K, n_mean=n_mean, n_sd=n_sd, covars_fix=covars_fix,
                            covars_rand=covars_rand, lin=lin,
                            distribution=distribution, target_sample=target_sample)},
                 silent = T)
  if (inherits(results, "try-error")) {
    print(paste("Converging error on iteration", i))
  } else{
    save(results, file=here("Data", "NonNormal", paste(paste("results",seed,iteration,K,n_mean,n_sd,"nonnormal",
                                   lin,distribution,sep = "_"),
                             ".Rdata",sep="")))
  }
  
}

