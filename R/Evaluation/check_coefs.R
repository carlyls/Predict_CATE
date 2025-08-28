## Check coefficient deviations for simulations

library(tidyverse)

#read in data
target_sample_hetintercept <- readRDS("~/Library/CloudStorage/Box-Box/Home Folder cb644/Documents/Research/Carly Projects/Moderation/Predict_CATE/Data/target_sample_hetintercept.RDS")
target_sample_hetmain <- readRDS("~/Library/CloudStorage/Box-Box/Home Folder cb644/Documents/Research/Carly Projects/Moderation/Predict_CATE/Data/target_sample_hetmain.RDS")
target_sample_hetinteraction <- readRDS("~/Library/CloudStorage/Box-Box/Home Folder cb644/Documents/Research/Carly Projects/Moderation/Predict_CATE/Data/target_sample_hetinteraction.RDS")


#compare coefficients on their distributions
truths <- data.frame(
  intercept_tau = rnorm(10000, mean=0, sd=0.25),
  intercept_age = rnorm(10000, mean=0, sd=0.25),
  main_tau = rnorm(10000, mean=0, sd=0.5),
  main_age = rnorm(10000, mean=0, sd=0.25),
  interaction_tau = rnorm(10000, mean=0, sd=1),
  interaction_age = rnorm(10000, mean=0, sd=0.5)
)

ggplot(truths, aes(x = intercept_tau)) + 
  geom_density()


#calculate z scores
unique(target_sample_hetintercept$eps_tau)/0.25
unique(target_sample_hetintercept$eps_age)/0.25

unique(target_sample_hetmain$eps_tau)/0.5
unique(target_sample_hetmain$eps_age)/0.25

unique(target_sample_hetinteraction$eps_tau)/1
unique(target_sample_hetinteraction$eps_age)/0.5
