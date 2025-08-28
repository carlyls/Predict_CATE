### MA Stage 1: Linear Model ###

#calculate CATEs and variances
lm_pred <- function(lin_mod, df, K) {
  
  #calculate variances
  fc <- vcov(lin_mod) #covariance matrix of fixed effects
  
  #fixed effects
  beta <- coef(lin_mod)[grep("W", names(coef(lin_mod)))] %>% matrix() #beta-hat
  var_beta <- fc[grep("W", rownames(fc)),
                 grep("W", colnames(fc))] #Var(beta-hat)
  
  #calculate theta-hats
  X <- df %>%
    dplyr::select(W, "age") %>%
    mutate(W = 1) %>%
    as.matrix()
  
  mean_theta <- X %*% beta %>% c()
  vcov_theta <- X %*% var_beta %*% t(X)
  var_theta <- diag(vcov_theta) #Var(theta-hat)
  
  df <- df %>%
    mutate(tau_hat = mean_theta,
           var = var_theta)
  
  return(df)
}

#get estimates for all studies
lm_stage1 <- function(study_dat, target_dat, K, covars) {
  
  #set formula
  formula <- as.formula(paste0("Y ~ W + ", paste(covars, collapse = " + "),
                               " + W:age"))
  
  #iterate through all studies to predict on target data
  target_dat_pred <- data.frame()
  for (i in 1:length(study_dat)) {
    lin_mod <- lm(formula, data=study_dat[[i]])
    df <- target_dat %>%
      mutate(S = i)
    target_dat_pred <- bind_rows(target_dat_pred,
                                 lm_pred(lin_mod, df, K))

  }
  
  return(target_dat_pred)
}
