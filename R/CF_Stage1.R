### MA Stage 1: Causal Forest ###

#get estimates for all studies
cf_stage1 <- function(study_dat, target_dat, K, covars, honesty = T) {
  
  target_dat_pred <- data.frame()
  for (i in 1:length(study_dat)) {
    
    #filter to study i
    train_df <- study_dat[[i]]
    
    #fit forest
    feat <- dplyr::select(train_df, all_of(covars))
    tau_forest <- grf::causal_forest(X=feat, Y=train_df$Y, W=train_df$W, 
                                     num.threads=3, honesty=honesty, num.trees=1000)
    
    #predict on target
    target_feat <- dplyr::select(target_dat, all_of(covars))
    target_pred <- predict(tau_forest, newdata = target_feat, estimate.variance = T)
    
    target_df <- target_dat %>%
      mutate(S = i,
             tau_hat = target_pred$predictions,
             var = target_pred$variance.estimates)
    target_dat_pred <- bind_rows(target_dat_pred,
                                 target_df)
  }
  
  return(target_dat_pred)
}
