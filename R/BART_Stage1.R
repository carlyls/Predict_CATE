## Fitting BART to Simulated Data


#get estimates for all studies
sbart_stage1 <- function(study_dat, target_dat, K, covars) {
  
  target_dat_pred <- data.frame()
  for (i in 1:length(study_dat)) {
    
    #filter to study i
    train_df <- study_dat[[i]]
    
    #fit forest
    feat <- dplyr::select(train_df, c(W, all_of(covars)))
    feat_cf <- feat %>%
      mutate(W = as.numeric(W == 0))
    y <- as.numeric(train_df$Y)
    
    sbart <- dbarts::bart(x.train=as.matrix(feat), 
                          y.train=y, 
                          x.test=as.matrix(feat_cf), 
                          keeptrees=T)
    
    #predict on target
    target_feat <- dplyr::select(target_dat, c(W, all_of(covars)))
    target_feat_cf <- target_feat %>%
      mutate(W = as.numeric(W == 0))

    target_pred <- predict(sbart, target_feat)
    target_pred_cf <- predict(sbart, target_feat_cf)
    
    #get mean and variances of predictions
    mean_obs <- apply(target_pred, 2, mean)
    mean_cf <- apply(target_pred_cf, 2, mean)
    var_obs <- apply(target_pred, 2, var)
    var_cf <- apply(target_pred_cf, 2, var)
    w_fac_new <- ifelse(target_dat$W == 0, -1, 1)
    mean_cate <- w_fac_new*(mean_obs - mean_cf)
    var_cate <- var_obs + var_cf
    
    #add to target
    target_df <- target_dat %>%
      mutate(S = i,
             tau_hat = mean_cate,
             var = var_cate)
    target_dat_pred <- bind_rows(target_dat_pred,
                                 target_df)
  }
  
  return(target_dat_pred)
}
