
VMP_fit_objs <- function(sim_dataset, N, domain){
  list_y = map(1:N, function(x){
    return(sim_dataset$Y[x,])
  })
  list_t = rep(list(domain), N)
  return(list(Yl = list_y, Tl = list_t))
}

VMP_extract <- function(vmp_mod, domain){
  EF_est = FPC_df(vmp_mod$list_Psi_hat, domain)
  
  Score_est = data.frame(t(vmp_mod$Zeta_hat))
  colnames(Score_est) = paste0("Curve ", 1:ncol(Score_est))
  Score_est$FPC = paste0("FPC ", 1:nrow(Score_est))
  Score_est = Score_est %>%
    pivot_longer(-c(FPC), names_to = "Curve", values_to = "Est")
  
  # Calculate CI using normal approximation
  N = length(vmp_mod$Cov_zeta_hat)
  Score_SD = map(1:N, function(x){
    sds = sqrt(diag(vmp_mod$Cov_zeta_hat[[x]]))
    return(data.frame(Curve = paste0("Curve ", x), 
                      FPC = paste0("FPC ", 1:length(sds)), 
                      SD = sds))
  }) %>% list_rbind()
  Score_est = left_join(Score_est, Score_SD, by = c("Curve", "FPC")) %>%
    mutate(LB = Est - 1.96*SD, UB = Est + 1.96*SD) %>%
    select(-c(SD))
  
  Mu_est = data.frame(Arg = domain, Est = vmp_mod$mu_hat)
  
  Smooth_est = map(1:N, function(x){
    output = data.frame(Est = vmp_mod$Y_hat[[x]], 
                        LB = vmp_mod$Y_low[[x]], 
                        UB = vmp_mod$Y_upp[[x]], 
                        Arg = domain, 
                        Curve = paste0("Curve ", x))
    return(output)
  }) %>% list_rbind()
    
  return(list(EF = EF_est, Score = Score_est, Mu = Mu_est, Smooth = Smooth_est))
}
