# Compare spanned space
phi_span <- function(trueM, estM, M_val){
  for_rb = map(estM, function(est){return(est/sqrt(M_val))})
  samples_riem = riemfactory(for_rb, name = "grassmann")
  true_riem = riemfactory(list(trueM/sqrt(M_val)), name = "grassmann")
  dists = rbase.pdist2(true_riem, samples_riem)[1,]
  return(mean(dists))
}

# Summary of FE performance
FE_inf <- function(trueFE_df, sampleFE_df, method){
  jdf = inner_join(trueFE_df, sampleFE_df, by = c("Arg")) %>%
    mutate(Cov = LB <= Func & UB >= Func, 
           SE = (Func - Est)^2) %>%
    summarize(Cov = mean(Cov), MSE = sum(wei * SE))
  jdf$Method = method
  return(jdf)
}

# Summary of EF performance
EF_inf <- function(trueEF_df, sampleEF_df, method){
  jdf = inner_join(trueEF_df, sampleEF_df, by = c("Arg", "FPC_Num")) %>%
    mutate(Cov = LB <= Func & UB >= Func, 
           SE = (Func - FPC_Val)^2) %>%
    group_by(FPC_Num) %>%
    summarize(Cov = mean(Cov), MSE = sum(wei * SE))
  jdf$Method = method
  return(jdf)
}

# Summary of score performance
Score_inf <- function(trueScore, sampleScore, method){
  jdf = inner_join(sampleScore, trueScore, by = c("FPC", "Curve")) %>%
    mutate(Cov = Score >= LB & Score <= UB, 
           SE = (Est - Score)^2) %>%
    group_by(FPC) %>%
    summarize(Cov = mean(Cov), MSE = mean(SE))
  jdf$Method = method
  return(jdf)
}
