# General helper functions
source("../Gen_Funcs/Libs.R")
source("../Gen_Funcs/Bases.R")
source("../Gen_Funcs/Generate.R")
source("../Gen_Funcs/PostProcess.R")
source("../Gen_Funcs/Comparisons.R")

# Specific to each method
source("../FAST/FAST_Help.R")

where_directory = "Results"

create_dir <- function(direct){
  if (!file.exists(direct)){
    dir.create(direct, recursive = T)
  }
}

args = commandArgs(trailingOnly=TRUE)

out_dir = paste0(args[1], "_Alpha", args[4])

message("Simulation: ", args[1])
message("Output Directory: ", out_dir)
message("Output file: ", args[2])
message("Number of simulations: ", args[3])
message("Alpha value for penalty: ", args[4])

# Create output directories
{
  FE_dir = paste0(where_directory, "/FE/", out_dir)
  EF_dir = paste0(where_directory, "/EF/", out_dir)
  Score_dir = paste0(where_directory, "/Score/", out_dir)
  
  create_dir(FE_dir)
  create_dir(EF_dir)
  create_dir(Score_dir) 
}

# Begin simulations
n_sim = as.numeric(args[3])
N = 50
M = 50
Q = 20
K = 3
alpha = as.numeric(args[4])
weight_Obj = gaussLegendre(M, 0, 1)
sim_domain = weight_Obj$x

FComps = Funcs[[args[1]]]()
FPCs = FComps[-length(FComps)]
MuFun = FComps[[length(FComps)]]
EVals = Lambdas[[args[1]]]
sigma2 = Sigma2[[args[1]]]

V_out = data.frame(Method = c(), Spread = c(), Sample = c(), Alpha = c())
FE_out = data.frame(Method = c(), Mean_Cov = c(), MSE = c(), Sample = c(), Alpha = c())
EF_out = data.frame(Method = c(), EFNum = c(), Mean_Cov = c(), MSE = c(), Sample = c(), Alpha = c())
Score_out = data.frame(Method = c(), EFNum = c(), Mean_cov = c(), MSE = c(), Sample = c(), Alpha = c())

for(x in 1:n_sim){
  print(paste0("Iteration: ", x))
  
  # Generate the data, evaluate at desired points
  sim_dataset = gen_FPCA(FPCs, EVals, sim_domain, N, MuFun, sqrt(sigma2))
  {
    mObjs = true_fcs(MuFun, FPCs, weight_Obj$w, sim_domain)
    
    trueScore = data.frame(t(sim_dataset$Score_true))
    colnames(trueScore) = paste0("Curve ", 1:ncol(trueScore))
    trueScore$FPC = paste0("FPC ", 1:nrow(trueScore))
    trueScore = trueScore %>%
      pivot_longer(-c(FPC), names_to = "Curve", values_to = "Score")
    
    obsData = data.frame(t(sim_dataset$Y))
    colnames(obsData) = paste0("Curve ", 1:ncol(obsData))
    obsData$Arg = sim_domain
    obsData = obsData %>%
      pivot_longer(-c(Arg), names_to = "Curve", values_to = "Observed")
  }
  
  # FAST
  {
    # Format input data
    data_list = FAST_datalist(sim_dataset$Y, N, K, Q, sim_domain, weight_Obj$w,
                              basis_type = "Splinet", scale = F, alpha = alpha)
    
    # Fit model
    fast_mod = stan(
      file = "../FAST/FAST_SL.stan",
      data = data_list,
      chains = 4,
      cores = 4,
      warmup = 1000,
      iter = 1500,
      control = list(max_treedepth = 12), 
      verbose = F,
      refresh = 0
    )
    
    # Extract
    objects = FAST_extract(fast_mod, data_list$B, sim_domain, data_list)
    align = align_weights(objects$Weights, objects$Score, data_list$B, mObjs$mFPC)
    scores_samples = out_Score(align$Score)
    
    EF_CI = CI_EF(align$EF, sim_domain)
    EF_est = Psi_SVD_FPC(align$Weights, data_list$B, sim_domain, K, mObjs$mFPC)
    Mu_df = out_FE(objects$Mu, sim_domain) %>%
      group_by(Arg) %>%
      summarize(Est = mean(Mu), 
                LB = quantile(Mu, probs = c(0.025)), 
                UB = quantile(Mu, probs = c(0.975)))
    
    # Accuracy/validity measures
    FE_df =  FE_inf(mObjs$Mu, Mu_df, "FAST") %>%
      mutate(Sample = x, Alpha = alpha)
    EF_df = EF_inf(mObjs$FPC, inner_join(EF_CI, EF_est), "FAST") %>%
      mutate(Sample = x, Alpha = alpha)
    Score_df = Score_inf(trueScore, scores_samples, "FAST") %>%
      mutate(Sample = x, Alpha = alpha)
  }
  
  # Add most recent result
  FE_out = rbind(FE_out, FE_df)
  EF_out = rbind(EF_out, EF_df)
  Score_out = rbind(Score_out, Score_df)
  
  # Write to storage
  write.csv(FE_out, paste0(FE_dir, "/", args[2], ".csv"))
  write.csv(EF_out, paste0(EF_dir, "/", args[2], ".csv"))
  write.csv(Score_out, paste0(Score_dir,  "/", args[2], ".csv"))
}


