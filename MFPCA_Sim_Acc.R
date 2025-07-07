# General helper functions
source("Gen_Funcs/Libs.R")
source("Gen_Funcs/Bases.R")
source("Gen_Funcs/Generate.R")
source("Gen_Funcs/PostProcess.R")
source("Gen_Funcs/Comparisons.R")

# Specific to each method
source("GFSR/GFSR_Help.R")
source("FAST/FAST_Help.R")

where_directory = "Results"

create_dir <- function(direct){
  if (!file.exists(direct)){
    dir.create(direct, recursive = T)
  }
}

args = commandArgs(trailingOnly=TRUE)
  
out_dir = args[1]

message("Simulation: ", args[1])
message("Output Directory: ", out_dir)
message("Output file: ", args[2])
message("Number of simulations: ", args[3])

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
I = 50
n_each = 5
ID = rep(1:I, each = n_each)
N = I*n_each
M = 50
weight_Obj = gaussLegendre(M, 0, 1)
sim_domain = weight_Obj$x

FPC1 = list(F1 = function(x) {
              return(sqrt(2) * sin(2*pi*x))
            },
            F2 = function(x) {
              return(sqrt(2) * cos(2*pi*x))
            }, 
            F3 = function(x) {
              return(sqrt(2) * sin(4*pi*x))
            },
            F4 = function(x) {
              return(sqrt(2) * cos(4*pi*x))
            })
FPC2 = list(F1 = function(x) {
              return(rep(1, length(x)))
            },
            F2 = function(x) {
              return(sqrt(3) * (2 * x - 1))
            }, 
            F3 = function(x) {
              return(sqrt(5) * (6*x^2 - 6*x + 1))
            },
            F4 = function(x) {
              return(sqrt(7)*(20*x^3 - 30*x^2 + 12*x - 1))
            })
MuFun = function(x){
  return(rep(0, length(x)))
}
K1 = length(FPC1)
K2 = length(FPC2)
EV1 = 0.5^(1:4 - 1)
EV2 = 0.5^(1:4 - 1)
sigma2 = 1
Q = 20 # Sufficiently rich bases, allow penalty to control smoothness

FE_out = data.frame(Method = c(), Mean_Cov = c(), MSE = c(), Sample = c())
EF_out = data.frame(Method = c(), EFNum = c(), Mean_Cov = c(), MSE = c(), Level = c(), Sample = c())
Score_out = data.frame(Method = c(), EFNum = c(), Mean_cov = c(), MSE = c(), Level = c(), Sample = c())

for(x in 1:n_sim){
  print(paste0("Iteration: ", x))
  
  # Generate the data, evaluate at desired points
  sim_dataset = gen_MFPCA(FPC1, FPC2, EV1, EV2, sim_domain, N, ID, MuFun, sqrt(sigma2))
  {
    mObjs = ML_true_fcs(MuFun, FPC1, FPC2, weight_Obj$w, sim_domain)
    
    tS1 = data.frame(t(sim_dataset$Scores1))
    colnames(tS1) = paste0("Curve ", 1:ncol(tS1))
    tS1$FPC = paste0("FPC ", 1:nrow(tS1))
    tS1 = tS1 %>%
      pivot_longer(-c(FPC), names_to = "Curve", values_to = "Score")
    
    tS2 = data.frame(t(sim_dataset$Scores2))
    colnames(tS2) = paste0("Curve ", 1:ncol(tS2))
    tS2$FPC = paste0("FPC ", 1:nrow(tS2))
    tS2 = tS2 %>%
      pivot_longer(-c(FPC), names_to = "Curve", values_to = "Score")
  }
  
  # FAST
  {
    # Format input data
    data_list = FAST_DL_ML(sim_dataset$Y, N, ID, K1, K2, Q, sim_domain,
                              basis_type = "Splinet", scale = F, alpha = 0.1)
    
    # Fit model
    fast_mod = stan(
      file = "FAST/FAST_ML.stan",
      data = data_list,
      chains = 4,
      cores = 4,
      warmup = 1000,
      iter = 2500,
      control = list(max_treedepth = 12),
      verbose = F,
      refresh = 0
    )
    
    # Extract
    objects = FAST_ML_extract(fast_mod, data_list$B, sim_domain, data_list)
    align1 = align_weights(objects$W1, objects$S1, data_list$B, mObjs$mFPC1)
    align2 = align_weights(objects$W2, objects$S2, data_list$B, mObjs$mFPC2)
    
    scores1 = out_Score(align1$Score)
    scores2 = out_Score(align2$Score)
    
    EF1_CI = CI_EF(align1$EF, sim_domain)
    EF1_est = Psi_SVD_FPC(align1$Weights, data_list$B, sim_domain, K1, mObjs$mFPC1)
    EF2_CI = CI_EF(align2$EF, sim_domain)
    EF2_est = Psi_SVD_FPC(align2$Weights, data_list$B, sim_domain, K2, mObjs$mFPC2)
    
    Mu_df = out_FE(objects$Mu, sim_domain) %>%
      group_by(Arg) %>%
      summarize(Est = mean(Mu), 
                LB = quantile(Mu, probs = c(0.025)), 
                UB = quantile(Mu, probs = c(0.975)))
    
    # Accuracy/validity measures
    fast_FE = FE_inf(mObjs$Mu, Mu_df, "FAST")
    fast_EF = rbind(
      EF_inf(mObjs$FPC1, inner_join(EF1_CI, EF1_est), "FAST") %>%
        mutate(Level = "Level 1"), 
      EF_inf(mObjs$FPC2, inner_join(EF2_CI, EF2_est), "FAST") %>%
        mutate(Level = "Level 2"))
    fast_S = rbind(
      Score_inf(tS1, scores1, "FAST") %>%
        mutate(Level = "Level 1"), 
      Score_inf(tS2, scores2, "FAST") %>%
        mutate(Level = "Level 2"))
  }
  message("FAST complete")
  
  # GFSR
  {
    # Format input data
    data_list = GFSR_DL_ML(sim_dataset$Y, I, n_each, ID, K1, K2, Q, sim_domain)
    
    # Fit the model
    regress_mod = stan(
      file = "GFSR/GFSR_ML.stan",
      data = data_list,
      chains = 4, 
      cores = 4,
      warmup = 1000,
      iter = 2500,
      control = list(max_treedepth = 12),
      verbose = F,
      refresh = 0
    )
    
    # Extract - higher resolution
    objects = GFSR_ML_extract(regress_mod, data_list$BS, sim_domain, data_list)
    align1 = align_FPCs(objects$EF1, objects$S1, mObjs$mFPC1)
    align2 = align_FPCs(objects$EF2, objects$S2, mObjs$mFPC2)
    
    scores1 = out_Score(align1$Score)
    scores2 = out_Score(align2$Score)
    
    EF1_CI = CI_EF(align1$EF, sim_domain)
    EF1_est = Phi_SVD_FPC(align1$EF, sim_domain, K1, mObjs$mFPC1)
    EF2_CI = CI_EF(align2$EF, sim_domain)
    EF2_est = Phi_SVD_FPC(align2$EF, sim_domain, K2, mObjs$mFPC2)
    
    Mu_df = out_FE(objects$Mu, sim_domain) %>%
      group_by(Arg) %>%
      summarize(Est = mean(Mu), 
                LB = quantile(Mu, probs = c(0.025)), 
                UB = quantile(Mu, probs = c(0.975)))
    
    # Accuracy/validity measures
    gfsr_FE =  FE_inf(mObjs$Mu, Mu_df, "GFSR")
    gfsr_EF = rbind(
      EF_inf(mObjs$FPC1, inner_join(EF1_CI, EF1_est), "GFSR") %>%
        mutate(Level = "Level 1"), 
      EF_inf(mObjs$FPC2, inner_join(EF2_CI, EF2_est), "GFSR") %>%
        mutate(Level = "Level 2"))
    gfsr_S = rbind(
      Score_inf(tS1, scores1, "GFSR") %>%
        mutate(Level = "Level 1"), 
      Score_inf(tS2, scores2, "GFSR") %>%
        mutate(Level = "Level 2"))
  }
  message("GFSR complete")
  
  # Collate accuracy measures
  {
    FE_df = rbind(fast_FE, gfsr_FE) %>% 
      mutate(Sample = x)
    EF_df = rbind(fast_EF, gfsr_EF) %>%
      mutate(Sample = x)
    Score_df = rbind(fast_S, gfsr_S) %>%
      mutate(Sample = x)
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


