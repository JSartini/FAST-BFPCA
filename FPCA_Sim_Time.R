# General helper functions
source("Gen_Funcs/Libs.R")
source("Gen_Funcs/Bases.R")
source("Gen_Funcs/Generate.R")
source("Gen_Funcs/PostProcess.R")
source("Gen_Funcs/Comparisons.R")

# Specific to each method
source("POLAR/POLAR_Help.R")
source("GFSR/GFSR_Help.R")
source("FAST/FAST_Help.R")
source("VMP/VMP_Help.R")

where_directory = "Timing_Results"

create_dir <- function(direct){
  if (!file.exists(direct)){
    dir.create(direct, recursive = T)
  }
}

args = commandArgs(trailingOnly=TRUE)
out_dir = paste0(args[1], "_N", args[4], "_M", args[5])

message("Simulation: ", args[1])
message("Output Directory: ", out_dir)
message("Output file: ", args[2])
message("Number of simulations: ", args[3])
message("Number of time series: ", args[4])
message("Density of grid: ", args[5])

# Create output directories
{
  Time_dir = paste0(where_directory, "/Time/", out_dir)
  create_dir(Time_dir)
}

# Begin simulations
n_sim = as.numeric(args[3])
N = as.numeric(args[4])
M = as.numeric(args[5])
weight_Obj = gaussLegendre(M, 0, 1)
sim_domain = weight_Obj$x

FComps = Funcs[[args[1]]]()
FPCs = FComps[-length(FComps)]
MuFun = FComps[[length(FComps)]]
EVals = Lambdas[[args[1]]]
sigma2 = Sigma2[[args[1]]]
Q = 20 # Sufficiently rich bases, allow penalty to control smoothness
K = length(FPCs)

Time_out = data.frame(Method = c(), Seconds = c(), Sample = c())

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
    start.time = Sys.time()

    # Format input data
    data_list = FAST_datalist(sim_dataset$Y, N, K, Q, sim_domain, weight_Obj$w,
                              basis_type = "Splinet", scale = F, alpha = 0.1)

    # Fit model
    fast_mod = stan(
      file = "FAST/FAST_SL.stan",
      data = data_list,
      chains = 4,
      cores = 4,
      warmup = 1000,
      iter = 1500,
      control = list(max_treedepth = 12),
      verbose = F,
      refresh = 0
    )

    # Extract and align
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

    end.time = Sys.time()
    fastS_time = difftime(end.time, start.time, units = "secs") %>% as.numeric()
  }
  message("FAST complete")
  
  # GFSR
  {
    start.time = Sys.time()

    # Format input data
    data_list = GFSR_datalist(sim_dataset$Y, N, K, Q, sim_domain)

    # Fit the model
    regress_mod = stan(
      file = "GFSR/GFSR_SL.stan",
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
    objects = GFSR_extract(regress_mod, data_list$BS, M, K)
    align = align_FPCs(objects$EF, objects$Score, mObjs$mFPC)
    scores_samples = out_Score(align$Score)

    EF_CI = CI_EF(align$EF, sim_domain)
    EF_est = Phi_SVD_FPC(align$EF, sim_domain, K, mObjs$mFPC)
    Mu_df = out_FE(objects$Mu, sim_domain) %>%
      group_by(Arg) %>%
      summarize(Est = mean(Mu),
                LB = quantile(Mu, probs = c(0.025)),
                UB = quantile(Mu, probs = c(0.975)))

    end.time = Sys.time()
    gfsr_time = difftime(end.time, start.time, units = "secs") %>% as.numeric()
  }
  message("GFSR complete")
  
  # POLAR
  {
    start.time = Sys.time()

    # Format input data
    data_list = POLAR_datalist(2, 1/6, K, sim_dataset$Y, N, M, sim_domain)

    # Fit model
    polar_mod = stan(
      file = "POLAR/POLAR.stan",
      data = data_list,
      chains = 4,
      cores = 4,
      warmup = 1000,
      iter = 1500,
      control = list(max_treedepth = 12),
      verbose = F,
      refresh = 0
    )

    # Extract results
    objects = POLAR_extract(polar_mod)
    align = align_FPCs(objects$EF, objects$Score, anchor = mObjs$mFPC)

    scores_samples = out_Score(align$Score)
    EF_CI = CI_EF(align$EF, sim_domain)
    EF_est = Smooth_SVD_FPC(out_Devs(align$EF, align$Score), sim_domain, K, mObjs$mFPC)

    end.time = Sys.time()
    polar_time = difftime(end.time, start.time, units = "secs") %>% as.numeric()
  }
  message("POLAR complete")
  
  # VMP
  {
    start.time = Sys.time()

    # Format input data
    vfit = VMP_fit_objs(sim_dataset, N, sim_domain)

    # Fit model - higher resolution
    vmp_mod = run_vmp_fpca(time_obs = vfit$Tl, Y = vfit$Yl, time_g = sim_domain,
                           L = K, verbose = F, Psi_g = mObjs$mFPC)

    # Extract results
    vmp_ests = VMP_extract(vmp_mod, sim_domain)
    EF_est = vmp_ests$EF
    EF_CI = EF_est %>%
      mutate(LB = NA, UB = NA) %>%
      select(-c(FPC_Val))
    scores_samples = vmp_ests$Score
    Mu_df = vmp_ests$Mu %>%
      mutate(LB = NA, UB = NA)

    end.time = Sys.time()
    vmp_time = difftime(end.time, start.time, units = "secs") %>% as.numeric()
  }
  message("VMP complete")
  
  # Collate timing
  {
    Time_df = data.frame(Seconds = c(fastS_time, vmp_time,
                                      gfsr_time, polar_time), 
                         Method = c("FAST", "VMP", "GFSR", "POLAR"),
                         Sample = x)
  }
  
  # Add most recent result
  Time_out = rbind(Time_out, Time_df)
  
  # Write to storage
  write.csv(Time_out, paste0(Time_dir, "/", args[2], ".csv"))
}