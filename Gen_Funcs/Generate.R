# Function to evaluate FPC bases
eval_Funcs <- function(func_list, domain){
  M = length(domain)
  N = length(func_list)
  out = matrix(0, nrow = M, ncol = N)
  for(n in 1:N){
    out[,n] = func_list[[n]](domain)
  }
  return(out)
}

# Function for GP deviations
gen_dev <- function(FPC_func, Lambdas, domain, N){
  M = length(domain)
  K = length(FPC_func)
  
  Phi = eval_Funcs(FPC_func, domain)
  Xi = matrix(0, nrow = N, ncol = K)
  for(k in 1:K){
    Xi[,k] = rnorm(N, sd = sqrt(Lambdas[k]))
  }
  
  return(list(Delta = Xi %*% t(Phi), Scores = Xi))
}

# One level of RE
gen_FPCA <- function(FPC_func, Lambdas, domain, N, mu_func, sigma_eps){
  M = length(domain)
  
  # Latent smooths
  GP_objs = gen_dev(FPC_func, Lambdas, domain, N)
  delta = GP_objs$Delta
  smooth = delta + matrix(rep(mu_func(domain), N), nrow = N, byrow = T)
  
  # Add noise
  outcome = smooth + matrix(rnorm(N*M, sd = sigma_eps), nrow = N, ncol = M)
  
  return(list(Y = outcome, Y_true = smooth, Score_true = GP_objs$Scores))
}

# Two levels of RE
gen_MFPCA <- function(FPC1, FPC2, Lambda1, Lambda2, domain, N, Subj, mu_func, sigma_eps){
  M = length(domain)
  I = n_distinct(Subj)
  
  # Latent smooths
  GP1 = gen_dev(FPC1, Lambda1, domain, I)
  delta1 = GP1$Delta
  GP2 = gen_dev(FPC2, Lambda2, domain, N)
  delta2 = GP2$Delta
  smooth = delta1[Subj, ] + delta2 + matrix(rep(mu_func(domain), N), nrow = N, byrow = T)
  
  # Add noise
  outcome = smooth + matrix(rnorm(N*M, sd = sigma_eps), nrow = N, ncol = M)
  
  return(list(Y = outcome, Y_true = smooth,
              Smooths1 = delta1, Smooths2 = delta2, 
              Scores1 = GP1$Scores, Scores2 = GP2$Scores))
}

# Realistic simulation based on CGM (DASH4D)
CGM_Funcs = function(){
  LP = polynomial.functions(legendre.polynomials(3, normalized = F))
  return(list(FPC1 = function(x){
                return(LP[[1]](2*x-1))
              }, 
              FPC2 = function(x){
                return(sqrt(84/31)*(LP[[2]](2*x-1) - 0.5*LP[[4]](2*x-1)))
              }, 
              FPC3 = function(x){
                return(-sqrt(5)*LP[[3]](2*x-1))
              }, 
              Mu = function(x){
                return(145 - 20*sqrt(2/5)*LP[[3]](2*x-1))
              }))
}
CGM_lambda = c(2250, 450, 150)
CGM_sigma2 = 4

# Canonical simulation (Fourier bases)
Can_Funcs = function(){
  return(list(FPC1 = function(x){
                        return(sqrt(2)*(sin(2 * pi * x)))
                      },
              FPC2 = function(x){
                        return(sqrt(2)*(cos(4 * pi * x)))
                      }, 
              FPC3 = function(x){
                        return(sqrt(2)*(sin(4 * pi * x)))
                      }, 
              Mu = function(x){
                        return(rep(0, length.out = length(x)))
                      }))
}
Can_lambda = 2^(0:-2)
Can_sigma2 = sum(Can_lambda)/5

# Collapse functions/evals/sigma2 into objects to pull from
Funcs = list(CGM = CGM_Funcs, Can = Can_Funcs)
Lambdas = list(CGM = CGM_lambda, Can = Can_lambda)
Sigma2 = list(CGM = CGM_sigma2, Can = Can_sigma2)
rm(CGM_Funcs, CGM_lambda, CGM_sigma2, 
   Can_Funcs, Can_lambda, Can_sigma2)

# Function for extracting Truth at various granularities
true_fcs <- function(MuFun, FPCs, weights, domain){
  tEFM = data.frame(eval_Funcs(FPCs, domain))
  colnames(tEFM) = paste0("FPC ", 1:ncol(tEFM))
  tEFM$Arg = domain
  tEFM$wei = weights
  tEFM = tEFM %>%
    pivot_longer(-c(Arg, wei), names_to = "FPC_Num", values_to = "Func")
  
  tMUM = data.frame(Arg = domain, Func = MuFun(domain), wei = weights)
  return(list(FPC = tEFM, Mu = tMUM, mFPC = eval_Funcs(FPCs, domain)))
}

ML_true_fcs <- function(MuFun, FPC1, FPC2, weights, domain){
  tEF1 = data.frame(eval_Funcs(FPC1, domain))
  colnames(tEF1) = paste0("FPC ", 1:ncol(tEF1))
  tEF1$Arg = domain
  tEF1$wei = weights
  tEF1 = tEF1 %>%
    pivot_longer(-c(Arg, wei), names_to = "FPC_Num", values_to = "Func")
  
  tEF2 = data.frame(eval_Funcs(FPC2, domain))
  colnames(tEF2) = paste0("FPC ", 1:ncol(tEF2))
  tEF2$Arg = domain
  tEF2$wei = weights
  tEF2 = tEF2 %>%
    pivot_longer(-c(Arg, wei), names_to = "FPC_Num", values_to = "Func")
  
  tMu = data.frame(Arg = domain, Func = MuFun(domain), wei = weights)
  return(list(FPC1 = tEF1, FPC2 = tEF2, Mu = tMu, 
              mFPC1 = eval_Funcs(FPC1, domain), mFPC2 = eval_Funcs(FPC2, domain)))
}
