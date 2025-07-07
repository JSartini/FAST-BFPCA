
# Produce data/hyper-priors for GFSR
GFSR_datalist <- function(Y, N, K, Q, Domain, alpha = 0.1){
  
  # FE matrix - just a population intercept
  X_mat = matrix(0, nrow = N, ncol = 1)
  X_mat[,1] = 1
  
  # Splines
  BS = bSpline(Domain, df = Q, intercept = T)
  
  # Penalty matrix - taken from Goldsmith et al.
  D = length(Domain)
  diff0 = diag(1, D, D)
  diff2 = matrix(rep(c(1,-2,1, rep(0, D-2)), D-2)[1:((D-2)*D)], D-2, D, byrow = TRUE)
  P0 = t(BS) %*% t(diff0) %*% diff0 %*% BS
  P2 = t(BS) %*% t(diff2) %*% diff2 %*% BS
  P.mat = alpha * P0 + (1-alpha) * P2
  
  # Format for GFSR
  stan_dat = list(I = N, D = length(Domain), p = 1, Kt = Q, Kp = K, Y = Y,
                     X = X_mat, BS = BS, PenMat = P.mat)
  return(stan_dat)
}

GFSR_DL_ML <- function(Y, I, N_Each, IDs, K1, K2, Q, Domain, alpha = 0.1){
  
  # FE matrix - just a population intercept
  X_mat = matrix(0, nrow = I*N_Each, ncol = 1)
  X_mat[,1] = 1
  
  # Splines
  BS = bSpline(Domain, df = Q, intercept = T)
  
  # Penalty matrix - taken from Goldsmith et al.
  D = length(Domain)
  diff0 = diag(1, D, D)
  diff2 = matrix(rep(c(1,-2,1, rep(0, D-2)), D-2)[1:((D-2)*D)], D-2, D, byrow = TRUE)
  P0 = t(BS) %*% t(diff0) %*% diff0 %*% BS
  P2 = t(BS) %*% t(diff2) %*% diff2 %*% BS
  P.mat = alpha * P0 + (1-alpha) * P2
  
  # Format for GFSR
  stan_dat = list(I = I, J = N_Each, IJ = I * N_Each, D = length(Domain), p = 1, Kt = Q,
                  Kp1 = K1, Kp2 = K2, Y = Y, X = X_mat, BS = BS, PenMat = P.mat)
  return(stan_dat)
  
  return()
}

# Place EF and scores in standard form for post-processing
GFSR_extract <- function(mod_fit, BS, M, K){
  samples = extract(mod_fit)
  n_samp = length(samples$sigma2)
  
  smooth_list = map(1:n_samp, function(x){
    smooth_sample = BS %*% t(samples$beta_psi[x,,])
    return(svd(t(smooth_sample)))
  })
  
  EF_list = map(smooth_list, function(svd_x){
    return(svd_x$v[,1:K]*sqrt(M))
  })
  
  Score_list = map(1:n_samp, function(x){
    U = smooth_list[[x]]$u[,1:K]
    D = diag(smooth_list[[x]]$d[1:K])
    Score_sample = samples$c[x,,] %*% U %*% D / sqrt(M)
    return(Score_sample)
  })
  
  Mu_list = map(1:n_samp, function(x){
    return(BS %*% samples$beta[x,,])
  }) 
  
  return(list(EF = EF_list, Score = Score_list, Mu = Mu_list))
}

GFSR_ML_extract <- function(mod_fit, BS, Domain, DL){
  samples = extract(mod_fit)
  n_samp = length(samples$sigma2)
  
  smooths1 = map(1:n_samp, function(x){
    smooth_sample = BS %*% t(samples$beta_psi1[x,,])
    return(svd(t(smooth_sample)))
  })
  
  smooths2 = map(1:n_samp, function(x){
    smooth_sample = BS %*% t(samples$beta_psi2[x,,])
    return(svd(t(smooth_sample)))
  })
  
  EFs1 = map(smooths1, function(svd_x){
    return(svd_x$v[,1:DL$Kp1]*sqrt(DL$D))
  })
  
  EFs2 = map(smooths2, function(svd_x){
    return(svd_x$v[,1:DL$Kp2]*sqrt(DL$D))
  })
  
  Scores1 = map(1:n_samp, function(x){
    U = smooths1[[x]]$u[,1:DL$Kp1]
    D = diag(smooths1[[x]]$d[1:DL$Kp1])
    Score_sample = samples$c1[x,,] %*% U %*% D / sqrt(DL$D)
    return(Score_sample)
  })
  
  Scores2 = map(1:n_samp, function(x){
    U = smooths2[[x]]$u[,1:DL$Kp2]
    D = diag(smooths2[[x]]$d[1:DL$Kp2])
    Score_sample = samples$c2[x,,] %*% U %*% D / sqrt(DL$D)
    return(Score_sample)
  })
  
  Mu_list = map(1:n_samp, function(x){
    return(BS %*% samples$beta[x,,])
  }) 
  
  return(list(EF1 = EFs1, EF2 = EFs2, S1 = Scores1, S2 = Scores2, Mu = Mu_list))
}

# Extract EF and Scores by-chain
GFSR_byChain <- function(mod_fit, N, M, K, B){
  chain_samples = extract(mod_fit, permuted = F)
  dimen_names = dimnames(chain_samples)$parameters
  
  n_samp = dim(chain_samples)[1]
  n_chain = dim(chain_samples)[2]
  
  Psi_idx = grepl("beta_psi", dimen_names)
  Score_idx = grepl("c", dimen_names)
  
  Phi_byChain = list()
  Score_byChain = list()
  
  for(j in 1:n_chain){
    Phi_byChain[[j]] = map(1:n_samp, function(i){
      Psi = chain_samples[i,j,Psi_idx]
      dim(Psi) = c(Q, K)
      Phi = svd(t(B %*% Psi))$v[, 1:K]
      return(Phi * sqrt(M))
    })
    Score_byChain[[j]] = map(1:n_samp, function(i){
      Psi = chain_samples[i,j,Psi_idx]
      dim(Psi) = c(Q, K)
      svd_objs = svd(t(B %*% Psi))
      
      Xi = chain_samples[i,j,Score_idx] 
      dim(Xi) = c(N, K)
      Scores = Xi %*% svd_objs$u[,1:K] %*% diag(svd_objs$d[1:K])
      return(Scores / sqrt(M))
    })
  }
  return(list(Phi = Phi_byChain, Score = Score_byChain))
}
