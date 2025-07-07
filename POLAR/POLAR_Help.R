
# Produce data/hyper-priors for POLAR
POLAR_datalist <- function(crossings, ratio, K, Y, N, p, domain){
  Yp = sweep(Y, 2, colMeans(Y), "-")
  
  expectation = 1/(2*pi*crossings)
  stddev = expectation*ratio
  
  solve_fn <- function(x){
    E <- x[2]/(x[1] - 1) - expectation
    V <- x[2]^2/((x[1] - 1)^2 * (x[1] - 2)) - stddev^2
    return(c(E, V))
  }
  
  smooth_priors = nleqslv(c(5, 5), solve_fn)$x
  alpha = smooth_priors[1]
  beta = smooth_priors[2]
  
  approx_comps = svd(Y, nv = K)
  best_approx = approx_comps$u[,1:K] %*% diag(approx_comps$d[1:K]) %*% t(approx_comps$v)
  nu = 1
  s2 = 3*var(as.numeric(Y - best_approx))
  
  tau2 = sum(diag(t(best_approx) %*% best_approx))/K
  return(list(n = N, p = p, Y = Yp, k = K, days = domain, nu = nu,
              s2 = s2, tau = sqrt(tau2), alpha = alpha, beta = beta))
}

# Place EF and scores in standard form for post-processing
POLAR_extract <- function(mod_fit){
  samples = extract(mod_fit)
  n_samp = length(samples$sig2)
  M = dim(samples$Phi_s)[2]
  
  EF_list = map(1:n_samp, function(x){
    EF_sample = samples$Phi_s[x,,]*sqrt(M)
    return(EF_sample)
  })
  
  Score_list = map(1:n_samp, function(x){
    Score_sample = samples$Scores[x,,]/sqrt(M)
    return(Score_sample)
  })
  
  return(list(EF = EF_list, Score = Score_list))
}

# Extract EF and Scores by-chain
POLAR_byChain <- function(mod_fit, N, M, K){
  chain_samples = extract(mod_fit, permuted = F)
  dimen_names = dimnames(chain_samples)$parameters
  
  n_samp = dim(chain_samples)[1]
  n_chain = dim(chain_samples)[2]
  
  Phi_idx = grepl("Phi_s", dimen_names)
  Score_idx = grepl("Scores", dimen_names)
  
  Phi_byChain = list()
  Score_byChain = list()
  
  for(j in 1:n_chain){
    Phi_byChain[[j]] = map(1:n_samp, function(i){
      Phi = chain_samples[i,j,Phi_idx]
      dim(Phi) = c(M, K)
      return(Phi)
    })
    Score_byChain[[j]] = map(1:n_samp, function(i){
      Xi = chain_samples[i,j,Score_idx]
      dim(Xi) = c(N, K)
      return(Xi)
    })
  }
  return(list(Phi = Phi_byChain, Score = Score_byChain))
}
