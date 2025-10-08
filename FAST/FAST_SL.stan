data {
  int N;    // Number of time series
  int M;    // Number of observations
  int Q;    // Number of spline bases
  int K;    // Number of Eigenfunctions
  
  matrix[N, M] Y;         // Original data
  matrix[M, Q] B;         // Orthogonal functional basis
  matrix[Q, Q] P_alpha;   // Penalty matrix
}

transformed data {
  real tr_P = trace(P_alpha);  // Trace of penalty
}

parameters {
  real<lower=0> sigma2; // Error in observation
  
  // Fixed-effect components
  vector[Q] w_mu;               // Population mean parameters
  real<lower=0> h_mu;           // Population mean smoothing parameter
  
  // Components/weights
  positive_ordered[K] lambda;        // Eigenvalues
  vector<lower=0>[K] H;              // EF Smoothing parameters
  matrix[Q, K] X;                    // Unconstrained EF weights (X)
  matrix[N, K] Scores;               // EF scores (xi)
}

transformed parameters{
  // Fixed effects - population mean
  row_vector[M] mu = (B * w_mu)';
  
  // Polar decomposition
  matrix[Q, K] Psi;
  {
    matrix[K,K] evec_XtX = eigenvectors_sym(X'*X); 
    vector[K] eval_XtX = eigenvalues_sym(X'*X);
    Psi = X*evec_XtX*diag_matrix(1/sqrt(eval_XtX))*evec_XtX';
  }
}

model {
  
  // Smoothing weight priors 
  H ~ gamma(0.01, tr_P/2 + 0.01);
  h_mu ~ gamma(0.001, 0.001);
  
   // Variance component priors
  lambda ~ inv_gamma(0.001, 0.001); 
  sigma2 ~ inv_gamma(0.001, 0.001);
  
  // Smoothing additions to the target density
  target += Q/2.0 * log(h_mu) - h_mu / 2.0 * w_mu' * P_alpha * w_mu;
  for(i in 1:K){
    target += Q/2.0 * log(H[i]) - H[i] / 2.0 * Psi[,i]' * P_alpha * Psi[,i]; 
  }
  
  // Uniform priors through matrix normals
  to_vector(X) ~ normal(0, 1);
  
  // Score priors 
  for(k in 1:K){
    to_vector(Scores[,k]) ~ normal(0, sqrt(lambda[k])); 
  }
  
  // Likelihood
  {
    // Smooth component estimates
    matrix[N, M] Theta = Scores * (B * Psi)';
    
    to_vector(Y) ~ normal(to_vector(rep_matrix(mu, N) + Theta), sqrt(sigma2));
  }
}
