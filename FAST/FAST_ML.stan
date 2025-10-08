data {
  int N;         // Total number of functions
  int I;         // Number of groups
  int ID[N];     // Index-based ID per function
  int M;         // Number of observations along the domain

  int Q;    // Dimension of spline bases
  int K1;   // Number of Eigenfunctions at group level
  int K2;   // Number of Eigenfunctions at visit level
  
  matrix[N, M] Y;     // Original data
  matrix[M, Q] B;     // Orthogonalized basis
  matrix[Q, Q] P_alpha;     // Penalty matrix for splines
}

transformed data {
  real tr_P = trace(P_alpha);  // Trace of penalty
}

parameters {
  real<lower=0> sigma2; // Error in observation
  
  // Fixed-effect components
  vector[Q] w_mu;       // B-spline weights for overall mean
  real<lower=0> H_mu;   // Smoothness parameter for overall mean
  
  // Person-specific components
  positive_ordered[K1] lambda_1;  // Variances of each eigen-function
  vector<lower=0>[K1] H_1;        // Smoothing multipliers
  matrix[Q, K1] X_1;              // Full-rank representation decomposed into group EF
  matrix[I, K1] Scores1;      // Raw group scores
  
  // Person-visit components
  positive_ordered[K2] lambda_2;  // Variances of each eigen-function
  vector<lower=0>[K2] H_2;        // Smoothing multiplier
  matrix[Q, K2] X_2;              // Full-rank representation decomposed into visit EF
  matrix[N, K2] Scores2;      // Raw instance scores
}

transformed parameters{
  // Fixed effects
  row_vector[M] mu = (B * w_mu)';
  
  // Polar decomposition
  matrix[Q, K1] Psi_1;
  matrix[Q, K2] Psi_2;
  {
    vector[K1] trans_eval_V1;   // Inverse sqrt of eigenvalues of X_1'*X_1
    matrix[K1,K1] evec_V1;      // Eigenvectors of X_1'*X_1
    
    trans_eval_V1 = 1/sqrt(eigenvalues_sym(X_1'*X_1));
    evec_V1 = eigenvectors_sym(X_1'*X_1); 
    Psi_1 = X_1*evec_V1*diag_matrix(trans_eval_V1)*evec_V1';
    
    vector[K2] trans_eval_V2;   // Inverse sqrt of eigenvalues of X_2'*X_2
    matrix[K2,K2] evec_V2;      // Eigenvectors of X_2'*X_2
    
    trans_eval_V2 = 1/sqrt(eigenvalues_sym(X_2'*X_2));
    evec_V2 = eigenvectors_sym(X_2'*X_2); 
    Psi_2 = X_2*evec_V2*diag_matrix(trans_eval_V2)*evec_V2';
  }
}

model {
  // Smoothing weight priors 
  H_mu ~ gamma(0.001, 0.001); 
  H_1 ~ gamma(0.01, tr_P/2 + 0.01);
  H_2 ~ gamma(0.01, tr_P/2 + 0.01);
  
  // Variance component priors
  lambda_1 ~ inv_gamma(0.001, 0.001);
  lambda_2 ~ inv_gamma(0.001, 0.001);
  sigma2 ~ inv_gamma(0.001, 0.001);
  
  // Smoothing additions to the target density
  target += Q/2.0 * log(H_mu) - H_mu / 2.0 * w_mu' * P_alpha * w_mu;
  for(i in 1:K1){
    target += Q/2.0 * log(H_1[i]) - H_1[i] / 2.0 * Psi_1[,i]' * P_alpha * Psi_1[,i];
  }
  for(i in 1:K2){
    target += Q/2.0 * log(H_2[i]) - H_2[i] / 2.0 * Psi_2[,i]' * P_alpha * Psi_2[,i];
  } 
  
  // Uniform priors through matrix normals
  to_vector(X_1) ~ normal(0, 1);
  to_vector(X_2) ~ normal(0, 1);
  
  // Score priors
  for(k in 1:K1){
    to_vector(Scores1[,k]) ~ normal(0, sqrt(lambda_1[k]));
  }
  for(k in 1:K2){
    to_vector(Scores2[,k]) ~ normal(0, sqrt(lambda_2[k]));
  }
  
  // Likelihood
  {
    matrix[I, M] UDV_1 = Scores1 * (B * Psi_1)';
    matrix[N, M] UDV_2 = Scores2 * (B * Psi_2)';
    to_vector(Y) ~ normal(to_vector(rep_matrix(mu, N) + UDV_1[ID] + UDV_2), sqrt(sigma2));
  }
}
