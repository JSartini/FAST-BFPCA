////////////////////////////////////////////////////////////
// February 1 2014
// Jeff Goldsmith
//
// This file contains code for cross-sectional bayesian 
// generalized function-on-scalar regression with FPCA for 
// functional residuals. 
////////////////////////////////////////////////////////////

data {
	int<lower=0> I;                           // number of subjects
	int<lower=0> D;                           // grid length
	int<lower=0> p;                           // number of fixed effects

	int Kt;                                   // number of spline basis functions
	int Kp;                                   // number of PC basis functions

	matrix[I,D] Y;                            // outcome matrix
	vector[p] X[I];                           // fixed effect design matrix
	matrix[D,Kt] BS;                          // B-spline evaluation matrix
	
  cov_matrix[Kt] PenMat;                    // prior precision matrix for spline effects (penalty matrix)
}

transformed data {
  vector[Kt] mu_beta;                       // prior mean for spline effects

	for (k in 1:Kt) {
		mu_beta[k] = 0;
	}
}


parameters {
	matrix[p,Kt] beta;	                       // matrix of fixed effect spline coefficients
	matrix[Kp,Kt] beta_psi;	                   // matrix of PC spline coefficients
	vector[Kp] c[I];	                         // matrix of PC loadings
	vector<lower=0,upper=100>[p] beta_sig;     // tuning parameter
	vector<lower=0,upper=100>[Kp] psi_sig;     // tuning parameter
	real<lower=0> sigma2;                      // noise variance
}

transformed parameters {
	vector<lower=0.01>[p] beta_tau2;    // tuning parameter
	vector<lower=0.01>[Kp] psi_tau2;    // tuning parameter

	for(pcur in 1:p) {
		beta_tau2[pcur] = pow(beta_sig[pcur], -1);
	}

	for(kcur in 1:Kp) {
		psi_tau2[kcur] = pow(psi_sig[kcur], -1);
	}	
}


model {
	/////////////////////////////////////
	// Prior distributions
	/////////////////////////////////////
  
  // Prior for noise variance
  sigma2 ~ inv_gamma(0.001, 0.001);
  
	// Prior for variance components controlling smoothness in beta
	for(pcur in 1:p) {
		beta_sig[pcur] ~ inv_gamma(.001,.001);
//		beta_sig[pcur] ~ uniform(0,25); 
	}
	
	// Prior for variance components controlling smoothness in psi
	for(kcur in 1:Kp) {
		psi_sig[kcur] ~ inv_gamma(.001,.001);
//		psi_sig[kcur] ~ uniform(0,25); 
	}

	// Prior for spline coefficients for beta
	for (pcur in 1:p) {
		(beta[pcur])' ~ multi_normal_prec(mu_beta, beta_tau2[pcur] * PenMat);
	}

	// Prior for spline coefficients for psi
	for (kcur in 1:Kp) {
		(beta_psi[kcur])' ~ multi_normal_prec(mu_beta, psi_tau2[kcur] * PenMat);
	}
	
	// Prior for curve-level PC scores
	for (i in 1:I) {
		c[i] ~ normal(0.0, 1.0);
	}

	/////////////////////////////////////
	// Outcome likelihood
	/////////////////////////////////////

	for (i in 1:I) {
		Y[i] ~ normal((BS * beta') * X[i] +                     // fixed effects 
		              (BS * beta_psi') * c[i], sqrt(sigma2));   // PC effects
	}
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
