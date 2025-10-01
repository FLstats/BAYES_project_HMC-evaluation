/*	Shared slopes (no grouping variable) (NCP)
*	Grouped intercepts (grouping variable = job) (NCP)
*/

data {
	int N;
	int P;
	int<lower=1> N_g1;								// Number of levels in covariate g1
	array[N] int<lower=1,upper=N_g1> g1_id;			// Mapping from obs to level in g1
	int<lower=0,upper=1> y[N];
	matrix[N,P] X;									// (no intercept)
}

parameters {
	// Shared slopes
	vector[P] beta_raw;								// Auxilliary slope parameters							
	real<lower=0> tau;								// Shared slope sd	
	
	// Grouped intercepts
	vector[N_g1] alpha_raw_g1;						// Auxilliary group intercept parameters
	real mu_alpha_g1;								// Group intercept mean
	real<lower=0> sigma_alpha_g1;					// Group intercept sd
}

transformed parameters {
	// NCP for grouped intercepts
	vector[N_g1] alpha_g1	= mu_alpha_g1
							+ alpha_raw_g1
							* sigma_alpha_g1;
	
	// NCP for shared slopes
	vector[P] beta	= beta_raw
					* tau;
}

model {
	// Log-likelihood
	target += bernoulli_logit_lpmf(y | alpha_g1[g1_id] + X * beta);
	
	// Log-priors
	target += normal_lpdf(beta_raw | 0, 1);
	target += student_t_lpdf(tau | 3, 0, 2.5);
		
	target += normal_lpdf(alpha_raw_g1 | 0, 1);
	target += normal_lpdf(mu_alpha_g1 | 0, 10);
	target += student_t_lpdf(sigma_alpha_g1 | 3, 0, 2.5);	
}
