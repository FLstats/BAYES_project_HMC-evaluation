/*	Shared slopes (CP)
*	Grouped intercepts (job) (NCP)
*/

data {
	int N;
	int P;
	int<lower=1> N_C1;									// Number of levels in covariate C1
	array[N] int<lower=1,upper=N_C1> C1_id;				// Mapping from obs to level in C1
	int<lower=0,upper=1> y[N];
	matrix[N,P] X;										// (no intercept)
}

parameters {
	// Shared slopes
	vector[P] beta;										// Shared slope vector for all obs
	real<lower=0> sigma2;								// Slope variance
	
	// Grouped intercepts
	vector[N_C1] alpha_raw_C1;						// One intercept per level in C1
	real mu_alpha_C1;
	real<lower=0> sigma_alpha_C1;
}

transformed parameters {
	// NCP for grouped intercepts
	vector[N_C1] alpha_C1	= mu_alpha_C1
							+ alpha_raw_C1
							* sigma_alpha_C1;
}

model {
	// Log-likelihood
	target += bernoulli_logit_lpmf(y | alpha_C1[C1_id] + X * beta);
	
	// Log-priors
	target += student_t_lpdf(beta | 3, 0, sqrt(sigma2));
		
	target += normal_lpdf(alpha_raw_C1 | 0, 1);
	target += normal_lpdf(mu_alpha_C1 | 0, 10);
	target += student_t_lpdf(sigma_alpha_C1 | 3, 0, 2.5);

	// Hyperprior
	target += exponential_lpdf(sigma2 | 0.01);
}
