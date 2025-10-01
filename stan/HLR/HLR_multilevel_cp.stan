data {
	int N;
	int P;
	int<lower=1> N_C1;									// Number of levels in covariate C1
	array[N] int<lower=1,upper=N_C1> C1_id;					// Mapping from obs to level in C1
	int<lower=0,upper=1> y[N];
	matrix[N,P] X;
}

parameters {
	vector[P] beta;										// Shared slope vector for all obs
	
	vector[N_C1] alpha_C1;
	real mu_alpha_C1;
	real<lower=0> sigma_alpha_C1;
	
	real<lower=0> sigma2;
}

transformed parameters {
	vector[N] alpha_obs;							// Every obs mapped to its intercept level
	for(n in 1:N) {
		alpha_obs[n] = alpha_C1[C1_id[n]];
	}
}

model {
	// Log-likelihood
	target += bernoulli_logit_lpmf(y | alpha_obs + X * beta);
	
	// Log-priors
	target += student_t_lpdf(beta | 3, 0, sqrt(sigma2));
		
	target += normal_lpdf(alpha_C1 | mu_alpha_C1, sigma_alpha_C1);
	target += normal_lpdf(mu_alpha_C1 | 0, 10);
	target += student_t_lpdf(sigma_alpha_C1 | 3, 0, 2.5);

	// Hyperprior
	target += exponential_lpdf(sigma2 | 0.01);
}
