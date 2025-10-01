data {
	int N;
	int P;
	int<lower=0,upper=1> y[N];
	matrix[N,P] X;
}

parameters {
	real alpha;
	vector[P] beta;
	real<lower=0> sigma2;
}

model {
	// Log-likelihood
	target += bernoulli_logit_lpmf(y | alpha + X * beta);
	
	// Log-priors
	target += normal_lpdf(alpha | 0, sqrt(sigma2));
	target += normal_lpdf(beta | 0, sqrt(sigma2));
	
	// Hyperprior
	target += exponential_lpdf(sigma2 | 0.01);
}
