data {
  int N;
  int P;
  int<lower=1> N_g;  // Number of levels in covariate g
  array[N] int<lower=1,upper=N_g> g_id;  // Mapping from obs to level in g
  array[N] int<lower=0,upper=1> y;
  matrix[N,P] X;
}

parameters {
  // Shared slope vector for all obs
  vector[P] beta;
  real<lower=0> sigma2_beta;

  // Group intercepts
  vector[N_g] alpha_g;
  real<lower=0> sigma2_alpha_g;
}

model {
  // Log-likelihood
  target += bernoulli_logit_lpmf(y | alpha_g[g_id] + X * beta);

  // Log-priors
  target += normal_lpdf(beta | 0, sqrt(sigma2_beta));
  target += exponential_lpdf(sigma2_beta | 0.01);
	
  target += normal_lpdf(alpha_g | 0, sqrt(sigma2_alpha_g));
  target += exponential_lpdf(sigma2_alpha_g | 0.01);
}
