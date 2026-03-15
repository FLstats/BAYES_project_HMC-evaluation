/* Shared slopes (no grouping variable) (NCP)
*  Grouped intercepts (grouping variable = job) (NCP)
*/

data {
  int N;
  int P;
  int<lower=1> N_g;  // Number of levels in covariate g
  array[N] int<lower=1,upper=N_g> g_id;  // Mapping from obs to level in g
  array[N] int<lower=0,upper=1> y;
  matrix[N,P] X;  // (no intercept)
}

parameters {
  // Shared slope vector for all obs
  vector[P] beta_raw;							
  real<lower=0> sigma2_beta;

  // Group intercepts
  vector[N_g] alpha_g_raw;
  real<lower=0> sigma2_alpha_g;
}

transformed parameters {
  vector[P] beta	= beta_raw
					* sigma2_beta;

  vector[N_g] alpha_g	= alpha_g_raw
						* sigma2_alpha_g;
}

model {
  // Log-likelihood
  target += bernoulli_logit_lpmf(y | alpha_g[g_id] + X * beta);

  // Log-priors
  target += normal_lpdf(beta_raw | 0, 1);
  target += exponential_lpdf(sigma2_beta | 0.01);
	
  target += normal_lpdf(alpha_g_raw | 0, 1);
  target += exponential_lpdf(sigma2_alpha_g | 0.01);	
}
