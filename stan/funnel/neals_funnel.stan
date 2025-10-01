data {
	int N;							// Number of observations
	int P;							// Number of parameters (v & x's)
	matrix[N, P] X;
}

parameters {
	real<lower=0> v;
	real<lower=0> sigma_v;
}

model {
	/*	Log likelihood
	*	X[n] is a row vector. Stan will go through each element (each x_k).
	*/
	
	for(n in 1:N) {
		target += normal_lpdf(X[n] | 0, sqrt(exp(v)));
	}
	
	// Log priors
	target += normal_lpdf(v | 0, sigma_v);
	target += student_t_lpdf(sigma_v | 3, 0, 2.5);
}
