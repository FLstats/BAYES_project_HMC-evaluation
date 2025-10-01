data {
	int<lower=1> N;
	int<lower=1> P;
	int<lower=0> sigma_v;
}

generated quantities {
	matrix[N, P+1] X;
	
	for(i in 1:N) {
		real v = normal_rng(0, sigma_v);
		X[i, 1] = v;
		
		for(k in 2:P+1) {
			X[i, k] = normal_rng(0, sqrt(exp(v)));
		}
	}
}
