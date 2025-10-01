data  {
	int N;
	int<lower=2> ni;
	int<lower=1> nj;
	int<lower=1> n_dim;				// Total dimension
	matrix[N, n_dim] X;				// Observed (simulated) data
}

parameters {
	real mu;
	real<lower=0> a;
	real<lower=0> b;
}

model {
	/* For each observation n, 
	*	First evaluate the likelihood for x1
	*	Then evaluate the likelihood for the following variables in the block.
	*
	*	Each observation n is a vector of length n_dim.
	*/

	// Log Likelihood
	for(n in 1:N) {
		// Root variable x_[j,1]
		target += normal_lpdf(X[n, 1] | mu, inv_sqrt(2 * a));
		
		// For each block j, the first variable depends on x1. 
		for(j in 1:nj) {
			int idx_first = 2 + (j-1) * (ni-1);		// index of x_[j,2]
			target += normal_lpdf(X[n, idx_first] | square(X[n, 1]), inv_sqrt(2 * b));
		
			// Then each variable depends on *previous variable in block* squared.
			for(i in 3:ni) {
				int idx_cur = idx_first + (i-2);		// index of x_[j,i]
				int idx_prev = idx_cur - 1;		// index of x_[j,i-1]
				target += normal_lpdf(X[n, idx_cur] | square(X[n, idx_prev]), inv_sqrt(2 * b));
			}
		}
	}
	
	// Log Priors
	target += normal_lpdf(mu | 0, 10);
	target += student_t_lpdf(a | 3, 0, 2.5);
	target += student_t_lpdf(b | 3, 0, 2.5);
}
