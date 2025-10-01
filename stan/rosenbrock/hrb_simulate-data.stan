/*	Implements the R simulation of the n-dimensional Hybrid Rosenbrock.
* 	Sampling without parameters (all output quantities are produced
*		in the generated quantities block).
*/

data {
	int<lower=1> N;						// Number of samples
	int<lower=1> ni;					// Index "inside block"
	int<lower=1> nj;					// Index "block"
	real mu;							// Mean
	real<lower=0> a;					// sd-factor for root x1
	real<lower=0> b;					// Shared sd-factor for following variables
}

transformed data {
	int n_dim = (ni-1)*nj+1;			// Total dimension
	real sd1 = 1/sqrt(2*a);
	real sd_b = 1/sqrt(2*b);
}

generated quantities {
	matrix[N, n_dim] X;
	
	for(k in 1:N) {
		int idx = 1;
		
		// Root variable
		X[k, idx] = normal_rng(mu, sd1);
		real x1 = X[k, idx];
		idx += 1;
		
		// Blocks
		for(j in 1:nj) {
			real x_prev = x1;			// Reset x_prev for each block
			
			// Variables in block
			for(i in 2:ni) {
				real x_new = normal_rng(square(x_prev), sd_b);
				X[k, idx] = x_new;
				x_prev = x_new;
				idx += 1;
			}
		}
	}
}
