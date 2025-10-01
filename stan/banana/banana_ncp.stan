data {
	real b;									// b<0: concave banana
}

parameters {
	real z1;
	real z2;
}

transformed parameters {
	real x1 = z1;
	real x2 = z2 + b*(z1^2 - 100);			// twist transformation
}

model {
	target += normal_lpdf(z1 | 0, 10);		// sd = 10
	target += normal_lpdf(z2 | 0, 1);
}
