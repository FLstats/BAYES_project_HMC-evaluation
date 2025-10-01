data {
	real mu;
	real<lower=0> a;
	real<lower=0> b;
}

parameters {
	real z1;
	real z2;
}

transformed parameters {
	real x1 = (z1 + mu) * sqrt(1/(2*a));
	real x2 = (z2 + square(x1)) * sqrt(1/(2*b));
}

model {
	target += normal_lpdf(z1 | 0, 1);
	target += normal_lpdf(z2 | 0, 1);
}
