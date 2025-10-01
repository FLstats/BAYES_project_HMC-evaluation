data {
	real mu;
	real<lower=0> a;
	real<lower=0> b;
}

parameters {
	real x1;
	real x2;
}

model {
	target += normal_lpdf(x1 | mu, sqrt( 1/(2*a) ) );
	target += normal_lpdf(x2 | square(x1), sqrt( 1/(2*b) ) );
}
