/* (In the likelihood)
*  to_vector(X) flattens the matrix X to a vector of length N x P.
*  [ X[1,1], X[1,2], ..., X[1,P],..., X[2,1],..., X[N,P] ].
*  The likelihood then evaluates all measurements (n,p) on the
*  normal(0, sqrt(exp(v))).
*/

data {
  int N;  // Number of observations
  int P;  // Number x-variables
  matrix[N, P] X;
}

parameters {
  real v;
}

model {
  // Log likelihood
  target += normal_lpdf(to_vector(X) | 0, sqrt(exp(v)));

  // Log priors
  target += normal_lpdf(v | 0, 3);
}
