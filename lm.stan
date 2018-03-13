// STAN model: Conn.scenario.1
data {
  int nobs; // number of data points
  real x[nobs]; // observed abundance values
  real y[nobs]; // observed sample errors for abundance values
}

parameters {
  real m; 
  real b; 
  real<lower=0> sigma; // process errors for each index
}

model {
  real ypred[nobs];
  // Priors
  m ~ normal(.2, 1);
  b ~ normal(4, 3);
  sigma ~ normal(0, 1);

  // Likelihood
  for (i in 1:nobs){
        ypred[i]=b+m*x[i];
        target += normal_lpdf(y[i] | ypred[i], sigma);
  }
}
//generated quantities {
 
// }
