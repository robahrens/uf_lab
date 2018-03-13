// STAN model: Conn.scenario.1
data {
  int<lower=0> n; // number of indices
  int<lower=0> y; // number of years
  real U[n,y]; // observed abundance values
  real samp_errors[n,y]; // observed sample errors for abundance values
}
transformed data {
  real log_U[n,y]; 

  log_U = log(U); //transform to log space
}
parameters {
  real log_mu[y]; 
  real log_q[n]; 
  real<lower=0> proc_errors[n]; // process errors for each index
}
model {
  real sigma[n,y];
  real MU[n,y];

  // Priors
  log_mu ~ normal(log(100), 1);
  log_q ~ normal(log(0.01), 0.5);
  proc_errors ~ normal(0, 1);

  // Likelihood
  for (i in 1:n)
      for (j in 1:y){
        sigma[i,j] = (samp_errors[i,j] + proc_errors[i])^0.5;
        MU[i,j] = log_mu[j]+log_q[i];
        
        target += normal_lpdf(log_U[i,j] | MU[i,j], sigma[i,j]);
      }      
}
generated quantities {
  real mu[y];
  real q[n]; 

  mu = exp(log_mu); //transform
  q = exp(log_q); //transform
}
