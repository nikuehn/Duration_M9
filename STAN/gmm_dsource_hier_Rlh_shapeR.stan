/* **************************************************************************
Subduction Duration Model

 ************************************************************************** */

#include functions.stan

data {
  int<lower=1> N;                   // number of records
  int<lower=1> NEQ;                 // number of events

  vector[N] R;                      // distance
  vector[N] Y;                      // log Duration values

  array[N] int<lower=1,upper=NEQ> eq;     // event index
}

transformed data {
  real delta_r = 1;
  //real<lower=0> R_b = 25;
  vector[N] Y_exp = exp(Y);

}

parameters {
  real<lower=0> rp_0;
  real<lower=0> rp_1;

  real<lower=0> r1;
  real<lower=0> r2;
  real<lower=0> R_b;

  real mu_ln_c1;
  real<lower=0> sigma_ln_c1;

  vector[NEQ] ln_c1;
}

transformed parameters {
  vector[NEQ] c1 = exp(ln_c1);
}

model {
  vector[N] shape_path = rp_0 + rp_1 * R;
  rp_0 ~ normal(0, 5);
  rp_1 ~ normal(0, 1);

  // priors for 
  r1 ~ normal(0,0.2);
  r2 ~ normal(0,0.2);
  R_b ~ normal(60, 40);

  mu_ln_c1 ~ normal(3.2, 2);
  sigma_ln_c1 ~ normal(0,2);

  ln_c1 ~ normal(mu_ln_c1, sigma_ln_c1);

  vector[N] mu;
  for(i in 1:N) {
    mu[i] = logistic_hinge(R[i], R_b, c1[eq[i]] + r1 * R_b, r1, r2, delta_r);
  }
  Y_exp ~ gamma(shape_path, shape_path ./ mu);
}

generated quantities{
  vector[N] log_lik;
  vector[N] resid;
  vector<lower=0>[N] y_rep;
  vector[NEQ] dln_c1 = ln_c1 - mu_ln_c1;

   //Likelihood
  for(i in 1:N) {
    real shape_path = rp_0 + rp_1 * R[i];
    real mu = logistic_hinge(R[i], R_b, c1[eq[i]] + r1 * R_b, r1, r2, delta_r); 

    log_lik[i] = gamma_lpdf(Y_exp[i]| shape_path, shape_path / mu);
    resid[i] = Y[i] - log(mu);
    y_rep[i] = gamma_rng(shape_path, shape_path / mu);
  }
   
}

