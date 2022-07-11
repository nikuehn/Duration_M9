/* **************************************************************************
Subduction Duration Model

uses functional form of Bahrampouri et al., (EQS 2021)
without Z1 term

 ************************************************************************** */

#include functions.stan

data {
  int<lower=1> N;                   // number of records
  int<lower=1> NEQ;                 // number of events

  vector[NEQ] MEQ;                  // magnitude
  vector[N] R;                      // distance
  vector[N] Y;                      // log Duration values

  array[N] int<lower=1,upper=NEQ> eq;     // event index

}

transformed data {
  real delta_r = 1;
  //real<lower=0> R_b = 25;
  vector[N] Y_exp = exp(Y);
  real R_ref = 300;

}

parameters {
  real<lower=0> rp_0;
  real<lower=-rp_0> rp_1;
  real<lower=0> rp_s;

  real<lower=0> r1;
  real<lower=0> r2;
  real<lower=0> R_b;
  real<lower=0> R_s;

  real<lower=0> c1;
  real<lower=0> sigma_source;

  vector[NEQ] eqterm;
}

model {
  vector[N] sigma_path = rp_0 + rp_1 * inv_logit(rp_s * (R - R_s));
  rp_0 ~ normal(0, 1);
  rp_1 ~ normal(0, 1);
  rp_s ~ normal(1, 2);
  sigma_source ~ normal(0, 1);

  // priors for 
  r1 ~ normal(0,0.2);
  r2 ~ normal(0,0.2);
  R_b ~ normal(60, 40);
  R_s ~ normal(60, 40);

  c1 ~ normal(25, 5);

  eqterm ~ normal(0, sigma_source);

  vector[N] mu;
  for(i in 1:N) {
    mu[i] = log(logistic_hinge(R[i], R_b, c1 + r1 * R_b, r1, r2, delta_r)) + eqterm[eq[i]];
  }
  Y ~ normal(mu, sigma_path);
}

generated quantities{
  vector[N] log_lik;
  vector[N] resid;
  vector<lower=0>[N] y_rep;

  //Likelihood
  for(i in 1:N) {
    real sigma_path = rp_0 + rp_1 * inv_logit(rp_s * (R[i] - R_s));
    real mu = log(logistic_hinge(R[i], R_b, c1 + r1 * R_b, r1, r2, delta_r)) + eqterm[eq[i]];

    log_lik[i] = lognormal_lpdf(Y_exp[i]| mu, sigma_path);
    resid[i] = Y[i] - mu;
    y_rep[i] = lognormal_rng(mu, sigma_path);
  }

}

