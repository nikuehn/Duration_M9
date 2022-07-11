/* **************************************************************************
Subduction Duration Model

uses functional form of Bahrampouri et al., (EQS 2021)
without Z1 term

 ************************************************************************** */

#include functions.stan

data {
  int<lower=1> N;                   // number of records
  int<lower=1> NEQ;                 // number of events
   int<lower=0> NTEST;

  vector[N] R;                      // distance
  vector[NTEST] R_test;
  vector[N] Y;                      // log Duration values
  vector[NTEST] Y_test;                      // log Duration values

  array[N] int<lower=1,upper=NEQ> eq;     // event index
}

transformed data {
  real delta_r = 1;
  //real<lower=0> R_b = 25;
  vector[N] Y_exp = exp(Y);
  vector[NTEST] Y_test_exp = exp(Y_test);
  int K = 1000;

}

parameters {
  real<lower=0> rp_0;
  real<lower=-rp_0> rp_1;
  real<lower=0> rp_s;

  real<lower=0> r1;
  real<lower=0> r2;
  real<lower=0> R_b;
  real<lower=0> R_s;

  real mu_ln_c1;
  real<lower=0> sigma_ln_c1;

  vector[NEQ] ln_c1;
}

transformed parameters {
  vector[NEQ] c1 = exp(ln_c1);
}

model {
  vector[N] sigma_path = rp_0 + rp_1 * inv_logit(rp_s * (R - R_s));
  rp_0 ~ normal(0, 1);
  rp_1 ~ normal(0, 1);
  rp_s ~ normal(1, 2);

  // priors for 
  r1 ~ normal(0,0.2);
  r2 ~ normal(0,0.2);
  R_b ~ normal(60, 40);
  R_s ~ normal(60, 40);

  mu_ln_c1 ~ normal(3.2, 2);
  sigma_ln_c1 ~ normal(0,2);

  ln_c1 ~ normal(mu_ln_c1, sigma_ln_c1);

  vector[N] mu;
  for(i in 1:N) {
    mu[i] = log(logistic_hinge(R[i], R_b, c1[eq[i]] + r1 * R_b, r1, r2, delta_r));
  }
  Y ~ normal(mu, sigma_path);
}

generated quantities{
  vector[NTEST] log_lik_test;

   //Likelihood
  for(i in 1:NTEST) {
    real sigma_path = rp_0 + rp_1 * inv_logit(rp_s * (R_test[i] - R_s));
    
    vector[K] ll;
    for(j in 1:K) {
      real c1_sam = lognormal_rng(mu_ln_c1, sigma_ln_c1);
      real mu = log(logistic_hinge(R_test[i], R_b, c1_sam + r1 * R_b, r1, r2, delta_r));
      ll[j] =  lognormal_lpdf(Y_test_exp[i]| mu, sigma_path);
    }

    log_lik_test[i] = log_sum_exp(ll) - log(K);
  }
   
}

