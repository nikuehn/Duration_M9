/* **************************************************************************
Subduction Duration Model

 ************************************************************************** */

#include functions.stan

data {
  int<lower=1> N;                   // number of records
  int<lower=1> NEQ;                 // number of events

  vector[NEQ] MEQ;
  vector[N] R;                      // distance
  vector[N] Y;                      // log Duration values

  array[N] int<lower=1,upper=NEQ> eq;     // event index
}

transformed data {
  real delta_r = 1;
  //real<lower=0> R_b = 25;
  vector[N] Y_exp = exp(Y);

  real ln_denom = log(3.2) + log(4.9) + log(10^6);
  vector[NEQ] ln_moment = log(10) * (1.5 * MEQ + 16.05);
  real ln_moment_M9 = log(10) * (1.5 * 9 + 16.05);
  real exp03 = -1./3;

  int K = 1000;
}

parameters {
  real<lower=0> rp_0;
  real<lower=0> rp_1;

  real<lower=0> r1;
  real<lower=0> r2;
  real<lower=0> R_b;

  real mu_ln_c1;
  real<lower=0> sigma_ln_c1;
}

generated quantities{
  vector[N] resid;

  //Likelihood
  for(i in 1:N) {
    real mu_ln_dsigma = (mu_ln_c1 + ln_denom) / exp03 + ln_moment_M9;
    //real sigma_ln_dsigma = sigma_ln_c1 * 3;
    //real ln_c1 = exp03 * (mu_ln_dsigma - ln_moment[eq[i]]) - ln_denom;
    //real mu = logistic_hinge(R[i], R_b, exp(ln_c1) + r1 * R_b, r1, r2, delta_r); 
    //resid_1[i] = Y[i] - log(mu);

    //vector[K] mub;
    //for(j in 1:K) {
    //  real c1b = lognormal_rng(ln_c1, sigma_ln_c1);
    //  mub[j] = logistic_hinge(R[i], R_b, c1b + r1 * R_b, r1, r2, delta_r); 
    //}
    //resid_2[i] = Y[i] - log(mean(mub));

    real c1c = exp(exp03 * (mu_ln_dsigma - ln_moment[eq[i]]) - ln_denom + 0.5 * square(sigma_ln_c1));
    real muc = logistic_hinge(R[i], R_b, c1c + r1 * R_b, r1, r2, delta_r);
    resid[i] = Y[i] - log(muc);
  }
   
}

