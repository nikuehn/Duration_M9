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

  real ln_denom = log(3.2) + log(4.9) + log(10^6);
  vector[NEQ] ln_moment = log(10) * (1.5 * MEQ + 16.05);
  real ln_moment_M9 = log(10) * (1.5 * 9 + 16.05);
  real exp03 = -1./3;

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
}

generated quantities{
  vector[N] resid;

  //Likelihood
  for(i in 1:N) {
    real mu_ln_dsigma = (log(c1) + ln_denom) / exp03 + ln_moment_M9;
    real ln_c1 = exp03 * (mu_ln_dsigma - ln_moment[eq[i]]) - ln_denom;

    real mu = log(logistic_hinge(R[i], R_b, exp(ln_c1) + r1 * R_b, r1, r2, delta_r));

    resid[i] = Y[i] - mu;
  }

}

