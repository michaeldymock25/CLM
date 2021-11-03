data{
  int<lower=1> J;         // number of trial arms
  int<lower=1> n[J];      // number of individuals on each arm
  int y[J];               // outcomes on each arm
  real prior_sd;          // prior standard deviation for coefficients
}

parameters{
  vector[J] beta;         // coefficients for arms
}

model{
  y ~ binomial_logit(n, beta);
  beta ~ normal(0, prior_sd);
}
