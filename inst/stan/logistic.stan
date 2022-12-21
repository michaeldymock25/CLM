data{
  int<lower=1> J;          // number of trial arms
  array<lower=1>[J] int n; // number of participants on each arm
  array<lower=1>[J] int y; // outcomes on each arm
  real prior_mean;         // prior mean for coefficients
  real prior_sd;           // prior standard deviation for coefficients
}

parameters{
  vector[J] beta;          // coefficients for arms
}

model{
  y ~ binomial_logit(n, beta);
  beta ~ normal(prior_mean, prior_sd);
}
