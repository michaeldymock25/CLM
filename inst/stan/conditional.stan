data{
  int<lower=1> J;            // number of trial arms
  int<lower=1> T;            // number of follow up times
  array<lower=0>[J,T] int n; // number of participants on each arm at each follow up that have not yet reached the endpoint
  array<lower=0>[J,T] int y; // outcomes on each arm at each follow up for participants that have not yet reached the endpoint
  real prior_mean;           // prior mean for coefficients
  real prior_sd;             // prior standard deviation for coefficients
}

parameters{
  vector[J] beta[T];         // coefficients for follow up times and arms
}

model{
  for(t in 1:T){
    y[,t] ~ binomial_logit(n[,t], beta[t]);
    beta[t] ~ normal(prior_mean, prior_sd);
  }
}
