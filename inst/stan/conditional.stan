data{
  int<lower=1> J;            // number of trial arms
  int<lower=1> T;            // number of follow up times
  array[J,T] int<lower=0> n; // number of participants on each arm at each follow up that have not yet reached the endpoint
  array[J,T] int<lower=0> y; // outcomes on each arm at each follow up for participants that have not yet reached the endpoint
  real prior_mean;           // prior mean for coefficients
  real prior_sd;             // prior standard deviation for coefficients
}

parameters{
  array[T] vector[J] beta;   // coefficients for follow up times and arms
}

model{
  for(t in 1:T){
    y[,t] ~ binomial_logit(n[,t], beta[t]);
    beta[t] ~ normal(prior_mean, prior_sd);
  }
}
