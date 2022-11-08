data{
  int<lower=1> J;            // number of trial arms
  int<lower=1> T;            // number of follow up times
  int<lower=0> n[J,T];       // number of participants on each arm at each follow up 
  int y[J,T];                // outcomes on each arm at each follow up
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
