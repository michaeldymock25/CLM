data{
  int<lower=1> J;             // number of trial arms
  int<lower=2> T;             // number of follow up times 
  int<lower=0> n_inc[J,T-1];  // number of participants with incomplete observations
  int<lower=0> n_star[J,T-1]; // number of completed participants that had not observed the endpoint at time tau     
  int y_star[J,T-1];          // outcomes of completed participants that had not observed the endpoint at time tau   
  real prior_mean;            // prior mean for coefficients
  real prior_sd;              // prior standard deviation for coefficients
}

parameters{
  vector[J] beta_star[T-1];   // coefficients for transition probabilities
}

transformed parameters{
  vector[J] pis_star[T-1] = inv_logit(beta_star);
}

model{
  for(t in 1:(T-1)){
    y_star[,t] ~ binomial_logit(n_star[,t], beta_star[t]);
    beta_star[t] ~ normal(prior_mean, prior_sd);
  }
}

generated quantities{
  int y_inc[J,T-1];             // imputed incomplete observations
  int y_inc_sum[J];
  
  for(t in 1:(T-1)) y_inc[,t] = binomial_rng(n_inc[,t], pis_star[t]);
  for(j in 1:J) y_inc_sum[j] = sum(y_inc[j,]);
}
