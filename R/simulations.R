## Simulations to compare the trial operating characteristics and performance of the conditional model against
## the simple logistic model (that ignores participants with incomplete data) and the Markov-transition model that
## uses the posterior predictive distribution to estimate the unobserved data and then uses the simple logistic model

## Load packages

#devtools::install_github("michaeldymock25/CLM") ## if package has not yet been installed
library(CLM)
library(cmdstanr)
library(data.table)

## Set global parameters

rm(list = ls())             ## remove objects in environment
set.seed(65757)             ## set seed to ensure reproducible results

## read in STAN models

conditional_mod <- cmdstan_model(write_stan_file(readLines(url(
                     "https://raw.githubusercontent.com/michaeldymock25/CLM/main/inst/stan/conditional.stan"
                   ))))

logistic_mod <- cmdstan_model(write_stan_file(readLines(url(
                  "https://raw.githubusercontent.com/michaeldymock25/CLM/main/inst/stan/logistic.stan"
                ))))

## define configurations to simulate

cfg <- CJ(n = list(1000, 10000),             ## trial sample sizes of 1000 and 10000
          pi_1 = list(0.05, 0.20, 0.50),     ## control arm probabilities of 0.05, 0.20 and 0.50
          OR = list(1, 1.2, 1.5),            ## odds ratios of 1, 1.2 and 1.5
          sorted = FALSE)

## run the simulations for each configuration

for(z in 1:nrow(cfg)){
  pi_1 <- unlist(cfg[z][["pi_1"]])                                          ## extract control arm probability
  pi_2 <- inv_odds(unlist(cfg[z][["OR"]])*odds(unlist(cfg[z][["pi_1"]])))   ## compute treatment arm probability
  prior_means <- c(conditional = qlogis(pi_1/4),                            ## prior mean is 1/4 control arm probability for 4 follow up periods
                   logistic = qlogis(pi_1),                                 ## prior mean is control arm probability
                   transition = qlogis(pi_1))                               ## prior mean is control arm probability
  prior_sds <- c(conditional = 1,                                           ## prior standard deviations of 1 (weakly informative)
                 logistic = 1,
                 transition = 1)
  res <- simulate_trials(nsim = 10000,                                      ## number of simulations to run
                         n = unlist(cfg[z][["n"]]),                         ## trial sample size
                         J = 2,                                             ## number of trial arms
                         p = c(pi_1, pi_2),                                 ## true arm probabilities
                         recruit_period = 365.25,                           ## 12 months recruitment period
                         endpoint_time = 365.25/3,                          ## endpoint collected at 4 months
                         follow_up_times = 365.25/3*seq(1/4, 1, by = 1/4),  ## follow up times monthly until 4 months
                         analysis_times = 365.25*seq(2/3, 4/3, by = 1/3),   ## scheduled analyses at 8 months, 12 months and 16 months
                         prior_means = prior_means,
                         prior_sds = prior_sds,
                         thresholds = seq(0.90, 0.995, by = 0.005),         ## range of decision thresholds for probability of superiority
                         num_cores = 20,                                    ## number of cores to run in parallel
                         simplify_output = TRUE)                            ## do not need to store full posterior distributions
  saveRDS(list(cfg = cfg[z], res = res), paste0("sims_cfg_", formatC(z, width = 2, flag = "0"), ".rds"))
}
