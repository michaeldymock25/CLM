# CLM R package

This R package contains the following:

- Models:
  - Conditional model
  - Simple logistic model
  - Markov-transition model
- Functions:
  - odds(): Computes odds from probability
  - inv_odds(): Computes probability from odds
  - gen_data(): Generates simulated data in an appropriate format for analysis
  - agg_conditional(): Helper function for analyses using the conditional model (aggregation of data)
  - conditional_analysis(): Conducts an analysis using the conditional model
  - logistic_analysis(): Conducts an analysis using the logistic model
  - transition_analysis(): Conducts an analysis using the transition model
  - run_trial(): Runs a trial using data simulated by gen_data()
  - simulate_trials(): Simulates multiple trials using run_trial()

## Installation

```r
devtools::install_github("michaeldymock25/CLM")
```

## Compilation of STAN models

Unfortunately, the STAN models are not pre-compiled due to issues with package installation (with no immediate or obvious solution). Therefore, each time you load the package you will need to compile the models yourself using the below:

```r
logistic_mod <- cmdstanr::cmdstan_model(write_stan_file(readLines(url("https://raw.githubusercontent.com/michaeldymock25/CLM/main/inst/stan/logistic.stan"))))
conditional_mod <- cmdstanr::cmdstan_model(write_stan_file(readLines(url("https://raw.githubusercontent.com/michaeldymock25/CLM/main/inst/stan/conditional.stan"))))
```

Alternatively, the STAN models will be compiled during the first run of either the conditional_analysis(), logistic_analysis() or transition_analysis() functions. The first run of these functions will therefore be slower. I apologise for any inconvenience.
