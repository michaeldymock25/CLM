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
