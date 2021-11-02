\name{CLM}
\alias{CLM}
\title{Fits Conditional Logistic Models for Adaptive Trials}
\usage{
CLM(data,
    T_follow,
    T_int,
    prior_sd,
    plot_it = TRUE,
    chains = 8,
    warmup_iter = 500,
    sampling_iter = 2000,
    n_cores = 1)
}
\arguments{
  \item{data}{a dataframe with at least the following columns:
             t_enrol - individual enrollment time
             t_obs - individual observation time (NA if endpoint not yet observed)
             arm - individual arm allocation
             y - individual binary endpoint}
  \item{T_follow}{vector of follow up time points (including the final)}
  \item{T_int}{interim analysis time point}
  \item{prior_sd}{standard deviation for the normal prior distribution}
  \item{plot_it}{if TRUE will plot posterior distributions}
  \item{chains}{number of chains to run}
  \item{warmup_iter}{number of warmup iterations}
  \item{sampling_iter}{number of sampling iterations}
  \item{n_cores}{number of cores if parallel processing is used}
}
\value{
  \item{tab_beta}{tabulated summary of the posterior distributions for the beta parameters}
  \item{tab_pi}{tabulated summary of the (transformed) posterior distributions for the pi parameters}
}
\description{
  This function fits a conditional logistic model to a potentially incomplete data set at a prespecified interim analysis for an adaptive trial. 
  The modelling specifications can be altered via the function's inputs. 
}
\examples{
CLM(data, T_follow = c(5, 10, 15, 20), T_int = 30, prior_sd = 5)
}
