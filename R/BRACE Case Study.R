## Demonstration of the conditional model for the BRACE study

### use simulation to compare the trial operating characteristics under the original model (simple logistic),
### the conditional model and the transition model

## Study Overview

### phase III, two group, multicentre, randomised placebo-controlled trial to determine if BCG vaccination
### reduces the incidence and severity of COVID-19 disease among healthcare workers
### endpoint is severe COVID-19 disease within 6 months post randomisation
### two arms are control = flu vaccine; intervention = BCG vaccine
### assumed 4% of control group will reach the endpoint within 6 months
### powered at 90% with 5% type one error with a trial sample size of 8,062 to detect one third reduction
### assumed constant recruitment over 6 months and endpoint available immediately if reached
### single interim analysis with stopping rule for superiority occurs once 100 endpoints reached
### simple logistic model

## Additional assumptions

### follow up participants every two months up to six months endpoint
### weakly informative prior mean and standard deviation on the log odds scale because we know the probabilities will be approximately <10%

## Load packages

#devtools::install_github("michaeldymock25/CLM") ## if package has not yet been installed
library(CLM)
library(cmdstanr)
library(data.table)
library(ggplot2)

## Set global parameters

rm(list = ls())                              ## remove objects in environment
set.seed(65757)                              ## set seed to ensure reproducible results
nsim <- 10000                                ## number of trial simulations
num_cores <- 12                              ## number of cores to run in parallel
J <- 2                                       ## number of trial arms
recruit_period <- 365.25/2 	                 ## length of recruitment period in days
endpoint_time <- 365.25/2                    ## time endpoint is measured
follow_up_times <- round(365.25/6*c(1:3))    ## follow up every two months
n <- 8062 	                                 ## maximum trial sample size
events_at_interim <- 100                     ## events required for interim analysis
p0 <- 0.04 			                             ## probability for control arm
p1 <- 2/3*0.04                               ## probability for intervention arm
p_null <- c(p0, p0)                          ## vector of event probabilities for null scenario
p_power <- c(p0, p1) 		                     ## vector of event probabilities for powered scenario
prior_means <- c(conditional = qlogis(p0)/3, ## centre prior mean on control arm probability/3 as there are three follow up periods
                 logistic = qlogis(p0),      ## centre prior mean on control arm probability
                 transition = qlogis(p0))    ## centre prior mean on control arm probability
prior_sds <- c(conditional = 1,              ## prior standard deviations of 1 on log odds scale (weakly informative)
               logistic = 1,
               transition = 1)
thresholds <- seq(0.90, 0.995, by = 0.005)   ## range of decision thresholds for probability of superiority
nsets <- 10                                  ## number of (predicted) data sets to compute for the transition model
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#117733", "#332288", "#AA4499", "#44AA99",
                             "#999933", "#882255", "#661100", "#6699CC", "#888888")

## read in STAN models

conditional_mod <- cmdstan_model(write_stan_file(readLines(url(
                     "https://raw.githubusercontent.com/michaeldymock25/CLM/main/inst/stan/conditional.stan"
                   ))))

logistic_mod <- cmdstan_model(write_stan_file(readLines(url(
                  "https://raw.githubusercontent.com/michaeldymock25/CLM/main/inst/stan/logistic.stan"
                ))))

## need to define new run_trial() function as the timing of the interim analysis is dependent not on time but on the accrued events

run_trial <- function(J, dat, p_true, events_at_interim, follow_up_times, prior_means, prior_sds, thresholds, nsets){
  t_interim <- ifelse(dat[, sum(event)] >= events_at_interim,
                      dat[min(which(cumsum(event) == events_at_interim)), t_end],
                      dat[, max(t_end)])

  interim_logistic <- logistic_analysis(J = J, dat = dat[t_end <= t_interim], prior_mean = prior_means["logistic"],
                                        prior_sd = prior_sds["logistic"])
  final_logistic <- logistic_analysis(J = J, dat = dat, prior_mean = prior_means["logistic"], prior_sd = prior_sds["logistic"])

  tmp_interim <- agg_conditional(J = J, dat = dat, follow_up_times = follow_up_times, analysis_time = t_interim)
  tmp_final <- agg_conditional(J = J, dat = dat, follow_up_times = follow_up_times, analysis_time = max(dat$t_end))

  interim_conditional <- conditional_analysis(J = J, follow_up_times = follow_up_times, n = tmp_interim$n, y = tmp_interim$y,
                                              prior_mean = prior_means["conditional"], prior_sd = prior_sds["conditional"])
  final_conditional <- conditional_analysis(J = J, follow_up_times = follow_up_times, n = tmp_final$n, y = tmp_final$y,
                                            prior_mean = prior_means["conditional"], prior_sd = prior_sds["conditional"])

  t_q_interim <- dat[, lapply(follow_up_times, function(tt) t_interim - t_recruit > tt)][, .(t = sum(.SD)), by = seq_len(nrow(dat))]$t
  t_q_final <- dat[, lapply(follow_up_times, function(tt) max(dat$t_end) - t_recruit > tt)][, .(t = sum(.SD)), by = seq_len(nrow(dat))]$t

  dat_tmp_interim <- dat[t_recruit <= t_interim, c("arm", paste("day", follow_up_times)), with = FALSE]
  dat_tmp_final <- dat[t_recruit <= max(dat$t_end), c("arm", paste("day", follow_up_times)), with = FALSE]

  for(t in 1:length(follow_up_times)){
    dat_tmp_interim[t_q_interim < t, paste("day", follow_up_times[t]) := NA]
    dat_tmp_final[t_q_final < t, paste("day", follow_up_times[t]) := NA]
  }

  interim_transition <- transition_analysis(J = J, dat = dat_tmp_interim, follow_up_times = follow_up_times,
                                            analysis_time = t_interim, n = tmp_interim$n, y = tmp_interim$y, t_q = t_q_interim,
                                            prior_mean = prior_means["transition"], prior_sd = prior_sds["transition"], nsets = nsets)
  final_transition <- transition_analysis(J = J, dat = dat_tmp_final, follow_up_times = follow_up_times, analysis_time = max(dat$t_end),
                                          n = tmp_final$n, y = tmp_final$y, t_q = t_q_final, prior_mean = prior_means["transition"],
                                          prior_sd = prior_sds["transition"], nsets = nsets)

  pi_draws <- rbindlist(list(rbindlist(list(interim_logistic$pi_draws, final_logistic$pi_draws), idcol = "analysis"),
                             rbindlist(list(interim_conditional$pi_draws, final_conditional$pi_draws), idcol = "analysis"),
                             rbindlist(list(interim_transition$pi_draws, final_transition$pi_draws), idcol = "analysis")),
                        idcol = "model")
  rmse <- rbindlist(lapply(1:3, function(mod)
                    rbindlist(lapply(1:2, function(anlys)
                                 RMSE(p_true = p_true, pi_draws = pi_draws[model == mod & analysis == anlys, -c("model", "analysis")])),
                              idcol = "analysis")),
                    idcol = "model")
  rmse[, `:=`(analysis = factor(analysis, labels = c("interim", "final")),
              model = factor(model, labels = c("logistic", "conditional", "transition")),
              variable = factor(variable, labels = c("pi_1", "pi_2"))
  )]
  supr <- rbindlist(lapply(1:3, function(mod)
                    rbindlist(lapply(1:2, function(anlys)
                                 superiority(pi_draws = pi_draws[model == mod & analysis == anlys, -c("model", "analysis")],
                                             thresholds = thresholds, base_var = "pi_1", comp_var = "pi_2", dir = "lesser")),
                              idcol = "analysis")),
                    idcol = "model")
  supr_tmp <- supr[, .(supr = any(supr)), by = .(model, threshold)]
  supr_tmp[, analysis := 3]
  supr <- dplyr::bind_rows(supr, supr_tmp)
  supr[, `:=`(analysis = factor(analysis, labels = c("interim", "final", "experiment-wise")),
              model = factor(model, labels = c("logistic", "conditional", "transition")),
              threshold = factor(threshold, labels = thresholds)
  )]
  gc()
  return(list(rmse = rmse, supr = supr))
}

## need to define new simulate_trials() function

simulate_trials <- function(nsim, n, J, p, recruit_period, endpoint_time, follow_up_times, prior_means,
                            prior_sds, thresholds, nsets, num_cores){
  dat <- parallel::mclapply(1:nsim, function(i) gen_data(n = n, J = J, p = p, recruit_period = recruit_period,
                                                         endpoint_time = endpoint_time, follow_up_times = follow_up_times),
                            mc.cores = num_cores)
  out <- parallel::mclapply(dat, function(d) run_trial(J = J, dat = d, p_true = p, follow_up_times = follow_up_times,
                                                       events_at_interim = events_at_interim, prior_means = prior_means,
                                                       prior_sds = prior_sds, thresholds = thresholds, nsets = nsets),
                            mc.cores = num_cores)
  rmse <- rbindlist(lapply(out, function(x) x$rmse), idcol = "sim")
  supr <- rbindlist(lapply(out, function(x) x$supr), idcol = "sim")
  return(list(rmse = rmse, supr = supr))
}

## run and save simulations

res_null <- simulate_trials(nsim = nsim, n = n, J = J, p = p_null, recruit_period = recruit_period, endpoint_time = endpoint_time,
                            follow_up_times = follow_up_times, prior_means = prior_means, prior_sds = prior_sds, thresholds = thresholds,
                            nsets = nsets, num_cores = num_cores)

saveRDS(res_null, "sims_null.rds")

res_power <- simulate_trials(nsim = nsim, n = n, J = J, p = p_power, recruit_period = recruit_period, endpoint_time = endpoint_time,
                             follow_up_times = follow_up_times, prior_means = prior_means, prior_sds = prior_sds, thresholds = thresholds,
                             nsets = nsets, num_cores = num_cores)

saveRDS(res_power, "sims_power.rds")

## summarise results

res_null <- readRDS("sims_null.rds")
res_power <- readRDS("sims_power.rds")

### type 1 error and power

supr <- dplyr::bind_rows(res_null$supr, res_power$supr, .id = "scenario")
supr <- supr[, .(supr = mean(supr)),by = .(scenario, model, analysis, threshold)]
supr[, `:=`(
            scenario = factor(scenario, labels = c("Null Scenario", "Powered Scenario")),
            analysis = factor(analysis, labels = c("Interim Analysis", "Final Analysis", "Experiment-Wise")),
            threshold = as.numeric(as.character(threshold))
)]

ggplot(supr, aes(x = threshold, y = supr, colour = model)) +
  facet_grid(scenario ~ analysis, scales = "free_y") +
  geom_line() +
  geom_hline(data = data.table(scenario = c("Null Scenario", "Powered Scenario"), yint = c(0.05, 0.80)),
             mapping = aes(yintercept = yint), linetype = "dashed") +
  scale_y_continuous("Probability of Superiority", n.breaks = 7) +
  scale_x_continuous("Decision Threshold", breaks = seq(0.90, 0.99, 0.01)) +
  scale_colour_manual(name = "Model",
                      limits = c("conditional", "logistic", "transition"),
                      values = safe_colorblind_palette) +
  theme(legend.position = "bottom")

### RMSE

rmse <- dplyr::bind_rows(res_null$rmse, res_power$rmse, .id = "scenario")
rmse[, `:=`(
            scenario = factor(scenario, labels = c("Null Scenario", "Powered Scenario")),
            analysis = factor(analysis, labels = c("Interim Analysis", "Final Analysis"))
)]

ggplot(rmse, aes(x = RMSE, colour = model, linetype = variable)) +
  facet_grid(scenario ~ analysis, scales = "free") +
  geom_density() +
  ylab("Density") +
  scale_colour_manual(name = "Model",
                      limits = c("conditional", "logistic", "transition"),
                      values = safe_colorblind_palette) +
  scale_linetype(name = "Parameter", labels = c(latex2exp::TeX(" $\\pi_1$"), latex2exp::TeX(" $\\pi_2$"))) +
  theme(legend.position = "bottom")
