## odds()
# requires p = probability to transform
# returns probability transformed to odds (numeric)

odds <- function(p) p/(1-p)

## inv_odds()
# requires odds = odds to transform
# returns odds transformed to probability (numeric)

inv_odds <- function(odds) odds/(1+odds)

## gen_data()
# requires n = trial sample size
#          J = number of trial arms
#          p = vector of event probabilities
#          recruit_period = number of days the trial will be recruiting
#          endpoint_time = number of days between randomisation and the time the endpoint is collected
#          follow_up_times = vector of follow up times in days (must include endpoint_time as the final follow up)
# randomises to arms with equal allocation probabilities
# generates ordered recruitment times in days
# generates event times in days (temporarily assumes all participants reach the endpoint)
# generates endpoint collection times
# generates events according to event probabilities p
# sets all event times for participants not reaching the endpoint to infinity
# generates new columns containing indicator for the endpoint being reached prior to each follow up time
# returns dataset (data.table) in a format ready for analysis

gen_data <- function(n, J, p, recruit_period, endpoint_time, follow_up_times){
  if(length(p) != J) stop("Ensure that length(p) = J (i.e, there is one event probability per arm)")
  if(round(last(follow_up_times)) != round(endpoint_time)) stop("Ensure the final follow up time is equal to endpoint_time")
  dat <- data.table(arm = sample(x = 1:J, size = n, replace = TRUE),
                    t_recruit = sort(round(runif(n = n, min = 0, max = recruit_period) + 0.5)))
  dat[, `:=`(t_event = t_recruit + round(runif(n, min = 0, max = endpoint_time) + 0.5),
             t_end = t_recruit + round(endpoint_time + 0.5),
             event = rbinom(n = n, size = 1, prob = p[arm])
  )]
  dat[, t_event := ifelse(event == 1, t_event, Inf)]
  new_cols <- paste("day", follow_up_times)
  for(t in 1:length(follow_up_times)) dat[, (new_cols[t]) := as.integer(t_event <= t_recruit + follow_up_times[t])]
  return(dat)
}

## agg_conditional()
# requires J = number of trial arms
#          dat = dataset in the format produced by gen_data()
#          follow_up_times = vector of follow up times in days (must include endpoint_time as the final follow up)
#          analysis_time = real time the analysis is conducted (must be later than endpoint_time)
# for each arm and follow up time computes the sample size and number of events conditional on the endpoint not yet being reached
# returns n (matrix) and y (matrix) ready for conditional_analysis()

agg_conditional <- function(J, dat, follow_up_times, analysis_time){
  if(analysis_time < last(follow_up_times)) stop("analysis_time must be later than the final follow up time")
  out <- lapply(1:length(follow_up_times), function(t){
           tmp <- dat[t_recruit < analysis_time - follow_up_times[t]]
           if(nrow(tmp) == 0) return(data.table(arm = 1:J, n = 0, y = 0))
           if(t > 1) for(tt in 1:(t-1)) tmp <- tmp[get(paste("day", follow_up_times[tt])) == 0]
           if(nrow(tmp) == 0) return(data.table(arm = 1:J, n = 0, y = 0))
           tmp[, .(n = .N, y = sum(get(paste("day", follow_up_times[t])))), keyby = arm]})
  return(list(n = sapply(out, function(x) x$n), y = sapply(out, function(x) x$y)))
}

## conditional_analysis()
# requires J = number of trial arms
#          follow_up_times = vector of follow up times in days (must include endpoint_time as the final follow up)
#          n = matrix with sample sizes by arms (rows) and follow up times (columns)
#              after first, subsequent columns must exclude participants with event times prior to the follow up times (see agg_conditional())
#          y = matrix with events by arms (rows) and follow up times (columns)
#              after first, subsequent columns must exclude participants with event times prior to the follow up times (see agg_conditional())
#          prior_mean = mean for the prior distribution on beta
#          prior_sd = standard deviation for the prior distribution on beta
#          ... = additional optional parameters for modelling
# loads conditional model if required
# sets up data required for conditional model
# samples from conditional model
# transforms beta parameter estimates to pi_t parameter estimates
# transforms pi_t parameter estimates to pi parameter estimates
# returns beta parameter estimates (data.table), pi_t parameter estimates (data.table) and pi parameter estimates (data.table)

conditional_analysis <- function(J, follow_up_times, n, y, prior_mean, prior_sd, ...){
  if(!exists("conditional_mod")) conditional_mod <- cmdstan_model(stan_file = system.file("stan", "conditonal.stan", package = "CLM"))
  mod_data <- list(J = J,
                   `T` = length(follow_up_times),
                   n = n,
                   y = y,
                   prior_mean = prior_mean,
                   prior_sd = prior_sd)
  drp <- utils::capture.output(fit <- conditional_mod$sample(data = mod_data, refresh = 0, ...))

  beta_draws <- data.table(posterior::as_draws_matrix(fit$draws("beta")))
  colnames(beta_draws) <- paste("beta", as.vector(outer(1:length(follow_up_times), 1:J, FUN = "paste", sep = "_")), sep = "_")
  pi_t_draws <- beta_draws[, lapply(.SD, plogis)]
  colnames(pi_t_draws) <- paste("pi", as.vector(outer(1:length(follow_up_times), 1:J, FUN = "paste", sep = "_")), sep = "_")
  pi_draws <- data.table(sapply(1:J, function(arm){
    pis <- as.matrix(pi_t_draws[,matrix(1:ncol(pi_t_draws), ncol = J)[,arm], with = FALSE])
    pis_temp <- cbind(1, 1, 1 - pis)
    rowSums(sapply(1:ncol(pis), function(t) pis[,t]*matrixStats::rowProds(pis_temp[,1:(t+1)])))}))
  colnames(pi_draws) <- paste("pi", 1:J, sep = "_")

  beta_draws <- data.table::melt(beta_draws,
                                 measure.vars = paste("beta", as.vector(outer(1:length(follow_up_times), 1:J, FUN = "paste",
                                                                              sep = "_")), sep = "_"),
                                 value.name = "sample")
  pi_t_draws <- data.table::melt(pi_t_draws,
                                 measure.vars = paste("pi", as.vector(outer(1:length(follow_up_times), 1:J, FUN = "paste",
                                                                            sep = "_")), sep = "_"),
                                 value.name = "sample")
  pi_draws <- data.table::melt(pi_draws, measure.vars = paste("pi", 1:J, sep = "_"), value.name = "sample")
  return(list(beta_draws = beta_draws, pi_t_draws = pi_t_draws, pi_draws = pi_draws))
}

## logistic_analysis()
# requires J = number of trial arms
#          dat = dataset in the format produced by gen_data() with at least "arm" and "event" columns
#          prior_mean = mean for the prior distribution on beta
#          prior_sd = standard deviation for the prior distribution on beta
#          ... = additional optional parameters for modelling
# loads logistic model if required
# sets up data required for logistic model
# samples from logistic model
# transforms beta parameter estimates to pi parameter estimates
# returns beta parameter estimates (data.table) and pi parameter estimates (data.table)

logistic_analysis <- function(J, dat, prior_mean, prior_sd, ...){
  if(!exists("logistic_mod")) logistic_mod <- cmdstan_model(stan_file = system.file("stan", "logistic.stan", package = "CLM"))
  mod_data <- list(J = J,
                   n = as.vector(dat[, table(arm)]),
                   y = as.vector(dat[event == 1, table(arm)]),
                   prior_mean = prior_mean,
                   prior_sd = prior_sd)
  drp <- utils::capture.output(fit <- logistic_mod$sample(data = mod_data, refresh = 0, ...))

  beta_draws <- data.table(posterior::as_draws_matrix(fit$draws("beta")))
  colnames(beta_draws) <- paste("beta", 1:J, sep = "_")
  pi_draws <- beta_draws[, lapply(.SD, plogis)]
  colnames(pi_draws) <- paste("pi", 1:J, sep = "_")

  beta_draws <- data.table::melt(beta_draws, measure.vars = paste("beta", 1:J, sep = "_"), value.name = "sample")
  pi_draws <- data.table::melt(pi_draws, measure.vars = paste("pi", 1:J, sep = "_"), value.name = "sample")
  return(list(beta_draws = beta_draws, pi_draws = pi_draws))
}

## transition_analysis()
# requires J = number of trial arms
#          dat = dataset in the format produced by gen_data() with at least "arm" follow up time point columns (with NA for yet to be observed)
#          follow_up_times = vector of follow up times in days (must include endpoint_time as the final follow up)
#          analysis_time = real time the analysis is conducted (must be later than endpoint_time)
#          n = matrix with sample sizes by arms (rows) and follow up times (columns)
#              after first, subsequent columns must exclude participants with event times prior to the follow up times (see agg_conditional())
#          y = matrix with events by arms (rows) and follow up times (columns)
#              after first, subsequent columns must exclude participants with event times prior to the follow up times (see agg_conditional())
#          t_q = last follow up time point for each participant
#          prior_mean = mean for the prior distribution on beta
#          prior_sd = standard deviation for the prior distribution on beta
#          nsets = if using the transition model, the number of predicted data sets to generate (defaults to 10)
#          ... = additional optional parameters for modelling
# loads conditional model and logistic model if required
# sets up data required for conditional model
# samples from conditional model (equal to the number of sets to be produced)
# transforms beta parameter estimates to pi_t parameter estimates
# generates nsets predicted data sets using recursive formula
# for each data set computes sample size and number of events per arm
# runs logistic model for each set of data and combines posterior distributions
# returns beta parameter estimates (data.table), pi_t parameter estimates (data.table) and pi parameter estimates (data.table)

transition_analysis <- function(J, dat, follow_up_times, analysis_time, n, y, t_q, prior_mean, prior_sd, nsets = 10, ...){
  if(!exists("conditional_mod")) conditional_mod <- cmdstan_model(stan_file = system.file("stan", "conditonal.stan", package = "CLM"))
  if(!exists("logistic_mod")) logistic_mod <- cmdstan_model(stan_file = system.file("stan", "logistic.stan", package = "CLM"))
  mod_data <- list(J = J,
                   `T` = length(follow_up_times),
                   n = n,
                   y = y,
                   prior_mean = prior_mean,
                   prior_sd = prior_sd)
  drp <- utils::capture.output(fit <- conditional_mod$sample(data = mod_data, refresh = 0, chains = nsets, iter_sampling = 1))

  beta_draws <- data.table(posterior::as_draws_matrix(fit$draws("beta")))
  colnames(beta_draws) <- paste("beta", as.vector(outer(1:length(follow_up_times), 1:J, FUN = "paste", sep = "_")), sep = "_")
  pi_t_draws <- beta_draws[, lapply(.SD, plogis)]
  colnames(pi_t_draws) <- paste("pi", as.vector(outer(1:length(follow_up_times), 1:J, FUN = "paste", sep = "_")), sep = "_")

  pred_list <- list()
  for(k in 1:nsets){
    dat_pred <- copy(dat)
    for(t in sort(unique(t_q))){
      if(t != length(follow_up_times)){
        for(tt in t:(length(follow_up_times) - 1)){
          pred <- rbinom(n = sum(t_q == t),
                         size = 1,
                         prob = unlist(pi_t_draws[k, paste("pi", tt + 1, dat_pred[t_q == t]$arm, sep = "_"), with = FALSE]))
          if(tt > 0) pred <- ifelse(dat_pred[t_q == t, get(paste("day", follow_up_times[tt]))] == 1, 1, pred)
          dat_pred[t_q == t, paste("day", follow_up_times[tt + 1]) := pred]
        }
      }
    }
    pred_list[[k]] <- dat_pred[, .(n = .N, y = sum(get(paste("day", last(follow_up_times))))), keyby = arm]
  }
  pred <- rbindlist(pred_list, idcol = "set")

  beta_draws <- rbindlist(lapply(1:nsets, function(i){
                            mod_data <- list(J = J,
                                             n = pred[set == i]$n,
                                             y = pred[set == i]$y,
                                             prior_mean = prior_mean,
                                             prior_sd = prior_sd)
                            drp <- utils::capture.output(fit <- logistic_mod$sample(data = mod_data, refresh = 0, ...))
                            data.table(posterior::as_draws_matrix(fit$draws("beta")))}))
  colnames(beta_draws) <- paste("beta", 1:J, sep = "_")
  pi_draws <- beta_draws[, lapply(.SD, plogis), .SDcols = paste("beta", 1:J, sep = "_")]
  colnames(pi_draws) <- paste("pi", 1:J, sep = "_")

  beta_draws <- data.table::melt(beta_draws, measure.vars = paste("beta", 1:J, sep = "_"), value.name = "sample")
  pi_draws <- data.table::melt(pi_draws, measure.vars = paste("pi", 1:J, sep = "_"), value.name = "sample")
  return(list(beta_draws = beta_draws, pi_draws = pi_draws))
}

## run_trial()
# requires J = number of trial arms
#          dat = dataset in the format produced by gen_data()
#          follow_up_times = vector of follow up times in days (must include endpoint_time as the final follow up)
#          analysis_times = vector of analysis times in days (first must be later than endpoint_time)
#          model = one of "conditional", "logistic" or "transition" (defaults to "conditional")
#          prior_mean = mean for the prior distribution on beta
#          prior_sd = standard deviation for the prior distribution on beta
#          nsets = if using the transition model, the number of predicted data sets to generate (defaults to 10)
#          ... = additional optional parameters for modelling
# for each analysis_time runs an analysis with the provided model
# returns beta parameter estimates (data.table) and pi parameter estimates (data.table)

run_trial <- function(J, dat, follow_up_times, analysis_times, model = "conditional", prior_mean, prior_sd, nsets = 10, ...){
  if(min(analysis_times) < last(follow_up_times)) stop("First analysis time must be later than the final follow up time")
  pi_draws_list <- beta_draws_list <- list()
  for(t in 1:length(analysis_times)){
    cat("Running analysis", t, "\n")
    if(model == "conditional"){
      tmp <- agg_conditional(J = J, dat = dat, follow_up_times = follow_up_times, analysis_time = analysis_times[t])
      out <- conditional_analysis(J = J, follow_up_times = follow_up_times, n = tmp$n, y = tmp$y,
                                  prior_mean = prior_mean, prior_sd = prior_sd, ...)
    }
    if(model == "logistic"){
      out <- logistic_analysis(J = J, dat = dat[t_end <= analysis_times[t]], prior_mean = prior_mean, prior_sd = prior_sd, ...)
    }
    if(model == "transition"){
      tmp <- agg_conditional(J = J, dat = dat, follow_up_times = follow_up_times, analysis_time = analysis_times[t])
      t_q <- dat[, lapply(follow_up_times, function(tt) analysis_times[t] - t_recruit > tt)][, .(t = sum(.SD)), by = seq_len(nrow(dat))]$t
      dat_tmp <- dat[t_recruit <= analysis_times[t], c("arm", paste("day", follow_up_times)), with = FALSE]
      for(tt in 1:length(follow_up_times)) dat_tmp[t_q < tt, paste("day", follow_up_times[tt]) := NA]
      out <- transition_analysis(J = J, dat = dat_tmp, follow_up_times = follow_up_times, analysis_time = analysis_times[t],
                                 n = tmp$n, y = tmp$y, t_q = t_q, prior_mean = prior_mean, prior_sd = prior_sd, nsets = nsets, ...)
    }
    beta_draws_list[[t]] <- out$beta_draws
    pi_draws_list[[t]] <- out$pi_draws
  }
  beta_draws <- rbindlist(beta_draws_list, idcol = "analysis")
  pi_draws <- rbindlist(pi_draws_list, idcol = "analysis")
  return(list(beta_draws = beta_draws, pi_draws = pi_draws))
}

## simulate_trials()
# requires nsim = number of trials to simulate
#          n = trial sample size
#          J = number of trial arms
#          p = vector of event probabilities
#          recruit_period = number of days the trial will be recruiting
#          endpoint_time = number of days between randomisation and the time the endpoint is collected#          J = number of trial arms
#          follow_up_times = vector of follow up times in days (must include endpoint_time as the final follow up)
#          analysis_times = vector of analysis times in days (first must be later than endpoint_time)
#          prior_mean = mean for the prior distribution on beta
#          prior_sd = standard deviation for the prior distribution on beta
#          nsets = if using the transition model, the number of predicted data sets to generate (defaults to 10)
#          ... = additional optional parameters for modelling
# generates nsim data sets
# sets up data.table to contain run times for each model
# for each data set runs a trial using each model as specified
# returns beta parameter estimates (data.table) and pi parameter estimates (data.table)

simulate_trials <- function(nsim, n, J, p, recruit_period, endpoint_time, follow_up_times, analysis_times,
                            prior_mean, prior_sd, nsets = 10, ...){
  models <- c("conditional", "logistic", "transition")
  dat <- lapply(1:nsim, function(i) gen_data(n = n, J = J, p = p, recruit_period = recruit_period,
                                             endpoint_time = endpoint_time, follow_up_times = follow_up_times))
  run_time <- data.table(model = models, time = NA_real_)
  pi_draws_list <- beta_draws_list <- list()
  for(mod in c("conditional", "logistic", "transition")){
    start_time <- Sys.time()
    out <- parallel::mclapply(dat, function(d)
                         run_trial(J = J, dat = d, follow_up_times = follow_up_times, analysis_times = analysis_times,
                                   model = mod, prior_mean = prior_mean, prior_sd = prior_sd, nsets = nsets, ...),
                              mc.cores = num_cores)
    end_time <- Sys.time()
    run_times[model == mod, time := end_time - start_time]
    beta_draws_list[[which(mod == models)]] <- rbindlist(lapply(out, function(x) x$beta_draws), idcol = "sim")
    pi_draws_list[[which(mod == models)]] <- rbindlist(lapply(out, function(x) x$pi_draws), idcol = "sim")
  }

  beta <- rbindlist(beta_draws_list, idcol = "model")
  pi <- rbindlist(pi_draws_list, idcol = "model")
  return(list(beta = beta, pi = pi, run_time = run_time))
}
