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
# sets up data required for conditional model
# samples from conditional model
# transforms beta parameter estimates to pi_t parameter estimates
# transforms pi_t parameter estimates to pi parameter estimates
# returns beta parameter estimates (data.table), pi_t parameter estimates (data.table) and pi parameter estimates (data.table)

conditional_analysis <- function(J, follow_up_times, n, y, prior_mean, prior_sd, ...){
  if(!exists("conditional_mod")) conditional_mod <- cmdstan_model(write_stan_file(readLines(url(
                                                       "https://raw.githubusercontent.com/michaeldymock25/CLM/main/inst/stan/conditional.stan"))))
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
# sets up data required for logistic model
# samples from logistic model
# transforms beta parameter estimates to pi parameter estimates
# returns beta parameter estimates (data.table) and pi parameter estimates (data.table)

logistic_analysis <- function(J, dat, prior_mean, prior_sd, ...){
  if(!exists("logistic_mod")) logistic_mod <- cmdstan_model(write_stan_file(readLines(url(
                                                 "https://raw.githubusercontent.com/michaeldymock25/CLM/main/inst/stan/logistic.stan"))))
  mod_data <- list(J = J,
                   n = dat[, .N, keyby = arm]$N,
                   y = dat[, .(y = sum(event)), keyby = arm]$y,
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
#          nsets = if using the transition model, the number of predicted data sets to generate
#          ... = additional optional parameters for modelling
# sets up data required for conditional model
# samples from conditional model (equal to the number of sets to be produced)
# transforms beta parameter estimates to pi_t parameter estimates
# generates nsets predicted data sets using recursive formula
# for each data set computes sample size and number of events per arm
# runs logistic model for each set of data and combines posterior distributions
# returns beta parameter estimates (data.table), pi_t parameter estimates (data.table) and pi parameter estimates (data.table)

transition_analysis <- function(J, dat, follow_up_times, analysis_time, n, y, t_q, prior_mean, prior_sd, nsets, ...){
  if(!exists("conditional_mod")) conditional_mod <- cmdstan_model(write_stan_file(readLines(url(
                                                       "https://raw.githubusercontent.com/michaeldymock25/CLM/main/inst/stan/conditional.stan"))))
  if(!exists("logistic_mod")) logistic_mod <- cmdstan_model(write_stan_file(readLines(url(
                                                 "https://raw.githubusercontent.com/michaeldymock25/CLM/main/inst/stan/logistic.stan"))))
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
                                
## RMSE()
# requires p_true = true probabilities                                
#          pi_draws = data.table of posterior draws for pi
# returns data.table with RMSE 
                                
RMSE <- function(p_true, pi_draws){
  rbindlist(lapply(1:length(p_true), function(j) pi_draws[variable == paste0("pi_", j), .(RMSE = sqrt(sum((sample - p_true[j])^2)))]), idcol = "variable")
}  
                   
## superiority()
# requires pi_draws = data.table of posterior draws for pi
#          thresholds = range of decision thresholds 
#          base_var = variable name for base of comparision
#          comp_var = variable name for to compare to base_var
#          dir = direction of comparison (defaults to greater)
# returns data.table with superiority decision                   
       
superiority <- function(pi_draws, thresholds, base_var, comp_var, dir = "greater"){
  OR_sample <- odds(pi_draws[variable == comp_var]$sample)/odds(pi_draws[variable == base_var]$sample)
  if(dir == "greater"){
    rbindlist(lapply(thresholds, function(thr) data.table(supr = mean(OR_sample > 1) >= thr)), idcol = "threshold")
  } else if(dir == "lesser"){
    rbindlist(lapply(thresholds, function(thr) data.table(supr = mean(OR_sample < 1) >= thr)), idcol = "threshold")
  } else {
    stop("dir must be 'greater' or 'lesser'")
  }
}
                    
## run_trial()
# requires J = number of trial arms
#          dat = dataset in the format produced by gen_data()
#          p_true = true probabilities                                
#          follow_up_times = vector of follow up times in days (must include endpoint_time as the final follow up)
#          analysis_times = vector of analysis times in days (first must be later than endpoint_time)
#          model = one of "conditional", "logistic" or "transition"
#          prior_mean = mean for the prior distribution on beta
#          prior_sd = standard deviation for the prior distribution on beta
#          thresholds = range of decision thresholds                    
#          nsets = if using the transition model, the number of predicted data sets to generate
#          base_var = variable name for base of comparision
#          comp_var = variable name for to compare to base_var                  
#          ... = additional optional parameters for modelling
# for each analysis_time runs an analysis with the provided model
# returns beta parameter estimates (data.table), pi parameter estimates (data.table), rmse (data.table) and supr (data.table)

run_trial <- function(J, dat, p_true, follow_up_times, analysis_times, model, prior_mean, prior_sd, thresholds, nsets, base_var, comp_var, ...){
  if(min(analysis_times) < last(follow_up_times)) stop("First analysis time must be later than the final follow up time")
  supr_list <- rmse_list <- pi_draws_list <- beta_draws_list <- list()
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
      t_q <- rowSums(dat[t_recruit <= analysis_times[t], lapply(follow_up_times, function(tt) analysis_times[t] - t_recruit >= tt)])
      dat_tmp <- dat[t_recruit <= analysis_times[t], c("arm", paste("day", follow_up_times)), with = FALSE]
      for(tt in 1:length(follow_up_times)) dat_tmp[t_q < tt, paste("day", follow_up_times[tt]) := NA]
      out <- transition_analysis(J = J, dat = dat_tmp, follow_up_times = follow_up_times, analysis_time = analysis_times[t],
                                 n = tmp$n, y = tmp$y, t_q = t_q, prior_mean = prior_mean, prior_sd = prior_sd, nsets = nsets, ...)
    }
    beta_draws_list[[t]] <- out$beta_draws
    pi_draws_list[[t]] <- out$pi_draws
    rmse_list[[t]] <- RMSE(p_true = p_true, pi_draws = out$pi_draws)
    supr_list[[t]] <- superiority(pi_draws = out$pi_draws, thresholds = thresholds, base_var = base_var, comp_var = comp_var)
  }
  beta_draws <- rbindlist(beta_draws_list, idcol = "analysis")
  pi_draws <- rbindlist(pi_draws_list, idcol = "analysis")
  rmse <- rbindlist(rmse_list, idcol = "analysis") 
  supr <- rbindlist(supr_list, idcol = "analysis")                                                              
  gc()
  return(list(beta_draws = beta_draws, pi_draws = pi_draws, rmse = rmse, supr = supr))
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
#          thresholds = range of decision thresholds
#          nsets = if using the transition model, the number of predicted data sets to generate (defaults to 10)             
#          base_var = variable name for base of comparision (defaults to "pi_1")
#          comp_var = variable name for to compare to base_var (defaults to "pi_2")
#          num_cores = number of cores to run simulations in parallel
#          simplify_output = if TRUE, does not return draws to reduce output size (defaults to FALSE)                                                               
#          ... = additional optional parameters for modelling
# generates nsim data sets
# sets up data.table to contain run times for each model
# for each data set runs a trial using each model as specified
# returns beta parameter estimates (data.table) and pi parameter estimates (data.table)

simulate_trials <- function(nsim, n, J, p, recruit_period, endpoint_time, follow_up_times, analysis_times, prior_mean, prior_sd, 
                            thresholds, nsets = 10, base_var = "pi_1", comp_var = "pi_2", num_cores, simplify_output = FALSE, ...){
  models <- c("conditional", "logistic", "transition")
  dat <- parallel::mclapply(1:nsim, function(i) gen_data(n = n, J = J, p = p, recruit_period = recruit_period,
                                                         endpoint_time = endpoint_time, follow_up_times = follow_up_times),
                            mc.cores = num_cores)
  run_time <- data.table(model = models, time = NA_real_)
  supr_list <- rmse_list <- pi_draws_list <- beta_draws_list <- list()
  for(mod in models){
    start_time <- Sys.time()
    out <- parallel::mclapply(dat, function(d)
                         run_trial(J = J, dat = d, p_true = p, follow_up_times = follow_up_times, analysis_times = analysis_times,
                                   model = mod, prior_mean = prior_mean, prior_sd = prior_sd, thresholds = thresholds, nsets = nsets,
                                   base_var = base_var, comp_var = comp_var, ...),
                              mc.cores = num_cores)
    end_time <- Sys.time()
    tt <- end_time - start_time
    units(tt) <- "mins"
    run_time[model == mod, time := tt]
    beta_draws_list[[which(mod == models)]] <- rbindlist(lapply(out, function(x) x$beta_draws), idcol = "sim")
    pi_draws_list[[which(mod == models)]] <- rbindlist(lapply(out, function(x) x$pi_draws), idcol = "sim")
    rmse_list[[which(mod == models)]] <- rbindlist(lapply(out, function(x) x$rmse), idcol = "sim")
    supr_list[[which(mod == models)]] <- rbindlist(lapply(out, function(x) x$supr), idcol = "sim")
  }

  beta <- rbindlist(beta_draws_list, idcol = "model")
  pi <- rbindlist(pi_draws_list, idcol = "model")
  rmse <- rbindlist(rmse_list, idcol = "model")
  supr <- rbindlist(supr_list, idcol = "model")
  if(simplify_output){
    return(list(run_time = run_time, rmse = rmse, supr = supr))
  } else {
    return(list(beta = beta, pi = pi, run_time = run_time, rmse = rmse, supr = supr))
  }  
}
