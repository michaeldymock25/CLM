gen_data <- function(n, arms, probs, T_rec, T_end){
  if(!is.character(arms)) arms <- as.character(arms)
  out <- data.table(t_rec = runif(n = n, min = 0, max = T_rec),
                    arm = factor(sample(x = arms, size = n, replace = TRUE), levels = arms))
  out[, y := rbinom(n = n, size = 1, prob = probs[arm])]
  out[, t_obs := ifelse(y == 1, runif(n = n, min = 0, max = T_end), NA)]
  return(out)
}

ext_dat <- function(data, T_follow){
  data[, t_obs := replace_na(t_obs, Inf)]
  new_cols <- paste0("end_", T_follow)
  for(t in 1:length(T_follow)) data[, (new_cols[t]) := as.integer(t_obs < T_follow[t])]
  return(data)
}

agg_dat <- function(data, T_follow, T_int){
  out <- lapply(2:length(T_follow), function(t){
                  temp <- data[t_rec < T_int - T_follow[t]]
                  for(tt in 1:(t-1)) temp <- temp[get(paste0("end_", T_follow[tt])) == 0]
                  temp[, .(n = .N, y = sum(get(paste0("end_", T_follow[t])))), keyby = arm]})
  return(out)
}

inc_dat <- function(data, T_follow, T_int){
  out <- lapply(2:(length(T_follow)-1), function(t){
                  temp <- data[t_rec < T_int - T_follow[t] & t_rec > T_int - T_follow[t+1]]
                  for(tt in 1:(t-1)) temp <- temp[get(paste0("end_", T_follow[tt])) == 0]
                  temp[levels(arm), .N, by = .EACHI]})
  out <- Reduce(merge, out)[, -"arm"]
  return(out)
}

trans_agg_dat <- function(data, T_follow, T_int){
  out <- lapply(2:(length(T_follow)-1), function(t){
                  temp <- data[t_rec < T_int - T_follow[length(T_follow)]]
                  for(tt in 1:(t-1)) temp <- temp[get(paste0("end_", T_follow[tt])) == 0]
                  temp[, .(n_star = .N, y_star = sum(get(paste0("end_", T_follow[length(T_follow)])))), keyby = arm]})
  return(out)
}

summ_fun <- function(draws) c(mean = mean(draws), median = median(draws), sd = sd(draws), quantile(draws, c(0.05, 0.95)))

clm <- function(data, T_follow, T_int, prior_sd, plot_it = TRUE, chains = 8, warmup_iter = 500, sampling_iter = 2000, n_cores = 1){

  J <- length(levels(data$arm))
  stan_mod <- cmdstan_model("CLM.stan")
  T_follow <- c(0, T_follow)
  data <- ext_dat(data, T_follow)
  
  data_agg <- agg_dat(data, T_follow, T_int)
  mod_data <- list(J = J,
                   `T` = length(T_follow) - 1, 
                   n = sapply(data_agg, function(x) x$n), 
                   y = sapply(data_agg, function(x) x$y), 
                   prior_sd = prior_sd)
  fit <- stan_mod$sample(
            data = mod_data,
            chains = chains,
            parallel_chains = min(chains, n_cores),
            refresh = 0,
            iter_warmup = warmup_iter, 
            iter_sampling = sampling_iter)
  
  beta_draws <- data.table(apply(fit$draws("beta"), 3, function(x) x))
  colnames(beta_draws) <- paste("beta", as.vector(outer(1:(length(T_follow)-1), 1:J, FUN = "paste", sep = "_")), sep = "_")
  
  pi_tau_draws <- beta_draws[, lapply(.SD, plogis)]
  colnames(pi_tau_draws) <- paste("pi", as.vector(outer(1:(length(T_follow)-1), 1:J, FUN = "paste", sep = "_")), sep = "_")
  
  pi_draws <- data.table(sapply(1:J, function(arm){
                          pis <- as.matrix(pi_tau_draws[,matrix(1:ncol(pi_tau_draws), ncol = J)[,arm], with = FALSE])
                          pis_temp <- cbind(1, 1, 1 - pis)
                          rowSums(sapply(1:ncol(pis), function(i) pis[,i]*rowProds(pis_temp[,1:(i+1)])))}))
  colnames(pi_draws) <- paste("pi", 1:J, sep = "_")

  return(list(beta_draws = beta_draws, pi_tau_draws = pi_tau_draws, pi_draws = pi_draws))
}  

discard <- function(data, T_follow, T_int, prior_sd, plot_it = TRUE, chains = 8, warmup_iter = 500, sampling_iter = 2000, n_cores = 1){
  
  J <- length(levels(data$arm))
  stan_mod <- cmdstan_model("discard.stan")
  T_follow <- c(0, T_follow)
  data <- ext_dat(data, T_follow)
  
  data <- data[t_enrol + T_follow[length(T_follow)] < T_int]
  mod_data <- list(J = J, 
                   n = as.vector(table(data$arm)),
                   y = unlist(data[, sum(y), keyby = arm][, -"arm"]),
                   prior_sd = prior_sd)
  fit <- stan_mod$sample(
            data = mod_data,
            chains = chains,
            parallel_chains = min(chains, n_cores),
            refresh = 0,
            iter_warmup = warmup_iter, 
            iter_sampling = sampling_iter)
  
  beta_draws <- data.table(apply(fit$draws("beta"), 3, function(x) x))
  colnames(beta_draws) <- paste("beta", 1:J, sep = "_")
  
  pi_draws <- beta_draws[, lapply(.SD, plogis)]
  colnames(pi_draws) <- paste("pi", 1:J, sep = "_")
  
  return(list(beta_draws = beta_draws, pi_draws = pi_draws))
}

transition <- function(data, T_follow, T_int, prior_sd, plot_it = TRUE, chains = 8, rep_data_sets = 10, warmup_iter = 500, sampling_iter = 2000, n_cores = 1){
  
  J <- length(levels(data$arm))
  stan_mod <- cmdstan_model("transition.stan")
  T_follow <- c(0, T_follow)
  data <- ext_dat(data, T_follow)
  
  setkey(data, arm)
  n_inc <- inc_dat(data, T_follow, T_int)
  out_inc <- trans_agg_dat(data, T_follow, T_int)
  n_star <- sapply(out_inc, function(x) x[,n_star])
  y_star <- sapply(out_inc, function(x) x[,y_star])
  mod_data <- list(J = J, 
                   `T` = length(T_follow) - 1,
                   n_inc = as.matrix(n_inc), 
                   n_star = as.matrix(n_star), 
                   y_star = as.matrix(y_star),
                   prior_sd = prior_sd)
  fit <- stan_mod$sample(
            data = mod_data,
            chains = max(chains, rep_data_sets),
            parallel_chains = min(chains, n_cores),
            refresh = 0,
            iter_warmup = max(warmup_iter, 1000), 
            iter_sampling = 1)
  y_inc_sum_draws <- fit$draws("y_inc_sum")
  y_inc_sum_draws <- apply(y_inc_sum_draws, 3, function(x) x)
  temp_data <- data[t_enrol + T_follow[length(T_follow)] < T_int]
  discard_mod <- cmdstan_model("discard.stan")
  beta_draws <- rbindlist(lapply(1:nrow(y_inc_sum_draws), function(i){
                            mod_data <- list(J = J,
                                             n = as.vector(table(temp_data$arm) + rowSums(n_inc)),
                                             y = unlist(temp_data[levels(arm), sum(y), by = .EACHI][,-"arm"]) + y_inc_sum_draws[i,],
                                             prior_sd = prior_sd)
                            fit <- discard_mod$sample(
                                                      data = mod_data,
                                                      chains = chains,
                                                      parallel_chains = min(chains, n_cores),
                                                      refresh = 0,
                                                      iter_warmup = warmup_iter, 
                                                      iter_sampling = sampling_iter)
                            beta_draws <- data.table(apply(fit$draws("beta"), 3, function(x) x))}))
  colnames(beta_draws) <- paste("beta", 1:J, sep = "_")
  
  pi_draws <- beta_draws[, lapply(.SD, plogis)]
  colnames(pi_draws) <- paste("pi", 1:J, sep = "_")
  
  return(list(beta_draws = beta_draws, pi_draws = pi_draws))
}

summarise <- function(draws, clm = TRUE){
  tab_pi <- t(apply(draws$pi_draws, 2, summ_fun))
  names(dimnames(tab_pi)) <- c("beta_time_arm", "statistic")
  tab_beta <- t(apply(draws$beta_draws, 2, summ_fun))
  names(dimnames(tab_beta)) <- c("beta_time_arm", "statistic")
  if(clm){
    tab_pi_tau <- t(apply(draws$pi_tau_draws, 2, summ_fun))
    names(dimnames(tab_pi_tau)) <- c("Arm", "statistic")
    return(list(tab_pi = tab_pi, tab_pi_tau = tab_pi_tau, tab_beta = tab_beta))
  } else {return(list(tab_pi = tab_pi, tab_beta = tab_beta))}
}

plot_pis <- function(draws, type = "arm"){
  if(type == "arm"){
    dat_vis_arm <- data.frame(Arm = fct_inorder(rep(paste("Arm", 1:ncol(draws$pi_draws)), each = nrow(draws$pi_draws))),
                              x = unlist(draws$pi_draws))
    p_arm <- ggplot(dat_vis_arm, aes(x = x)) + 
                geom_freqpoly(bins = 25, aes(y = stat(count / sum(count)), colour = Arm), size = 1) +
                scale_colour_manual("Arm", values = c("#000000", "#E69F00", "#56B4E9", "#009E73",
                                                      "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), drop = TRUE) +
                xlab(expression(pi[j])) +
                ylab("Density")
    return(p_arm)
  } else if(type == "all"){
    if(!("pi_tau_draws" %in% names(draws))) stop("this plot is only available for CLM output")
    dat_vis_all <- data.frame(Arm = fct_inorder(rep(paste("Arm", 1:ncol(draws$pi_draws)), 
                                                    each = nrow(draws$pi_tau_draws)*ncol(draws$pi_tau_draws)/ncol(draws$pi_draws))),
                              x = unlist(draws$pi_tau_draws),
                              Time = fct_inorder(rep(rep(paste(c(0, T_follow)[-length(T_follow)], "-", c(0, T_follow)[-1]), 
                                                         each = nrow(draws$pi_tau_draws)), 
                                                     ncol(draws$pi_draws))))
    p_all <- ggplot(dat_vis_all, aes(x = x)) + 
      facet_wrap(~Arm) +
      geom_freqpoly(bins = 25, aes(y = stat(count / sum(count)),  colour = Time), size = 1) +
      scale_colour_manual("Time Period", values = c("#000000", "#E69F00", "#56B4E9", "#009E73",
                                                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
      xlab(expression(pi[paste(j,tau)])) +
      ylab("Density")
    return(p_all)
  } else {
    stop("type must be either 'arm' or 'all'")
  }
}
