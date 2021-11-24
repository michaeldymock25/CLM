
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
    if(nrow(temp) == 0) return(data.table(arm = 1:J, n = 0, y = 0))
    for(tt in 1:(t-1)) temp <- temp[get(paste0("end_", T_follow[tt])) == 0]
    temp[, .(n = .N, y = sum(get(paste0("end_", T_follow[t])))), keyby = arm]})
  return(out)
}

inc_dat <- function(data, J, T_follow, T_int){
  out <- lapply(2:(length(T_follow)-1), function(t){
    temp <- data[t_rec < T_int - T_follow[t] & t_rec > T_int - T_follow[t+1]]
    for(tt in 1:(t-1)) temp <- temp[get(paste0("end_", T_follow[tt])) == 0]
    temp[levels(arm), .N, by = .EACHI][,-"arm"]})
  return(matrix(unlist(out), nrow = J))
}

trans_agg_dat <- function(data, T_follow, T_int){
  out <- lapply(2:(length(T_follow)-1), function(t){
    temp <- data[t_rec < T_int - T_follow[length(T_follow)]]
    for(tt in 1:(t-1)) temp <- temp[get(paste0("end_", T_follow[tt])) == 0]
    temp[, .(n_star = .N, y_star = sum(get(paste0("end_", T_follow[length(T_follow)])))), keyby = arm]})
  return(out)
}

clm <- function(data, T_follow, T_int, prior_sd, plot_it = TRUE, chains = 8, warmup_iter = 500, sampling_iter = 2000, n_cores = 1){
  
  J <- length(levels(data$arm))
  stan_mod <- cmdstan_model("CLM.stan")
  T_follow <- c(0, T_follow)
  data <- ext_dat(data, T_follow)
  
  data_agg <- lapply(T_int, function(t_int) agg_dat(data, T_follow, T_int = t_int))
  names(data_agg) <- T_int
  
  beta_draws_list <- list()
  pi_t_draws_list <- list()
  pi_draws_list <- list()
  for(t in 1:length(T_int)){
    mod_data <- list(J = J,
                     `T` = length(T_follow) - 1, 
                     n = sapply(data_agg[[as.character(T_int[t])]], function(x) x$n), 
                     y = sapply(data_agg[[as.character(T_int[t])]], function(x) x$y), 
                     prior_sd = prior_sd)
    cat("Running interim", t, "\n")
    drp <- utils::capture.output(fit <- stan_mod$sample(
      data = mod_data,
      chains = chains,
      parallel_chains = min(chains, n_cores),
      refresh = 0,
      iter_warmup = warmup_iter, 
      iter_sampling = sampling_iter))
    
    beta_draws <- data.table(apply(fit$draws("beta"), 3, function(x) x))
    colnames(beta_draws) <- paste("beta", as.vector(outer(1:(length(T_follow)-1), 1:J, FUN = "paste", sep = "_")), sep = "_")
    pi_t_draws <- beta_draws[, lapply(.SD, plogis)]
    colnames(pi_t_draws) <- paste("pi", as.vector(outer(1:(length(T_follow)-1), 1:J, FUN = "paste", sep = "_")), sep = "_")
    pi_draws <- data.table(sapply(1:J, function(arm){
      pis <- as.matrix(pi_t_draws[,matrix(1:ncol(pi_t_draws), ncol = J)[,arm], with = FALSE])
      pis_temp <- cbind(1, 1, 1 - pis)
      rowSums(sapply(1:ncol(pis), function(i) pis[,i]*rowProds(pis_temp[,1:(i+1)])))}))
    colnames(pi_draws) <- paste("pi", 1:J, sep = "_")
    
    beta_draws_list[[t]] <- beta_draws
    pi_t_draws_list[[t]] <- pi_t_draws
    pi_draws_list[[t]] <- pi_draws
  }
  
  beta_draws <- rbindlist(beta_draws_list, idcol = "interim")
  pi_t_draws <- rbindlist(pi_t_draws_list, idcol = "interim")
  pi_draws <- rbindlist(pi_draws_list, idcol = "interim")
  
  return(list(beta_draws = beta_draws, pi_t_draws = pi_t_draws, pi_draws = pi_draws))
}  

discard <- function(data, T_follow, T_int, prior_sd, plot_it = TRUE, chains = 8, warmup_iter = 500, sampling_iter = 2000, n_cores = 1){
  
  J <- length(levels(data$arm))
  stan_mod <- cmdstan_model("discard.stan")
  T_follow <- c(0, T_follow)
  data <- ext_dat(data, T_follow)
  
  beta_draws_list <- list()
  pi_draws_list <- list()
  for(t in 1:length(T_int)){
    data_tmp <- data[t_rec + T_follow[length(T_follow)] < T_int[t]]
    mod_data <- list(J = J, 
                    n = as.vector(table(data_tmp$arm)),
                    y = unlist(data_tmp[, sum(y), keyby = arm][, -"arm"]),
                    prior_sd = prior_sd)
    cat("Running interim", t, "\n")
    drp <- utils::capture.output(fit <- stan_mod$sample(
      data = mod_data,
      chains = chains,
      parallel_chains = min(chains, n_cores),
      refresh = 0,
      iter_warmup = warmup_iter, 
      iter_sampling = sampling_iter))
  
    beta_draws <- data.table(apply(fit$draws("beta"), 3, function(x) x))
    colnames(beta_draws) <- paste("beta", 1:J, sep = "_")
    pi_draws <- beta_draws[, lapply(.SD, plogis)]
    colnames(pi_draws) <- paste("pi", 1:J, sep = "_")
    beta_draws_list[[t]] <- beta_draws
    pi_draws_list[[t]] <- pi_draws
  }
  
  beta_draws <- rbindlist(beta_draws_list, idcol = "interim")
  pi_draws <- rbindlist(pi_draws_list, idcol = "interim")
  
  return(list(beta_draws = beta_draws, pi_draws = pi_draws))
}

transition <- function(data, T_follow, T_int, prior_sd, plot_it = TRUE, chains = 8, rep_data_sets = 10, warmup_iter = 500, sampling_iter = 2000, n_cores = 1){
  
  J <- length(levels(data$arm))
  stan_mod <- cmdstan_model("transition.stan")
  T_follow <- c(0, T_follow)
  data <- ext_dat(data, T_follow)
  
  setkey(data, arm)
  
  n_inc <- lapply(T_int, function(t_int) inc_dat(data, J = J, T_follow, T_int = t_int))
  names(n_inc) <- T_int
  
  out_inc <- lapply(T_int, function(t_int) trans_agg_dat(data, T_follow, T_int = t_int))
  names(out_inc) <- T_int
  
  beta_draws_list <- list()
  pi_draws_list <- list()
  for(t in 1:length(T_int)){
    n_star <- sapply(out_inc[[as.character(T_int[t])]], function(x) x[,n_star])
    y_star <- sapply(out_inc[[as.character(T_int[t])]], function(x) x[,y_star])
    mod_data <- list(J = J, 
                     `T` = length(T_follow) - 1,
                     n_inc = n_inc[[t]], 
                     n_star = as.matrix(n_star), 
                     y_star = as.matrix(y_star),
                     prior_sd = prior_sd)
    cat("Running interim", t, "\n")
    drp <- utils::capture.output(fit <- stan_mod$sample(
      data = mod_data,
      chains = max(chains, rep_data_sets),
      parallel_chains = min(chains, n_cores),
      refresh = 0,
      iter_warmup = max(warmup_iter, 1000), 
      iter_sampling = 1))
    y_inc_sum_draws <- fit$draws("y_inc_sum")
    y_inc_sum_draws <- apply(y_inc_sum_draws, 3, function(x) x)
    data_tmp <- data[t_rec + T_follow[length(T_follow)] < T_int[t]]
    discard_mod <- cmdstan_model("discard.stan")
    beta_draws <- rbindlist(lapply(1:nrow(y_inc_sum_draws), function(i){
                            mod_data <- list(J = J,
                                             n = as.vector(table(data_tmp$arm) + rowSums(n_inc[[as.character(T_int[t])]])),
                                             y = unlist(data_tmp[levels(arm), sum(y), by = .EACHI][,-"arm"]) + y_inc_sum_draws[i,],
                                             prior_sd = prior_sd)
                            drp <- utils::capture.output(fit <- discard_mod$sample(
                                                      data = mod_data,
                                                      chains = chains,
                                                      parallel_chains = min(chains, n_cores),
                                                      refresh = 0,
                                                      iter_warmup = warmup_iter, 
                                                      iter_sampling = sampling_iter))
                            data.table(apply(fit$draws("beta"), 3, function(x) x))}), idcol = "rep")
    colnames(beta_draws) <- c("rep", paste("beta", 1:J, sep = "_"))
    pi_draws <- beta_draws[, lapply(.SD, plogis), .SDcols = paste("beta", 1:J, sep = "_")]
    colnames(pi_draws) <- paste("pi", 1:J, sep = "_")
    beta_draws_list[[t]] <- beta_draws
    pi_draws_list[[t]] <- pi_draws
  }
  
  beta_draws <- rbindlist(beta_draws_list, idcol = "interim")
  pi_draws <- rbindlist(pi_draws_list, idcol = "interim")                      
  
  return(list(beta_draws = beta_draws, pi_draws = pi_draws))
}

summ_fun <- function(draws) c(mean = mean(draws), median = median(draws), sd = sd(draws), quantile(draws, c(0.05, 0.95)))

summarise <- function(draws, model = "CLM"){
  J <- ncol(draws$pi_draws) - 1
  rem <- "interim"
  tab_pi <- rbindlist(lapply(1:length(unique(draws$pi_draws$interim)), function(int) 
                          data.table(par = paste("pi", 1:J, sep = "_"), t(apply(draws$pi_draws[interim == int, -..rem], 2, summ_fun)))), 
                      idcol = "interim")
  if(model == "transition") rem <- c("rep", rem)
  tab_beta <- rbindlist(lapply(1:length(unique(draws$beta_draws$interim)), function(int) 
                            data.table(par = paste("beta", 1:J, sep = "_"), t(apply(draws$beta_draws[interim == int, -..rem], 2, summ_fun)))), 
                        idcol = "interim")
  if(model == "CLM"){
    tab_pi_t <- rbindlist(lapply(1:length(unique(draws$pi_t_draws$interim)), function(int) 
                              data.table(par = paste("pi", as.vector(outer(1:((ncol(draws$pi_t_draws) - 1)/J), 1:J, FUN = "paste", sep = "_")), sep = "_"), 
                                         t(apply(draws$pi_t_draws[interim == int, -"interim"], 2, summ_fun)))), 
                          idcol = "interim")
    return(list(tab_pi = tab_pi, tab_pi_t = tab_pi_t, tab_beta = tab_beta))
  } else {return(list(tab_pi = tab_pi, tab_beta = tab_beta))}
}

plot_pis <- function(draws, T_follow, type = "arm"){
  if(type == "arm"){
    dat_vis_arm <- data.frame(Arm = fct_inorder(rep(paste("Arm", 1:(ncol(draws$pi_draws)-1)), each = nrow(draws$pi_draws))),
                              x = unlist(draws$pi_draws[,-"interim"]))
    p_arm <- ggplot(dat_vis_arm, aes(x = x)) + 
      geom_freqpoly(bins = 25, aes(y = stat(count / sum(count)), colour = Arm), size = 1) +
      scale_colour_manual("Arm", values = c("#000000", "#E69F00", "#56B4E9", "#009E73",
                                            "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), drop = TRUE) +
      xlab(expression(pi[j])) +
      ylab("Density")
    return(p_arm)
  } else if(type == "all"){
    if(!("pi_t_draws" %in% names(draws))) stop("this plot is only available for CLM output")
    dat_vis_all <- data.frame(Arm = fct_inorder(rep(paste("Arm", 1:(ncol(draws$pi_draws)-1)), 
                                                    each = nrow(draws$pi_t_draws)*(ncol(draws$pi_t_draws)-1)/(ncol(draws$pi_draws)-1))),
                              x = unlist(draws$pi_t_draws[,-"interim"]),
                              Time = fct_inorder(rep(rep(paste(c(0, T_follow)[-length(T_follow)], "-", c(0, T_follow)[-1]), 
                                                         each = nrow(draws$pi_t_draws)), 
                                                     (ncol(draws$pi_draws)-1))))
    p_all <- ggplot(dat_vis_all, aes(x = x)) + 
      facet_wrap(~Arm) +
      geom_freqpoly(bins = 25, aes(y = stat(count / sum(count)),  colour = Time), size = 1) +
      scale_colour_manual("Time Period", values = c("#000000", "#E69F00", "#56B4E9", "#009E73",
                                                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
      xlab(expression(pi[paste(j,t)])) +
      ylab("Density")
    return(p_all)
  } else {
    stop("type must be either 'arm' or 'all'")
  }
}


