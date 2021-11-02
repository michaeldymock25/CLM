ext_dat <- function(data, T_follow){
  data$t_obs <- replace_na(data$t_obs, Inf)
  for(t in T_follow) data[,paste0("end_",t)] <- ifelse(data$t_obs < t, 1, 0)
  return(data)
}

agg_dat <- function(data, T_follow, T_int){
  out <- lapply(levels(data$arm), function(trt){
    out <- sapply(2:length(T_follow), function(t){
      temp <- data %>% filter(t_enrol < T_int - T_follow[t], arm == trt)
      for(tt in 1:(t-1)) temp <- temp %>% filter(get(paste0("end_", T_follow[tt])) == 0)
      c(n = nrow(temp), y = sum(temp %>% select(paste0("end_", T_follow[t]))))})
    colnames(out) <- paste("Follow up", T_follow[-length(T_follow)])
    out})
  names(out) <- paste("Arm", levels(data$arm))
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
                   n = t(sapply(data_agg, function(dat) dat["n",])),
                   y = t(sapply(data_agg, function(dat) dat["y",])),
                   prior_sd = prior_sd)
  fit <- stan_mod$sample(
      data = mod_data,
      chains = chains,
      parallel_chains = min(chains, n_cores),
      refresh = 0,
      iter_warmup = warmup_iter,
      iter_sampling = sampling_iter)
  beta_draws <- fit$draws("beta")
  tab_beta <- t(apply(beta_draws, 3, summ_fun))
  names(dimnames(tab_beta)) <- c("beta[time,arm]", "statistic")
  pis_draws <- plogis(apply(beta_draws, 3, function(x) x))
  pis <- sapply(1:length(levels(data$arm)), function(arm){
    pis <- pis_draws[,matrix(1:ncol(pis_draws), ncol = length(levels(data$arm)))[,arm]]
    pis_temp <- cbind(1, 1, 1 - pis)
    rowSums(sapply(1:ncol(pis), function(i) pis[,i]*rowProds(pis_temp[,1:(i+1)])))})
  
  tab_pi <- t(apply(pis, 2, summ_fun))
  dimnames(tab_pi)[[1]] <- levels(data$arm)
  names(dimnames(tab_pi)) <- c("Arm", "statistic")
  out <- list(tab_beta = tab_beta, tab_pi = tab_pi)

  if(plot_it){
    dat_vis_all <- data.frame(Arm = fct_inorder(rep(paste("Arm", levels(data$arm)), each = length(pis_draws)/length(levels(data$arm)))),
                              x = as.vector(pis_draws),
                              Time = fct_inorder(rep(rep(paste(T_follow[-length(T_follow)], "-", T_follow[-1]), each = nrow(pis_draws)), length(levels(data$arm)))))
    p_all <- ggplot(dat_vis_all, aes(x = x)) +
        facet_wrap(~Arm) +
        geom_freqpoly(bins = 25, aes(y = stat(count / sum(count)),  colour = Time), size = 1) +
        scale_colour_manual("Time Period", values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
        xlab(expression(pi[paste(j,tau)])) +
        ylab("Density")
      print(p_all)
    dat_vis_arm <- data.frame(Arm = fct_inorder(rep(paste("Arm", levels(data$arm)), each = nrow(pis))), x = as.vector(pis))
    p_arm <- ggplot(dat_vis_arm, aes(x = x)) +
      geom_freqpoly(bins = 25, aes(y = stat(count / sum(count)), colour = Arm), size = 1) +
      scale_colour_manual("Arm", values = c("#000000", "#E69F00", "#56B4E9", "#009E73",
                                            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
      xlab(expression(pi[j])) +
      ylab("Density")
    print(p_arm)
  }
  return(out)
}
