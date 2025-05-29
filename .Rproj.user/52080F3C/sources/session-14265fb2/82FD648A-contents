##############################################################################################
# this is the simulation code accompanying the manuscript
# "A novel approach to the design and sample size planning of animal experiments based on treatment effect estimation"
# by Dario Zocholl, Henrike Solveen, and Matthias Schmid

# authors of this simulation code: Dario Zocholl & Henrike Solveen
# corresponding author: Dario Zocholl, dario.zocholl@imbie.uni-bonn.de
# last modifications: May 29, 2025
##############################################################################################


# libraries
library(tidyverse)
library(doParallel)
library(foreach)
library(furrr)
library(rstan)
library(emmeans)
library(multcomp)
library(ggpubr)


##############################################################################################
# function to sample data based on the distribution and its parameters specified
sample_data <- function(distribution, parameters_new, no_groups, n_per_group){
  if(distribution == "normal"){
    sapply(1:no_groups, function(i) {
      rnorm(mean = parameters_new[[i]][1], sd = parameters_new[[i]][2], n = n_per_group)
    }) %>% as.data.frame()
  } else if(distribution == "lognormal"){
    sapply(1:no_groups, function(i) {
      rlnorm(meanlog = parameters_new[[i]][1], sdlog = parameters_new[[i]][2], n = n_per_group)
    }) %>% as.data.frame()
  } else if(distribution == "cauchy"){
    sapply(1:no_groups, function(i) {
      rcauchy(location = parameters_new[[i]][1], scale = parameters_new[[i]][2], n = n_per_group)
    }) %>% as.data.frame()
  }
}



##############################################################################################
# function to find the "best" pairwise comparison among the comparisons specified
find_best_comparison <- function(samples, type_of_comparison, custom_comparisons = NULL){
  
  if(type_of_comparison == "all_pairwise"){
    best_comparison <- samples %>% pivot_longer(., cols=everything(), names_to = "group") %>% 
      compare_means(value ~ group, data=., method = "t.test") %>% dplyr::select(group1, group2, p) %>% slice(which.min(p))
    best_groups <- as.numeric(gsub(pattern = "V", replacement = "", x = c(best_comparison$group1, best_comparison$group2)))
    results <- best_comparison
  } 
  
  else if(type_of_comparison == "all_vs_control"){
    glht_fit <- glht(aov(value ~ name, 
                         data = samples %>% pivot_longer(cols=everything()) %>% mutate(name = as.factor(name))), 
                     linfct = mcp(name = "Dunnett"))
    glht_fit_summary <- summary(glht_fit)
    # need to get the correct indices of the best comparison from glht-object
    best_comparison_name <- names(glht_fit_summary$test$tstat)[which.min(glht_fit_summary$test$pvalues)]
    best_groups <- as.numeric(gsub(pattern = "V", replacement = "", x = unlist(strsplit(best_comparison_name, split = "-"))))
    results <- glht_fit_summary
  } 
  
  else if(type_of_comparison == "custom"){# custom comparison
    # the trick is to:
    # 1. access the list of all possible contrasts,
    # 2. get the indices of the defined contrasts and transform into emmeans language,
    # 3. keep only the indices of the contrasts that were specified.
    custom_comparisons_emmeans <- gsub("(\\d+)", "V\\1", gsub(" vs\\. ", " - ", custom_comparisons))
    fit <- lm(value ~ name, 
              data = samples %>% pivot_longer(cols = everything()) %>% mutate(name = as.factor(name)))
    # if adjust = "none" is not specified, the Sidak method is used to adjust p-values
    all_contrasts <- summary(emmeans(fit, pairwise ~ name)$contrasts, adjust = "none") %>%
      filter(contrast %in% custom_comparisons_emmeans) 
    best_pvalue_dsg1 <- min(all_contrasts[,"p.value"])
    best_groups <- as.numeric(gsub(pattern = "V", replacement = "", x = unlist(strsplit(all_contrasts[which.min(all_contrasts[,"p.value"]),]$contrast, " - "))))
    results <- all_contrasts
  }
  
  return(list("best_groups"=best_groups, "results"=results))
}




##############################################################################################

estimate_effect <- function(data_stage1, data_stage2 = NULL, n_per_group, dsg, type_of_comparison, multiple_test = NULL, stanmodel, resampling_distributions){
  
  if((dsg == "standard" & is.null(multiple_test))){
    effect_est <- mean(data_stage1[,2] - data_stage1[,1]) * (-1)
    fit <- t.test(data_stage1[,1], data_stage1[,2])
    CI <- fit$conf.int[1:2]
    best_pvalue <- fit$p.value  
  }
  
  if((dsg == "repeat_and_pool" & is.null(multiple_test))){
    data <- rbind(data_stage1, data_stage2)
    effect_est <- mean(data[,2] - data[,1]) * (-1)
    fit <- t.test(data[,1], data[,2])
    CI <- fit$conf.int[1:2]
    best_pvalue <- fit$p.value  
  }
  
  if(dsg == "repeat_and_replace"){
    effect_est <- mean(data_stage2[,2] - data_stage2[,1]) * (-1)
    fit <- t.test(data_stage2[,1], data_stage2[,2])
    CI <- fit$conf.int[1:2]
    best_pvalue <- fit$p.value  
  }

  if(dsg == "bayes"){
    
    historical_data <- bind_rows(data.frame(y = data_stage1[,2], x = 0), data.frame(y = data_stage1[,1], x = 1))
    historical_fit <- lm(y ~ x, data = historical_data)
    
    fit <- rstan::sampling(stanmodel, 
                           data = list("df_beta" = historical_fit$df.residual, "y" = c(data_stage2[,2], data_stage2[,1]), 
                                       "group" = c(rep(0,length(data_stage2[,2])), rep(1,length(data_stage2[,1]))), "N" = length(c(data_stage2[,2], data_stage2[,1])), 
                                       "delta" = 0.5, "diff_prior" = as.numeric(historical_fit$coefficients[2]), "sigma_prior" = summary(historical_fit)$coefficients[2,2], 
                                       "diff_skeptical" = 0, "sigma_skeptical" = summary(historical_fit)$coefficients[2,2]), 
                           refresh = 0)      
    posterior_mixture <- rstan::extract(fit) %>% as_tibble()
    
    effect_est <- mean(posterior_mixture$beta)
    CI <- quantile(posterior_mixture$beta, probs = c(0.025, 0.975))
    best_pvalue <- min(sum(posterior_mixture$beta < 0)/nrow(posterior_mixture), sum(posterior_mixture$beta > 0)/nrow(posterior_mixture))*2
  }
  
  if(dsg == "resampled_shrinkage"){
    ttest <- t.test(data_stage1[,1], data_stage1[,2]) 
    backtransf_pval <- ecdf(resampling_distributions$pvalue_distr)(ttest$p.value)
    if(backtransf_pval == 0){backtransf_pval <- ecdf(resampling_distributions$pvalue_distr)(min(resampling_distributions$pvalue_distr))} 
    if(backtransf_pval == 1){backtransf_pval <- ecdf(resampling_distributions$pvalue_distr)(max(resampling_distributions$pvalue_distr))} 
    backtransf_est <- qt(p = backtransf_pval/2, df = ttest$parameter, lower.tail = FALSE) * ttest$stderr
    if(ttest$statistic > 0){backtransf_est <- (-1)*backtransf_est}
    ####
    effect_est <- backtransf_est
    CI <- c(backtransf_est + qt(p=0.025, df = ttest$parameter) * ttest$stderr, 
            backtransf_est + qt(p=0.975, df = ttest$parameter) * ttest$stderr)
    best_pvalue <- ecdf(resampling_distributions$pvalue_distr)(ttest$p.value)
  }
  
  return(list("effect_est" = effect_est, "CI_lwr" = CI[1], "CI_upr" = CI[2], "best_pvalue" = best_pvalue, "n_per_group" = n_per_group, "dsg" = dsg))
    
}




##############################################################################################

# obtain p-value distribution for resampled shrinkage estimator
obtain_resampling_distributions <- function(resample_n, type_of_comparison, distribution, no_groups, n_per_group, multiple_test = NULL){
  pvalue_distr <- c()
  est_distr <- c()
  for(i in 1:resample_n){
    samples <- sample_data(distribution = distribution, parameters_new = rep(list(c(0,1)), no_groups), 
                           no_groups = no_groups, n_per_group = n_per_group)
    best_comparison <- find_best_comparison(samples = samples, type_of_comparison = type_of_comparison)
    resampled_effect_estimates <- estimate_effect(data_stage1 = samples[, best_comparison$best_groups], dsg = "standard", type_of_comparison = type_of_comparison, n_per_group = n_per_group)
    pvalue_distr[i] <- resampled_effect_estimates$best_pvalue
    est_distr[i] <- resampled_effect_estimates$effect_est
  }
  return(list("pvalue_distr" = pvalue_distr, "est_distr" = est_distr))
}


##############################################################################################
# data sampling

get_data <- function(n_per_group, no_groups, distribution, parameters_new, type_of_comparison, custom_comparisons = NULL){
  
  ############################################################ 
  # data for single-stage designs
  samples <- sample_data(distribution = distribution, parameters_new = parameters_new, no_groups = no_groups, n_per_group = n_per_group) 
  
  # find best comparison
  results <- find_best_comparison(samples = samples, type_of_comparison = type_of_comparison)
  data_singlestage <- samples[, results$best_groups]
  
  
  ############################################################ 
  # data for two-stage designs
  n_per_group_stage1 <- floor((n_per_group * no_groups) / (no_groups+2)) 
  samples <- sample_data(distribution = distribution, parameters_new = parameters_new, no_groups = no_groups, n_per_group = n_per_group_stage1) 
  
  # find best comparison
  results <- find_best_comparison(samples = samples, type_of_comparison = type_of_comparison)
  data_stage1 <- samples[, results$best_groups]
  
  # the sample size for the second stage is the remaining sample size after first stage from the total sample size, which is n_per_group * no_groups
  n_per_group_stage2 = floor(((n_per_group * no_groups) - (n_per_group_stage1 * no_groups)) / 2)
  
  data_stage2 <- sample_data(distribution = distribution, 
                                 parameters_new = parameters_new[results$best_groups], 
                                 no_groups = 2, n_per_group = n_per_group_stage2) 
  colnames(data_stage2) <- colnames(samples[, results$best_groups])
  
  ############################################################ 
  return(list("data_singlestage" = data_singlestage, "data_stage1" = data_stage1, "data_stage2" = data_stage2))
}


##############################################################################################
########## Function wrapper for simulation ##########

run_simulation <- function(dsg, data_stage1, data_stage2, n_per_group, no_groups, distribution, type_of_comparison, stanmodel = NULL, resampling_distributions = NULL){
  
  # use pmap to iterate through rows of the grid
  data <- estimate_effect(
    data_stage1 = data_stage1, data_stage2 = data_stage2, n_per_group = n_per_group, dsg = dsg, type_of_comparison = type_of_comparison,
    multiple_test = NULL, stanmodel = stanmodel, resampling_distributions = resampling_distributions
  )
  return(data %>% bind_rows() %>% add_column("no_groups" = no_groups, "distribution" = distribution, "type_of_comparison" = type_of_comparison))
}

##############################################################################################
# RUN simulation
nsim <- 1e4

# define simulation scenarios
dsg <- c("standard", "repeat_and_pool", "repeat_and_replace", "bayes")
n_per_group <- c(6, 12, 18)
no_groups <- 2:10
distribution <- c("lognormal")
effect <- c(0,1,2)
type_of_comparison <- c("all_pairwise")

# create a grid containing all simulation scenarios
param_grid <- expand.grid(
  # dsg = dsg,
  n_per_group = n_per_group,
  no_groups = no_groups,
  distribution = distribution,
  type_of_comparison = type_of_comparison,
  effect = effect,
  stringsAsFactors = FALSE
)

# add sample sizes for two-stage designs
param_grid <- param_grid %>% 
  mutate(n_per_group_stage1 = floor((n_per_group * no_groups) / (no_groups+2))) %>% 
  rowwise() %>% mutate(n_per_group_stage2 = floor(((n_per_group * no_groups) - (n_per_group_stage1 * no_groups)) / 2))


##### preparations 
stanmodel <- stan_model("mixture_model_ttest.stan")

# define parameters for effect & number of groups 
generate_parameters <- function(effect, no_groups) {
  return(mapply(c, c(effect, rep(0, no_groups - 1)), rep(1, no_groups), SIMPLIFY = FALSE))
}

# add parameters to the grid 
param_grid$parameters_new <- mapply(generate_parameters, param_grid$effect, param_grid$no_groups)


# make a grid with all simulation runs
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...)) # expand.grid for data.frames
param_grid_sim <- expand.grid.df(param_grid, data.frame("simID" = 1:nsim)) %>% as_tibble()

param_grid_sim <- expand.grid.df(
  param_grid,
  expand.grid(simid = 1:nsim,
              dsg = dsg)
)


#####data##############################################################################
### effect estimation

(starttime <- Sys.time())
registerDoParallel(detectCores()-1)

data <- foreach(i = 1:nrow(param_grid), .combine = 'bind_rows', .packages = c("tidyverse", "rstan", "emmeans", "multcomp", "ggpubr")) %dopar% {

  sim_results <- list()
  for(sim in 1:nsim){
    generated_data <- get_data(n_per_group = param_grid$n_per_group[i], no_groups = param_grid$no_groups[i],
                               distribution = param_grid$distribution[i], parameters_new = param_grid$parameters_new[[i]],
                               type_of_comparison = param_grid$type_of_comparison[i], custom_comparisons = NULL) 
    dsg_results <- list()
    for(dsg_i in dsg){
      if(dsg_i %in% c("repeat_and_pool", "repeat_and_replace", "bayes")){
        data_stage1 <- generated_data$data_stage1
      } else {data_stage1 <- generated_data$data_singlestage}
      
      dsg_results[[which(dsg_i == dsg)]] <- run_simulation(dsg = dsg_i, n_per_group = param_grid$n_per_group[i], no_groups = param_grid$no_groups[i],
                                                           data_stage1 = data_stage1, data_stage2 = generated_data$data_stage2,
                                                           distribution = param_grid$distribution[i], type_of_comparison = param_grid$type_of_comparison[i],
                                                           stanmodel = stanmodel, resampling_distributions = NULL)
    }
    sim_results[[sim]] <- bind_rows(dsg_results) %>% add_column(true_effect = param_grid$effect[i])
  }
  return(bind_rows(sim_results))
}
stopImplicitCluster()
endtime <- Sys.time()
difftime(endtime, starttime)

saveRDS(data, file = paste("simulated_data_lognormal_", Sys.Date(), ".Rds", sep = ""))


##############################################################################################
