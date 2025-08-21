# function to estimate the effect from generated data
estimate_effect <- function(data_stage1, data_stage2 = NULL, n_per_group, dsg, 
                            type_of_comparison, multiple_test = NULL, stanmodel, delta = NULL){
  
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
    
    historical_data <- bind_rows(data.frame(y = data_stage1[,2], x = 0), 
                                 data.frame(y = data_stage1[,1], x = 1))
    historical_fit <- lm(y ~ x, data = historical_data)
    
    fit <- rstan::sampling(
      stanmodel, 
      data = list("df_beta" = historical_fit$df.residual, 
                  "y" = c(data_stage2[,2], data_stage2[,1]), 
                  "group" = c(rep(0,length(data_stage2[,2])), rep(1,length(data_stage2[,1]))), 
                  "N" = length(c(data_stage2[,2], data_stage2[,1])), 
                  "delta" = delta, 
                  "diff_prior" = as.numeric(historical_fit$coefficients[2]), 
                  "sigma_prior" = summary(historical_fit)$coefficients[2,2], 
                  "diff_skeptical" = 0, 
                  "sigma_skeptical" = summary(historical_fit)$coefficients[2,2]), 
      refresh = 0
      )      
    posterior_mixture <- rstan::extract(fit) %>% as_tibble()
    
    effect_est <- mean(posterior_mixture$beta)
    CI <- quantile(posterior_mixture$beta, probs = c(0.025, 0.975))
    best_pvalue <- min(sum(posterior_mixture$beta < 0)/nrow(posterior_mixture), 
                       sum(posterior_mixture$beta > 0)/nrow(posterior_mixture))*2
  }
  
  return(list("effect_est" = effect_est, "CI_lwr" = CI[1], "CI_upr" = CI[2], 
              "best_pvalue" = best_pvalue, "n_per_group" = n_per_group, "dsg" = dsg))
  
}

