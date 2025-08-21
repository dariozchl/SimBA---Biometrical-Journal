# run simulation
(starttime <- Sys.time())
function_in_parallel <- function(i){
  sim_results <- list()
  for(sim in 1:nsim){
    generated_data <- get_data(
      n_per_group = param_grid$n_per_group[i], no_groups = param_grid$no_groups[i],
      distribution = param_grid$distribution[i], parameters_new = param_grid$parameters_new[[i]],
      type_of_comparison = param_grid$type_of_comparison[i], custom_comparisons = NULL
      ) 
    dsg_results <- list()
    for(dsg_i in dsg){
      if(dsg_i %in% c("repeat_and_pool", "repeat_and_replace", "bayes")){
        data_stage1 <- generated_data$data_stage1
      } else {data_stage1 <- generated_data$data_singlestage}
      
      dsg_results[[which(dsg_i == dsg)]] <- run_simulation(
        dsg = dsg_i, n_per_group = param_grid$n_per_group[i], no_groups = param_grid$no_groups[i],
        data_stage1 = data_stage1, data_stage2 = generated_data$data_stage2,
        distribution = param_grid$distribution[i], type_of_comparison = param_grid$type_of_comparison[i],
        stanmodel = stanmodel, delta = param_grid$delta[i]
        )
    }
    sim_results[[sim]] <- bind_rows(dsg_results) %>% add_column(true_effect = param_grid$effect[i])
  }
  return(bind_rows(sim_results))
}


future::plan(multisession, workers = length(availableWorkers())-1)
sim_results_list <- future_lapply(1:nrow(param_grid), FUN = function_in_parallel, 
                                  future.seed = 19082025)
data <- bind_rows(sim_results_list)

endtime <- Sys.time()
print(difftime(endtime, starttime))

saveRDS(data, file = paste("simulated_data_", Sys.Date(), ".Rds", sep = ""))

