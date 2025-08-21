# function to sample data 
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
