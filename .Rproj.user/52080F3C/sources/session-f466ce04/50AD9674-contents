# function wrapper to run simulation
run_simulation <- function(dsg, data_stage1, data_stage2, n_per_group, no_groups, distribution, 
                           type_of_comparison, stanmodel = NULL, delta = NULL){
  
  # use pmap to iterate through rows of the grid
  data <- estimate_effect(
    data_stage1 = data_stage1, data_stage2 = data_stage2, n_per_group = n_per_group, 
    dsg = dsg, type_of_comparison = type_of_comparison, multiple_test = NULL, 
    stanmodel = stanmodel, delta = delta
  )
  return(data %>% bind_rows() %>% 
           add_column("no_groups" = no_groups, 
                      "distribution" = distribution, 
                      "type_of_comparison" = type_of_comparison))
}
