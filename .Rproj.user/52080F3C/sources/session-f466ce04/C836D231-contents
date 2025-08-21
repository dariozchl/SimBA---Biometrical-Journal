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