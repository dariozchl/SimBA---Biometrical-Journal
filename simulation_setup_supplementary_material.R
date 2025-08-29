############ simulation setup

# number of simulation runs per scenario
nsim <- 1e4

# define simulation scenarios
dsg <- c("bayes")
n_per_group <- c(6, 12, 18)
no_groups <- 2:10
distribution <- c("normal")
effect <- c(0,1,2)
type_of_comparison <- c("all_pairwise")
delta <- c(0.25, 0.75)

# create a grid containing all simulation scenarios
param_grid_suppl <- expand.grid(
  # dsg = dsg,
  n_per_group = n_per_group,
  no_groups = no_groups,
  distribution = distribution,
  type_of_comparison = type_of_comparison,
  effect = effect,
  delta = delta,
  stringsAsFactors = FALSE
)

# add sample sizes for two-stage designs
param_grid_suppl <- param_grid_suppl %>% 
  mutate(n_per_group_stage1 = floor((n_per_group * no_groups) / (no_groups+2))) %>% 
  rowwise() %>% 
  mutate(n_per_group_stage2 = floor(((n_per_group * no_groups) - (n_per_group_stage1 * no_groups)) / 2))


##### preparations 
stanmodel <- stan_model("mixture_model_ttest.stan")

# define parameters for effect & number of groups 
generate_parameters <- function(effect, no_groups) {
  return(mapply(c, c(effect, rep(0, no_groups - 1)), rep(1, no_groups), SIMPLIFY = FALSE))
}

# add parameters to the grid 
param_grid_suppl$parameters_new <- mapply(generate_parameters, param_grid_suppl$effect, param_grid_suppl$no_groups)


# make a grid with all simulation runs
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...)) # expand.grid for data.frames
param_grid_suppl_sim <- expand.grid.df(param_grid_suppl, data.frame("simID" = 1:nsim)) %>% as_tibble()

param_grid_suppl_sim <- expand.grid.df(
  param_grid_suppl,
  expand.grid(simid = 1:nsim,
              dsg = dsg)
)
