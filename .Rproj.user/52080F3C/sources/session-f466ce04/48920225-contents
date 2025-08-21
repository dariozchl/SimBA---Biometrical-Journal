data <- data %>% mutate(dsg = factor(dsg, 
                                     levels = c("standard", "repeat_and_pool", 
                                                "repeat_and_replace", "bayes"))) %>% 
  rowwise() %>% 
  mutate(squared_error = diff(c(true_effect, effect_est))^2) 


scenario_labeller <- labeller(`dsg` = c(`repeat_and_pool` = "Repeat and pool", 
                                        `repeat_and_replace` = "Repeat and replace", 
                                        `standard` = "Standard",
                                        `bayes` = "Bayesian robust mixture"),
                              `n_per_group` = c(`6` = "n = 6", 
                                                `12` = "n = 12", 
                                                `18` = "n = 18"),
                              `true_effect` = c(`0` = "d = 0", `1` = "d = 1", `2` = "d = 2")) 

scenario_labeller2 <- labeller(`no_groups` = c(`2` = "k = 2", `6` = "k = 6", `10` = "k = 10"),
                               `n_per_group` = c(`6` = "n = 6", `12` = "n = 12", `18` = "n = 18")) 

scenario_labeller3 <- labeller(`true_effect` = c(`0` = "d = 0", `1` = "d = 1", `2` = "d = 2"),
                               `dsg` = c(`repeat_and_pool` = "Repeat and pool", 
                                         `repeat_and_replace` = "Repeat and replace", 
                                         `standard` = "Standard", 
                                         `bayes` = "Bayesian robust mixture")) 
