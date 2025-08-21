# effect estimate distribution for all designs and hypotheses
data %>% 
  filter(no_groups == 10 & distribution == "normal") %>% 
  group_by(dsg, no_groups, n_per_group) %>% 
  #mutate(no_groups = as.factor(no_groups)) %>% 
  ggplot(.) + geom_histogram(aes(x = effect_est), bins = 100) + ylab("Number of simulations") + xlab("Effect estimate") +
  facet_grid(true_effect ~ dsg, labeller = scenario_labeller3) + theme_bw()
ggsave(filename = "figures/Figure 7.png", device = "png", width = 7, height = 2.5)
