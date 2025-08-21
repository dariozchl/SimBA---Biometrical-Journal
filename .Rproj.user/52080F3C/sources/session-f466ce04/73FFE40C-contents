# create Figure 9
ploteff0 <- data %>% filter(true_effect == 0 & distribution == "normal") %>% 
  group_by(dsg, no_groups, n_per_group, distribution, type_of_comparison, true_effect) %>% 
  summarise(MSE = mean(squared_error)) %>% 
  ggplot(.) + geom_point(aes(x = no_groups, y = MSE), size = 1) + 
  facet_nested(n_per_group ~ dsg , labeller = scenario_labeller) + 
  scale_y_continuous(breaks = c(0, 0.5, 1)) + expand_limits(y=0) +
  theme_bw() + ggtitle("True largest effect: d = 0") + xlab("Number of groups")
ploteff1 <- data %>% filter(true_effect == 1 & distribution == "normal") %>% 
  group_by(dsg, no_groups, n_per_group, distribution, type_of_comparison, true_effect) %>% 
  summarise(MSE = mean(squared_error)) %>% 
  ggplot(.) + geom_point(aes(x = no_groups, y = MSE), size = 1) +
  facet_nested(n_per_group ~ dsg , labeller = scenario_labeller) + 
  scale_y_continuous(breaks = c(0, 0.5, 1)) + expand_limits(y=0) +
  theme_bw() + ggtitle("True largest effect: d = 1") + xlab("Number of groups")
ploteff2 <- data %>% filter(true_effect == 2 & distribution == "normal") %>% 
  group_by(dsg, no_groups, n_per_group, distribution, type_of_comparison, true_effect) %>% 
  summarise(MSE = mean(squared_error)) %>% 
  ggplot(.) + geom_point(aes(x = no_groups, y = MSE), size = 1) +
  facet_nested(n_per_group ~ dsg , labeller = scenario_labeller) + 
  scale_y_continuous(breaks = c(0, 0.3, 0.6)) + expand_limits(y=0) +
  theme_bw() + ggtitle("True largest effect: d = 2") + xlab("Number of groups")

ggarrange(ploteff0, ploteff1, ploteff2, ncol=1)
ggsave(filename = "figures/Figure 9.png", device = "png", width = 7, height = 7)
