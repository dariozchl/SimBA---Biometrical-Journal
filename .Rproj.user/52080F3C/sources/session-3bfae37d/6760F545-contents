
main_data <- data

main_data <- main_data %>% 
  mutate(delta = case_when(dsg == "bayes" ~ 0.5,
                           TRUE ~ as.numeric(NA)))

data_supplementary_material$delta <- rep(param_grid_suppl$delta, each = nsim)

data <- bind_rows(data_supplementary_material, main_data)


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

# create Figure 9
ploteff0 <- data %>% filter(true_effect == 0 & distribution == "normal" & dsg == "bayes") %>% 
  group_by(delta, no_groups, n_per_group, distribution, type_of_comparison, true_effect) %>% 
  summarise(MSE = mean(squared_error)) %>% 
  ggplot(.) + geom_point(aes(x = no_groups, y = MSE), size = 1) + 
  facet_nested(n_per_group ~ delta) + 
  expand_limits(y=0) +
  theme_bw() + ggtitle("True largest effect: d = 0") + xlab("Number of groups")
ploteff1 <- data %>% filter(true_effect == 1 & distribution == "normal" & dsg == "bayes") %>% 
  group_by(delta, no_groups, n_per_group, distribution, type_of_comparison, true_effect) %>% 
  summarise(MSE = mean(squared_error)) %>% 
  ggplot(.) + geom_point(aes(x = no_groups, y = MSE), size = 1) +
  facet_nested(n_per_group ~ delta , labeller = scenario_labeller) + 
  expand_limits(y=0) +
  theme_bw() + ggtitle("True largest effect: d = 1") + xlab("Number of groups")
ploteff2 <- data %>% filter(true_effect == 2 & distribution == "normal" & dsg == "bayes") %>% 
  group_by(delta, no_groups, n_per_group, distribution, type_of_comparison, true_effect) %>% 
  summarise(MSE = mean(squared_error)) %>% 
  ggplot(.) + geom_point(aes(x = no_groups, y = MSE), size = 1) +
  facet_nested(n_per_group ~ delta , labeller = scenario_labeller) + 
  expand_limits(y=0) +
  theme_bw() + ggtitle("True largest effect: d = 2") + xlab("Number of groups")

ggarrange(ploteff0, ploteff1, ploteff2, ncol=1)
ggsave(filename = "figures/SupplMater_Figure1.png", device = "png", width = 7, height = 7)



#######

# traffic light system under alternative hypothesis

# for each n_per_group, the MSE should be calculated as percentage from the MSE of dsg=="standard design" and no_groups==2
plot1 <- data %>% filter(true_effect == 1, distribution == "normal" & type_of_comparison == "all_pairwise" & dsg != "resampled_shrinkage") %>% 
  group_by(dsg, no_groups, n_per_group, distribution, type_of_comparison) %>% 
  summarise(MSE = mean(squared_error)) %>% 
  ungroup() %>% group_by(n_per_group) %>% 
  mutate(MSE_ref = MSE/MSE[dsg=="standard" & no_groups==2]) %>% 
  mutate(ampelfarbe = case_when(MSE_ref <= MSE_ref[dsg=="standard" & no_groups==2]*1.05~ "green",
                                MSE_ref <= MSE_ref[dsg=="standard" & no_groups==3]*1.05 ~ "yellow",
                                TRUE ~ "red"),
         dsg = factor(dsg, levels = c("standard", "repeat_and_pool", "repeat_and_replace", "bayes"))) %>% 
  ggplot(.) + geom_point(aes(x = no_groups, y = MSE_ref)) +
  geom_rect(aes(xmin = no_groups-0.499, xmax = no_groups+0.501, ymin=0, ymax=Inf, fill = ampelfarbe), alpha = 0.5) +
  facet_grid(n_per_group ~ dsg, labeller = scenario_labeller) +
  scale_fill_manual(values = c("#009e73", "#f0e442", "#d55e00"), breaks = c("green", "yellow", "red"), guide="none") +
  scale_x_continuous(breaks = 2:10) + xlab("Number of groups") + ylab("MSE (relative to reference design)") +
  theme_bw() + theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  ggtitle("Traffic light system under alternative hypothesis d = 1")
plot1

ggsave(filename = "figures/SupplMater_Figure2.png", device = "png", width = 10, height = 3.5)


plot2 <- data %>% filter(true_effect == 2, distribution == "normal" & type_of_comparison == "all_pairwise" & dsg != "resampled_shrinkage") %>% 
  group_by(dsg, no_groups, n_per_group, distribution, type_of_comparison) %>% 
  summarise(MSE = mean(squared_error)) %>% 
  ungroup() %>% group_by(n_per_group) %>% 
  mutate(MSE_ref = MSE/MSE[dsg=="standard" & no_groups==2]) %>% 
  mutate(ampelfarbe = case_when(MSE_ref <= MSE_ref[dsg=="standard" & no_groups==2]*1.05~ "green",
                                MSE_ref <= MSE_ref[dsg=="standard" & no_groups==3]*1.05 ~ "yellow",
                                TRUE ~ "red"),
         dsg = factor(dsg, levels = c("standard", "repeat_and_pool", "repeat_and_replace", "bayes"))) %>% 
  ggplot(.) + geom_point(aes(x = no_groups, y = MSE_ref)) +
  geom_rect(aes(xmin = no_groups-0.499, xmax = no_groups+0.501, ymin=0, ymax=Inf, fill = ampelfarbe), alpha = 0.5) +
  facet_grid(n_per_group ~ dsg, labeller = scenario_labeller) +
  scale_fill_manual(values = c("#009e73", "#f0e442", "#d55e00"), breaks = c("green", "yellow", "red"), guide="none") +
  scale_x_continuous(breaks = 2:10) + xlab("Number of groups") + ylab("MSE (relative to reference design)") +
  theme_bw() + theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  ggtitle("Traffic light system under alternative hypothesis d = 2")
plot2

ggsave(filename = "figures/SupplMater_Figure3.png", device = "png", width = 10, height = 3.5)
