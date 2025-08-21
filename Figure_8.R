# confidence interval distribution for all designs and hypotheses
data %>% filter(n_per_group == 12 & no_groups == 10 & distribution == "normal") %>% 
  group_by(true_effect, dsg) %>% 
  arrange(no_groups, n_per_group, CI_lwr) %>% 
  mutate(rowid = 1:n()) %>% ungroup() %>% 
  mutate(CI_excludes_0 = case_when(CI_lwr < true_effect & CI_upr > true_effect ~ FALSE, TRUE ~ TRUE)) %>% 
  mutate(CI_excludes_0f = factor(CI_excludes_0, levels = c(TRUE, FALSE)),
         length_CI = CI_upr - CI_lwr) %>% 
  ggplot(.) + geom_errorbar(aes(y = rowid, xmin = CI_lwr, xmax = CI_upr, colour = CI_excludes_0f)) + 
  scale_colour_manual(values = c("red4", "grey50"), guide = "none") + 
  facet_grid(true_effect ~ dsg, labeller = scenario_labeller3) + 
  xlab("Effect estimate") + theme_bw() + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  geom_vline(aes(xintercept = true_effect), linetype = "dashed") +
  geom_text(aes(x = -2, y = 9000, label = coverage), size = 2.5, 
            data = . %>% group_by(true_effect, dsg) %>% 
              summarise(coverage = round(sum(!CI_excludes_0)/n(),2)) %>% 
              mutate(coverage = paste0("CP: ", coverage))
  ) +
  geom_text(aes(x = -2, y = 7000, label = al), size = 2.5, 
            data = . %>% group_by(true_effect, dsg) %>% 
              summarise(al = round(mean(length_CI),2)) %>% 
              mutate(al = paste0("AL: ", al))
  )

ggsave(filename = "figures/Figure 8.png", device = "png", width = 7, height = 2.5)
