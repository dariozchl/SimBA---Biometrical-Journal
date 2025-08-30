# confidence interval distribution
data %>% 
  filter(dsg == "standard" & true_effect == 0 & no_groups %in% c(2, 6, 10) & distribution == "normal") %>% 
  group_by(no_groups, n_per_group) %>% 
  arrange(CI_lwr) %>% 
  mutate(rowid = 1:n()) %>% ungroup() %>% 
  mutate(CI_excludes_0 = case_when(CI_lwr < 0 & CI_upr > 0 ~ FALSE, TRUE ~ TRUE)) %>% 
  mutate(CI_excludes_0f = factor(CI_excludes_0, levels = c(TRUE, FALSE)),
         length_CI = CI_upr - CI_lwr) %>% 
  ggplot(.) + geom_errorbar(aes(y = rowid, xmin = CI_lwr, xmax = CI_upr, colour = CI_excludes_0f)) + 
  scale_color_manual(values = c("red4", "grey50"), guide = "none") + 
  facet_grid(no_groups ~ n_per_group, labeller = scenario_labeller2) + 
  xlab("Effect estimate") + theme_bw() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  geom_text(aes(x = -3, y = 9000, label = coverage), size = 2.5, 
            data = . %>% group_by(no_groups, n_per_group) %>% 
              summarise(coverage = round(sum(!CI_excludes_0)/n(),2)) %>% 
              mutate(coverage = paste0("CP: ", coverage))
  ) +
  geom_text(aes(x = -3, y = 7000, label = al), size = 2.5, 
            data = . %>% group_by(no_groups, n_per_group) %>% 
              summarise(al = round(mean(length_CI),2)) %>% 
              mutate(al = paste0("AL: ", al))
  )

ggsave(filename = "figures/Figure 3.png", device = "png", width = 7, height = 2.5)
