# function to find the "best" pairwise comparison among the comparisons specified
find_best_comparison <- function(samples, type_of_comparison, custom_comparisons = NULL){
  
  if(type_of_comparison == "all_pairwise"){
    best_comparison <- samples %>% pivot_longer(., cols=everything(), names_to = "group") %>% 
      compare_means(value ~ group, data=., method = "t.test") %>% dplyr::select(group1, group2, p) %>% slice(which.min(p))
    best_groups <- as.numeric(gsub(pattern = "V", replacement = "", x = c(best_comparison$group1, best_comparison$group2)))
    results <- best_comparison
  } 
  
  else if(type_of_comparison == "all_vs_control"){
    glht_fit <- glht(aov(value ~ name, 
                         data = samples %>% pivot_longer(cols=everything()) %>% mutate(name = as.factor(name))), 
                     linfct = mcp(name = "Dunnett"))
    glht_fit_summary <- summary(glht_fit)
    # need to get the correct indices of the best comparison from glht-object
    best_comparison_name <- names(glht_fit_summary$test$tstat)[which.min(glht_fit_summary$test$pvalues)]
    best_groups <- as.numeric(gsub(pattern = "V", replacement = "", x = unlist(strsplit(best_comparison_name, split = "-"))))
    results <- glht_fit_summary
  } 
  
  else if(type_of_comparison == "custom"){# custom comparison
    # the trick is to:
    # 1. access the list of all possible contrasts,
    # 2. get the indices of the defined contrasts and transform into emmeans language,
    # 3. keep only the indices of the contrasts that were specified.
    custom_comparisons_emmeans <- gsub("(\\d+)", "V\\1", gsub(" vs\\. ", " - ", custom_comparisons))
    fit <- lm(value ~ name, 
              data = samples %>% pivot_longer(cols = everything()) %>% mutate(name = as.factor(name)))
    # if adjust = "none" is not specified, the Sidak method is used to adjust p-values
    all_contrasts <- summary(emmeans(fit, pairwise ~ name)$contrasts, adjust = "none") %>%
      filter(contrast %in% custom_comparisons_emmeans) 
    best_pvalue_dsg1 <- min(all_contrasts[,"p.value"])
    best_groups <- as.numeric(gsub(pattern = "V", replacement = "", x = unlist(strsplit(all_contrasts[which.min(all_contrasts[,"p.value"]),]$contrast, " - "))))
    results <- all_contrasts
  }
  
  return(list("best_groups"=best_groups, "results"=results))
}
