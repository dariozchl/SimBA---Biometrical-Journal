library(shiny)
library(shinyjs)

##################################################################################################
ui <- fluidPage(
  
  # App title ----
  titlePanel("Simulation of estimation error in complex animal experimental designs."),
  
  sidebarLayout(
    
    sidebarPanel(
      
      numericInput("no_groups", "Number of groups:",
                   min = 0, max = 100, value = 5),
      
      numericInput("n_per_group", "N per group:",
                   min = 3, max = 30, value = 10),
      
      radioButtons(inputId = "define_comparisons", label = "How should the comparisons be defined?", 
                   choices = c("All vs. control" = "all_vs_control",
                               "All possible pairwise comparisons" = "all_pairwise", 
                               "User-defined" = "custom")),
      
      uiOutput(outputId = "define_effect"),

      uiOutput(outputId = "custom_comparisons"),
      
      uiOutput(outputId = "choose_effect"),
      
      uiOutput(outputId = "choose_distribution"),
      
      uiOutput(outputId = "choose_distribution_custom"),

      uiOutput(outputId = "user_def_par1"),
      
      uiOutput(outputId = "user_def_par2"),
      
      radioButtons("multiple_test", "Adjust for multiple testing?",
                   choices = c("Yes", "No"), selected = "No"),
      
      uiOutput(outputId = "procedure_choices"),
      withMathJax(),
      
      numericInput("nsim", "How many simulation runs should be performed?",
                   min = 0, max = 100, value = 30),
      
      radioButtons(inputId = "choose_bayes", "Include Bayesian mixture model in the simulations? (increases computation time)",
                   choices = c("Yes", "No"), selected = "No"),
      
      uiOutput(outputId = "mixing_weight"),
      
      radioButtons(inputId = "traffic_light_system", "Show traffic light system?",
                   choices = c("Yes", "No")),
      
      uiOutput(outputId = "tls_cust_text"),
      
      uiOutput(outputId = "tls_no_groups"),
      
      uiOutput(outputId = "tls_n_per_group"),
      
      uiOutput(outputId = "tls_reference_green"),
      
      uiOutput(outputId = "tls_reference_yellow"),
      
      uiOutput(outputId = "check_tls"),
      
      uiOutput(outputId = "check_mcp"),
      
      useShinyjs(),
      
      actionButton("run", "Run")
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      plotOutput(outputId = "plot"),
      textOutput("rejection_rates")
    )  
  )
)



##################################################################################################
server <- function(input, output) {
  
  library(tidyverse)
  library(ggpubr)
  library(gridExtra)
  library(multcomp)
  library(emmeans)
  library(rstan)
  library(purrr)
  library(foreach)
  library(rlist)
  
  
  
  # Function to generate all possible pairwise comparisons
  generate_pairwise <- function(k) {
    combinations <- combn(1:k, 2, simplify = FALSE)
    # Convert each combination to a string
    combo_strings <- sapply(combinations, function(combo) paste(combo, collapse = " vs. "))
    return(combo_strings)
  }
  
  # Create a reactive expression to generate combinations based on number of groups k
  combos <- reactive({
    k <- input$no_groups
    if (is.null(k) | k < 1) return(NULL)
    generate_pairwise(k)
  })
  
  # Render the relevant group comparisons dynamically
  output$custom_comparisons <- renderUI({
    if (input$define_comparisons == "custom") {
      selectInput(inputId = "custom_comparisons", "Select all relevant group comparisons:", 
                  choices = combos(), multiple = TRUE)
    }
  })
  
  output$define_effect <- renderUI({
    if (input$define_comparisons == "all_vs_control") {
      radioButtons(inputId = "define_effect", label = "How should the assumed effect size be defined?", 
                   choices = c("All comparisons have the same assumed effect size" = "all_equal", 
                               "User-defined data generating process for each group" = "custom"))
    } 
  })
  
  # Render the effect size slider dynamically
  output$choose_effect <- renderUI({
    req(input$define_comparisons) # to avoid evaluation before an input is given
    req(input$define_effect)
      if (input$define_comparisons == "all_vs_control" & input$define_effect == "all_equal") {
      sliderInput(inputId = "choose_effect", label = "Assumed effect size (Cohen's d)", 
                  min = 0, max = 10, value = 0, step = 0.1)
    }
  })
  
  # Render the distribution radio button dynamically 
  # distributional parameters are not needed if the same effect size is assumed for all comparisons
  output$choose_distribution <- renderUI({
    req(input$define_comparisons)
    req(input$define_effect)
    if (input$define_comparisons == "all_vs_control" & input$define_effect == "all_equal") {
      radioButtons(inputId = "choose_distribution", label = "Which distribution is assumed for the data?", 
                   choices = c("Normal distribution" = "normal", 
                               "Lognormal distribution" = "lognormal",
                               "Cauchy distribution" = "cauchy"))
    }
  })
  output$choose_distribution_custom <- renderUI({
    req(input$define_comparisons)
    req(input$define_effect)
    if (!(input$define_comparisons == "all_vs_control" & input$define_effect == "all_equal")) {
      radioButtons(inputId = "choose_distribution_custom", label = "Which distribution is assumed for the data?", 
                   choices = c("Normal distribution (mean and SD)" = "normal", 
                               "Lognormal distribution (mean and SD on logscale)" = "lognormal",
                               "Cauchy distribution (location and scale)" = "cauchy"))
    }
  })
  
  # distributional parameters are not needed if the same effect size is assumed for all comparisons
  output$user_def_par1 <- renderUI({
    req(input$define_comparisons)
    req(input$define_effect)
    if (!(input$define_comparisons == "all_vs_control" & input$define_effect == "all_equal")) {
      textInput(inputId = "user_def_par1", 
                label = "First distributional parameter for each group delimited by comma", 
                value = paste(rep("0", input$no_groups), collapse = ", "))
    }
  })

  output$user_def_par2 <- renderUI({
    req(input$define_comparisons)
    req(input$define_effect)
    if (!(input$define_comparisons == "all_vs_control" & input$define_effect == "all_equal")) {
      textInput(inputId = "user_def_par2", 
                label = "Second distributional parameter for each group delimited by comma", 
                value = paste(rep("1", input$no_groups), collapse = ", "))
    }
  })
  
  output$procedure_choices <- renderUI({
    if (input$multiple_test == "No") return(NULL)
    
    if (input$define_comparisons == "all_vs_control") {
      choices <- c("Bonferroni", "Holm", "Hochberg", "Hommel", "Benjamini-Hochberg", "Dunnett", "Step-Down Dunnett")
    } else if  (input$define_comparisons == "all_pairwise") {
      choices <- c("Bonferroni", "Holm", "Hochberg", "Hommel", "Benjamini-Hochberg", "Tukey", "S2", "Step-Down Tukey")
    } else if (input$define_comparisons == "custom") {
      choices <- c("Bonferroni", "Holm", "Hochberg", "Hommel", "Benjamini-Hochberg", "S2")
    }
    
  output$mixing_weight <- renderUI({
    if (input$choose_bayes == "Yes"){
      sliderInput(inputId = "mixing_weight", "Determine the mixing weight:",
                  min = 0, max = 1, value = 0.5, step = 0.05)
    }
      
  })
  
  output$tls_cust_text <- renderUI({
    if (input$traffic_light_system == "Yes") {
      if (input$define_comparisons == "custom") {
        return(p("Note: the reference designs in the traffic light system apply all pairwise comparisons."))
      }
    }
  })
    
  output$tls_no_groups <- renderUI({
    if (input$traffic_light_system == "Yes") {
      if (input$define_comparisons == "custom") {
        max <- input$no_groups-1 # possible no_groups so there are no overlaps
        selectInput(inputId = "tls_no_groups", "Select number of groups to be included in the traffic light system:",
                    choices = c(2:max), multiple = TRUE)
      } else {
      selectInput(inputId = "tls_no_groups", "Select number of groups to be included in the traffic light system:",
                  choices = c(2:100), multiple = TRUE)
      }
    }
  }) 
  
  output$tls_n_per_group <- renderUI({
    if (input$traffic_light_system == "Yes") {
      
      selectInput(inputId = "tls_n_per_group", "Select group sizes to be included in the traffic light system:",
                  choices = c(3:30), multiple = TRUE)
    }
  })
  
  output$tls_reference_green <- renderUI({
    if (input$traffic_light_system == "Yes") {
      if (input$define_comparisons == "custom") {
        max <- input$no_groups-1 # possible no_groups so there are no overlaps
        selectInput(inputId = "tls_reference_green", "Select number of groups to be the reference for the upper green limit:",
                    choices = c(2:max), selected = 2)
      } else {
      selectInput(inputId = "tls_reference_green", "Select number of groups to be the reference for upper green limit:",
                  choices = c(2:100))
      
      }
    }
  })
  
  output$tls_reference_yellow <- renderUI({
    if (input$traffic_light_system == "Yes") {
      if (input$define_comparisons == "custom") {
        max <- input$no_groups-1 # possible no_groups so there are no overlaps
        selectInput(inputId = "tls_reference_yellow", "Select number of groups to be the reference for the upper yellow limit:",
                    choices = c(2:max), selected = 3)
      } else {
      selectInput(inputId = "tls_reference_yellow", "Select number of groups to be the reference for upper yellow limit:",
                  choices = c(2:100), selected = 3)
      }
      
    }
  })
  
  
    
  output$check_tls <- renderUI({
    if (input$traffic_light_system == "Yes") { 
      
      # if user defined parameters (means) exist, check whether they are equal (-> null hypothesis)
      if (!(is.null(input$user_def_par1))){
      # unlist user_def_par1 
      user_parameters <- as.numeric(unlist(strsplit(x = input$user_def_par1, split = ",")))
      
      if (length(unique(user_parameters)) > 1) {
        return(p("Traffic light system is only valid under the null hypothesis."))
        
      }}
      
      # if Cohen's d was given as an input, check whether it is 0 (-> null hypothesis)  
      req(input$choose_effect)
      if (input$choose_effect != 0) {
        return(p("Traffic light system is only valid under the null hypothesis."))
      }
    }
  })
  
  output$check_mcp <- renderUI({
    if (input$multiple_test == "Yes") {
      if (length(selected_methods()) == 0) {
        return(p("To adjust for multiple testing, at least one procedure has to be selected."))
      }
    }
  })
  
  # "run" button is disabled if 
  # 1. traffic light system is chosen but simulation is not under null hypothesis OR
  # 2. adjusting for multiple comparisons is chosen but no procedure to be applied is selected 
  observe({
    
    tls_invalid <- FALSE
    mcp_invalid <- FALSE
    
    if (input$traffic_light_system == "Yes") {
      if (!(is.null(input$user_def_par1))){
        user_parameters <- as.numeric(unlist(strsplit(x = input$user_def_par1, split = ",")))
        
        if (length(unique(user_parameters)) > 1) {
          tls_invalid <- TRUE
        }}
      
      if (!(is.null(input$choose_effect)) && input$choose_effect != 0){
        tls_invalid <- TRUE
      }
    }
    
    if (input$multiple_test == "Yes") {
      procedures <- selected_methods()
      if (length(procedures) == 0){
        mcp_invalid <- TRUE
      }
    }

    if (tls_invalid == TRUE || mcp_invalid == TRUE){
      shinyjs::disable("run")
    } else {
      shinyjs::enable("run")
    }
    
  })
  
    
    checkboxes_with_info <- lapply(choices, function(proc) {
      # for each procedure, a (fluid) row with three columns is created
      fluidRow(
        # checkbox 
        column(1, 
               checkboxInput(
                 inputId = paste0("check_", proc),
                 label = NULL,
                 value = FALSE
               )
        ),
        # label
        column(7, tags$label(`for` = paste0("check_", proc), proc)),
        # info icon
        column(1, actionLink(
          inputId = paste0("info_", proc),
          label = icon("info-circle")
        ))
      )
    })
    
    tagList(
      tags$label("Select at least one procedure to be applied:"),
      tags$div(style = "margin-bottom: 10px;", checkboxes_with_info)
    )
  })
  
  # return procedures that were chosen
  selected_methods <- reactive({
    if (input$define_comparisons == "all_vs_control") {
      choices <- c("Bonferroni", "Holm", "Hochberg", "Benjamini-Hochberg", "Hommel", "Dunnett", "Step-Down Dunnett")
    } else if  (input$define_comparisons == "all_pairwise") {
      choices <- c("Bonferroni", "Holm", "Hochberg", "Benjamini-Hochberg", "Hommel", "Tukey", "S2", "Step-Down Tukey")
    } else if (input$define_comparisons == "custom") {
      choices <- c("Bonferroni", "Holm", "Hochberg", "Benjamini-Hochberg", "Hommel", "S2")
    }
    
    Filter(function(proc) {
      isTRUE(input[[paste0("check_", proc)]])
    }, choices)
  })
  
  
  
  procedure_info_text <- function(procedure) {
    switch(procedure,
           "Bonferroni" = withMathJax("Single-step procedure that controls the FWER. 
                                      The unadjusted p-values \\(p_1, p_2, ..., p_m\\) are compared with the threshold \\(\\frac{\\alpha}{m}\\)."),
           "Holm" = withMathJax("Step-down version of the Bonferroni that controls the FWER and is at least as powerful as Bonferroni. 
                                The smallest p-value \\(p_1\\)
                                is tested first: if \\(p_1 > \\frac{\\alpha}{m}\\), the procedure stops without any rejections.
                                Otherwise, \\(H_{(1)}\\) is rejected and \\(H_{(2)}\\) is tested at the larger significance level \\(\\frac{\\alpha}{(m-1)}\\).
                                These steps are repeated until either the first non-rejection occurs or all hypotheses \\(H_{(1)},...,H_{(m)}\\) are rejected."),
           "Hochberg" = withMathJax("Step-up extension of the Simes test that can be considered a reversed Holm procedure. It is at least as powerful as Holm.
                                    The null hypothesis \\(H_{(m)}\\) associated with the largest p-value \\(p_m\\) is tested first:
                                    if \\(p_m \\leq {\\alpha}\\ \\), all hypotheses \\(H_{(1)},...,H_{(m)}\\) are rejected. 
                                    Otherwise, \\(H_{(m)}\\) is retained and \\(H_{(m-1)}\\) is tested at the smaller significance level \\(\\frac{\\alpha}{2}\\).
                                    If \\(p_{(m-1)} \\leq \\frac{\\alpha}{2}\\), the procedure stops and all hypotheses \\(H_{(1)},...,H_{(m-1)}\\) are rejected.
                                    Theses steps are repeated until either the first rejection occurs or all null hypotheses \\(H_{(1)},...H_{(m)}\\) are retained."),
           "Hommel" = withMathJax(" Improvement of and therefore more powerful than the Hochberg procedure.
                                    For a system of \\(m\\) hypotheses, let \\(S\\) denote the set of all \\(j \\in \\{1,\\ldots,m\\}\\) so that exactly \\(j\\) of the \\(m\\)
                                    hypotheses can simultaneously be true, while the remaining  \\((m-j)\\) hypotheses are false. 
                                    \\(j\\) is computed as \\(j = \\max\\{i \\in S: p_{(m-i+k)} > \\frac{k\\cdot \\alpha}{i}\\}\\) for \\(k = 1,\\ldots,i\\). 
                                    If this maximum does not exist, all \\(H_i\\) are rejected. Otherwise, all \\(H_{i}\\) with \\(p_i \\leq \\frac{\\alpha}{j}\\) 
                                    are rejected. 
                                    If all \\(m\\) tests are independent, the Hommel procedure controls the FWER."),
           "Benjamini-Hochberg" = withMathJax("The Benjamini-Hochberg procedure controls the FDR at level \\(q^*\\). Let \\(k\\) denote the largest \\(i\\) for which 
                                              \\(p_{(i)} \\leq \\frac{i}{m}q^*\\), then all \\(H_{(i)}, i = 1,2,\\ldots,k\\) can be rejected."),
           "Dunnett" = withMathJax("Standard procedure for many-to-one comparisons (e.g. comparing several treatments to a control). 
                                   Assuming balanced groups and homoscedasticity, pairwise t-tests are computed: 
                                   \\(t_i = \\frac{\\bar{y_i}-\\bar{y_0}}{s\\sqrt{\\frac{1}{n_i}+\\frac{1}{n_0}}}\\), \\(i = 1,\\ldots,m\\), where \\(\\bar{y_i}\\) 
                                   denotes the arithmetic mean of group \\(i = 0,\\ldots,m\\), \\(i = 0\\) denotes the control group and \\(s^2\\) is the pooled 
                                   variance estimate. Depending on the sideness of the test problem, the Dunnett test takes the minimum or the maximum of the univariate 
                                   t-distributed test statistics and compares it to critical values tabulated by Dunnett. The adjusted p-values are calculated from the 
                                   multivariate \\(t\\)-distribution of the vector \\(t\\) of test statistics \\(t_1,\\ldots,t_m\\), thereby accounting for correlations between 
                                   these statistics."),
           "Step-Down Dunnett" = withMathJax("Step-wise extension of Dunnett's test for comparing multiple treatment groups against a single control. 
                                              The null hypotheses \\( H_{(1)}, \\ldots, H_{(m)} \\) (corresponding to the ordered p-values 
                                               \\( p_{(1)} \\leq \\ldots \\leq p_{(m)} \\)) are tested sequentially, 
                                              and the adjusted p-values \\( q_{(i)} = \\max_{I : i \\in I} p_I \\) are obtained, 
                                              stopping as soon as \\( q_{(i)} > \\alpha \\) 
                                              (The index set \\( I \\) refers to intersections of elementary hypotheses in the closed testing procedure that include \\( i \\)). 
                                              The use of max-\\(t\\) tests allows for shortcuts to the full closure. 
                                              Therefore, power is higher than for Dunnett's test."),
           "Tukey" = withMathJax("Standard procedure for all pairwise comparison problems. For balanced groups and homoscedasticity, pairwise t-tests are 
                                 computed: \\(t_{i,j} = \\frac{\\bar{y_i}-\\bar{y_j}}{s \\sqrt{\\frac{2}{n}}}\\), \\(i,j \\in M\\), with \\(i\\neq j\\) for groups of size \\(n\\). 
                                 The maximum over the absolute values of these statistics is taken and compared to critical values from the studentised range distribution 
                                 which accounts for correlations between the statistics."),
           "S2" = withMathJax("A stepwise extension of the Bonferroni test.
                                    Due to the closed nature of the procedure, the FWER is maintained. 
                                    It makes use of logical interrelationships between the hypotheses, thereby effectively reducing the size 
                                    of the family of hypotheses to be tested in successive steps. Let \\(t_i\\) denote the family size in step \\(i\\). 
                                    It is determined based on which specific hypotheses have been rejected at steps \\(1,\\ldots,i-1\\).
                                    The procedure starts by testing the first null hypothesis analogously to Bonferroni and Holm: it is rejected 
                                    if \\(p_{(1)} \\leq \\frac{\\alpha}{m}\\). Having rejected \\(H_{(1)}\\), the cardinality of the largest family of 
                                    hypotheses that can simultaneously be true (\\(t_i\\)) is used as the Bonferroni multiplier in the next step. 
                                    If the free combination condition is satisfied, the S2 procedure reduces to the Holm procedure."),
           "Step-Down Tukey" = "A stepwise extension of the Tukey test using truncated closure. Relevant subsets are found by using
           rank conditions involving contrast matrices for the various subsets. The set of subsets considered
           in computing the adjusted p-values can often be reduced due to the restricted combination condition.
           The procedure is uniformly more powerful than Tukey and S2."
    )
  }
  
  # activate info icon/button
  observe({
    all_choices <- c("Bonferroni", "Holm", "Hochberg", "Benjamini-Hochberg", "Hommel", "Dunnett", 
                     "Step-Down Dunnett", "Tukey", "S2", "Step-Down Tukey")
    
    lapply(all_choices, function(proc) {
      observeEvent(input[[paste0("info_", proc)]], {
        showModal(modalDialog(
          title = paste(proc, "Information"),
          procedure_info_text(proc),
          easyClose = TRUE,
          footer = NULL
        ))
      })
    })
  })
  
  # UI names -> R names
  method_mapping <- c(
    "Bonferroni" = "bonferroni",
    "Holm" = "holm",
    "Hochberg" = "hochberg",
    "Benjamini-Hochberg" = "BH",
    "Hommel" = "hommel",
    "Dunnett" = "single-step", 
    "Tukey" = "single-step",
    "S2" = "Shaffer",
    "Step-Down Tukey" = "Westfall",
    "Step-Down Dunnett" = "free" 
  )

  
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


  ##############################################################################################
  find_best_comparison <- function(samples, type_of_comparison, custom_comparisons = NULL, multiple_test = NULL, dsg){
    
    # may be specific to the multiple testing procedure being selected  
    
    # how should the output be structured? 
    # 1. we need the group names for any further data usage
    # 2. we should store the results for standard design, to avoid doing the same computations again

   
    if (length(multiple_test) > 1) {
      results_list <- lapply(multiple_test, function(mt) {
        find_best_comparison(samples = samples, type_of_comparison = type_of_comparison, custom_comparisons = custom_comparisons, multiple_test = mt)
      })
      names(results_list) <- multiple_test
      return(results_list)
    }
    
    if((type_of_comparison == "all_pairwise" && is.null(multiple_test))){
      best_comparison <- samples %>% pivot_longer(., cols=everything(), names_to = "group") %>% 
        compare_means(value ~ group, data=., method = "t.test") %>% dplyr::select(group1, group2, p) %>% slice(which.min(p))
      best_groups <- as.numeric(gsub(pattern = "V", replacement = "", x = c(best_comparison$group1, best_comparison$group2)))

      # effect, p-value and CI boundaries will not be taken from here 
      best_pvalue_glht <- NULL
      effect_est_glht <- NULL
      CI_glht_lwr <- NULL
      CI_glht_upr <- NULL
      
      results <- best_comparison
    } 
    
    else if (!is.null(type_of_comparison) && type_of_comparison == "all_pairwise" && !is.null(multiple_test) && multiple_test %in% names(method_mapping)) {
      glht_fit <- glht(aov(value ~ name, 
                           data = samples %>% pivot_longer(cols = everything()) %>% mutate(name = as.factor(name))),
                       linfct = mcp(name = "Tukey"))
      r_method <- method_mapping[multiple_test]
      if (multiple_test == "single-step"){
        glht_fit_summary <- summary(glht_fit, test = adjusted(type = r_method), maxpts = 1e6)
      } else {
        glht_fit_summary <- summary(glht_fit, test = adjusted(type = r_method))
      }
      # need to get the correct indices of the best comparison from glht-object
      best_comparison_name <- names(glht_fit_summary$test$tstat)[which.min(glht_fit_summary$test$pvalues)]
      best_groups <- as.numeric(gsub(pattern = "V", replacement = "", x = unlist(strsplit(best_comparison_name, split = "-"))))
      best_pvalue_glht <- as.numeric(glht_fit_summary$test$pvalue[which.min(glht_fit_summary$test$pvalue)])
      
      # confidence interval
      best_index <- which.min(glht_fit_summary$test$pvalues)
      CI_all <- confint(glht_fit_summary)
      CI_glht <- CI_all$confint[best_index,]
      CI_glht_lwr <- CI_glht[2]
      CI_glht_upr <- CI_glht[3]
      
      # effect estimate
      effect_est_glht <- as.numeric(CI_glht[1])
      
      results <- glht_fit_summary
    }
    
    else if(type_of_comparison == "all_vs_control" && is.null(multiple_test)){
      glht_fit <- glht(aov(value ~ name, 
                           data = samples %>% pivot_longer(cols=everything()) %>% mutate(name = as.factor(name))), 
                       linfct = mcp(name = "Dunnett"))
      glht_fit_summary <- summary(glht_fit, test = adjusted(type = "none")) # no multiplicity adjustment
      # need to get the correct indices of the best comparison from glht-object
      best_comparison_name <- names(glht_fit_summary$test$tstat)[which.min(glht_fit_summary$test$pvalues)]
      best_groups <- as.numeric(gsub(pattern = "V", replacement = "", x = unlist(strsplit(best_comparison_name, split = "-"))))

      # effect, p-value and CI boundaries will not be taken from here 
      best_pvalue_glht <- NULL
      effect_est_glht <- NULL
      CI_glht_lwr <- NULL
      CI_glht_upr <- NULL
      
      results <- glht_fit_summary
    } 
    
    else if (!is.null(type_of_comparison) && type_of_comparison == "all_vs_control" && !is.null(multiple_test) && multiple_test %in% names(method_mapping)) {
      glht_fit <- glht(aov(value ~ name,
                           data = samples %>% pivot_longer(cols = everything()) %>% mutate(name = as.factor(name))),
                       linfct = mcp(name = "Dunnett"))
      r_method <- method_mapping[multiple_test]
      if (multiple_test == "single-step"){
        glht_fit_summary <- summary(glht_fit, test = adjusted(type = r_method), maxpts = 1e6)
      } else {
        glht_fit_summary <- summary(glht_fit, test = adjusted(type = r_method))
      }
      # need to get the correct indices of the best comparison from glht-object
      best_comparison_name <- names(glht_fit_summary$test$tstat)[which.min(glht_fit_summary$test$pvalues)]
      best_groups <- as.numeric(gsub(pattern = "V", replacement = "", x = unlist(strsplit(best_comparison_name, split = "-"))))
      best_pvalue_glht <- as.numeric(glht_fit_summary$test$pvalue[which.min(glht_fit_summary$test$pvalue)])
      
      # confidence interval 
      best_index <- which.min(glht_fit_summary$test$pvalues)
      CI_all <- confint(glht_fit_summary)
      CI_glht <- CI_all$confint[best_index,]
      CI_glht_lwr <- CI_glht[2]
      CI_glht_upr <- CI_glht[3]
      
      # effect estimate 
      effect_est_glht <- as.numeric(CI_glht[1])
      
      results <- glht_fit_summary
    }
    
    else if(type_of_comparison == "custom" && is.null(multiple_test)){# custom comparison
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
      
      # effect, p-value and CI boundaries will not be taken from here 
      best_pvalue_glht <- NULL
      effect_est_glht <- NULL
      CI_glht_lwr <- NULL
      CI_glht_upr <- NULL
      
      results <- all_contrasts
    }
    
    else if(!is.null(type_of_comparison) && type_of_comparison == "custom" && !is.null(multiple_test) && multiple_test %in% names(method_mapping)) {
      custom_comparisons_glht <- gsub("(\\d+)", "V\\1", gsub(" vs\\. ", " - ", custom_comparisons))
      no_groups <- ncol(samples)
      contrast_matrix <- lapply(custom_comparisons_glht, function(contrast) {
        parts <- as.numeric(gsub("V", "", unlist(strsplit(contrast, "-"))))
        vec <- rep(0, no_groups)
        vec[parts[1]] <- 1
        vec[parts[2]] <- -1 
        return(vec)
      }) %>% do.call(rbind,.)
      
      
      rownames(contrast_matrix) <- custom_comparisons_glht
      glht_fit <- glht(aov(value ~ name,
                           data = samples %>% pivot_longer(cols = everything()) %>% mutate(name = as.factor(name))),
                           linfct = contrast_matrix)
      r_method <- method_mapping[multiple_test]
      if (multiple_test == "single-step"){
        glht_fit_summary <- summary(glht_fit, test = adjusted(type = r_method), maxpts = 1e6)
      } else {
        glht_fit_summary <- summary(glht_fit, test = adjusted(type = r_method))
      }
      # need to get the correct indices of the best comparison from glht-object
      best_comparison_name <- names(glht_fit_summary$test$tstat)[which.min(glht_fit_summary$test$pvalues)]
      best_groups <- as.numeric(gsub("V", "", unlist(strsplit(best_comparison_name, "-"))))
      best_pvalue_glht <- as.numeric(glht_fit_summary$test$pvalue[which.min(glht_fit_summary$test$pvalue)])
      
      # confidence interval 
      best_index <- which.min(glht_fit_summary$test$pvalues)
      CI_all <- confint(glht_fit_summary)
      CI_glht <- CI_all$confint[best_index,]
      CI_glht_lwr <- CI_glht[2]
      CI_glht_upr <- CI_glht[3]

      # effect estimate 
      effect_est_glht <- as.numeric(CI_glht[1]) 
      
      results <- glht_fit_summary
    } 
    
    return(list("best_groups"=best_groups, "results"=results, "best_pvalue_glht"=best_pvalue_glht, "effect_est_glht"=effect_est_glht, "CI_glht_lwr"=CI_glht_lwr, "CI_glht_upr"=CI_glht_upr))
  }
  
  ##############################################################################################
  estimate_effect <- function(data_stage1, data_stage2 = NULL, n_per_group, dsg, type_of_comparison, multiple_test = NULL, stanmodel){
    # specific to each design
    # input: 
    # data of the best comparison, may include 2nd stage data
    # type of comparison
    # multiple testing correction, if applicable
    
    if((dsg == "standard" & is.null(multiple_test))){
      fit <- t.test(data_stage1[,1], data_stage1[,2])
      effect_est <- as.numeric(diff(fit$estimate)*(-1)) # effect estimate is calculated as mean(x)-mean(y)
      CI <- fit$conf.int[1:2]
      best_pvalue <- fit$p.value  
    }
    if((dsg == "repeat_and_pool" & is.null(multiple_test))){
      data <- rbind(data_stage1, data_stage2)
      fit <- t.test(data[,1], data[,2])
      effect_est <- as.numeric(diff(fit$estimate)*(-1))
      CI <- fit$conf.int[1:2]
      best_pvalue <- fit$p.value  
    }
    if(dsg == "repeat_and_replace"){
      fit <- t.test(data_stage2[,1], data_stage2[,2])
      effect_est <- as.numeric(diff(fit$estimate)*(-1))
      CI <- fit$conf.int[1:2]
      best_pvalue <- fit$p.value  
    }
    if(dsg == "bayes"){
      historical_data <- bind_rows(data.frame(y = data_stage1[,1], x = 0), data.frame(y = data_stage1[,2], x = 1))
      historical_fit <- lm(y ~ x, data = historical_data)
      fit <- rstan::sampling(stanmodel, 
                             data = list("df_beta" = historical_fit$df.residual, "y" = c(data_stage2[,1], data_stage2[,2]), 
                                         "group" = c(rep(0,length(data_stage2[,1])), rep(1,length(data_stage2[,1]))), "N" = length(c(data_stage2[,1], data_stage2[,2])), 
                                         "delta" = input$mixing_weight, "diff_prior" = as.numeric(historical_fit$coefficients[2]), "sigma_prior" = summary(historical_fit)$coefficients[2,2], 
                                         "diff_skeptical" = 0, "sigma_skeptical" = summary(historical_fit)$coefficients[2,2]), 
                             refresh = 0)      
      posterior_mixture <- rstan::extract(fit) %>% as_tibble()
      effect_est <- mean(posterior_mixture$beta)
      CI <- as.numeric(quantile(posterior_mixture$beta, probs = c(0.025, 0.975)))
      best_pvalue <- min(sum(posterior_mixture$beta < 0)/nrow(posterior_mixture), sum(posterior_mixture$beta > 0)/nrow(posterior_mixture))*2
    }
  
    return(list("effect_est" = effect_est, "CI_lwr" = CI[1], "CI_upr" = CI[2], "best_pvalue" = best_pvalue, "dsg" = dsg))
  }
  
  ##############################################################################################
  ###### wrap all functions together into a function to perform the simulations
  run_simulation <- function(n_per_group, no_groups, distribution, parameters_new, custom_comparisons = NULL,
                             type_of_comparison, dsg, stanmodel = NULL, multiple_test = NULL){
    
    ### initial data sampling based on user defined input parameters 
    samples <- sample_data(distribution = distribution, parameters_new = parameters_new, no_groups, n_per_group) 
    
    ### identify best comparison
    results <- find_best_comparison(samples = samples, type_of_comparison = type_of_comparison, custom_comparisons = custom_comparisons, multiple_test = multiple_test)
    #browser()
    # at least two multiple adjustments were chosen
    if (length(multiple_test) > 1) {
      effect_results_list <- lapply(names(results), function(mt) {
        # save best p-value (first stage)
        best_pvalue_glht <- results[[mt]]$best_pvalue_glht
        # sample data for second stage
        new_data_2stage <- sample_data(distribution = distribution, parameters_new = parameters_new[results[[mt]]$best_groups], 
                                       no_groups = 2, n_per_group = n_per_group)
        colnames(new_data_2stage) <- colnames(samples[, results[[mt]]$best_groups])
        
        
        effect_results <- estimate_effect(data_stage1 = samples[, results[[mt]]$best_groups], data_stage2 = new_data_2stage, n_per_group = n_per_group,
                                          dsg = dsg, type_of_comparison = type_of_comparison, stanmodel = stanmodel)
        
        effect_results$multiple_test <- mt
        
        effect_results$effect_est_glht <- results[[mt]]$effect_est_glht
        
        effect_results$CI_glht_lwr <- results[[mt]]$CI_glht_lwr
        
        effect_results$CI_glht_upr <- results[[mt]]$CI_glht_upr
        
        effect_results$best_pvalue_glht <- best_pvalue_glht
        
        effect_results$best_groups_1 <- results[[mt]]$best_groups[1]
        
        effect_results$best_groups_2 <- results[[mt]]$best_groups[2]
        
        effect_results$no_groups <- no_groups
        
        effect_results$n_per_group <- n_per_group
        
        return(effect_results)
      })

      return(do.call(rbind, lapply(effect_results_list, as.data.frame)))
    } else {
      # less than two multiple adjustments were chosen 
      # save best p-value (glht model)
      best_pvalue_glht <- results$best_pvalue_glht
      # sample data for second stage
      new_data_2stage <- sample_data(distribution = distribution, parameters_new = parameters_new[results$best_groups], 
                                     no_groups = 2, n_per_group = n_per_group)
      colnames(new_data_2stage) <- colnames(samples[, results$best_groups])
      
      effect_results <- estimate_effect(data_stage1 = samples[, results$best_groups], data_stage2 = new_data_2stage, n_per_group = n_per_group,
                                        dsg = dsg, type_of_comparison = type_of_comparison, stanmodel = stanmodel)
      
      effect_results$multiple_test <- ifelse(is.null(multiple_test), "none", multiple_test)
      
      effect_results$effect_est_glht <- results$effect_est_glht
      
      effect_results$CI_glht_lwr <- results$CI_glht_lwr
      
      effect_results$CI_glht_upr <- results$CI_glht_upr
      
      effect_results$best_pvalue_glht <- best_pvalue_glht
      
      effect_results$best_groups_1 <- results$best_groups[1]
      
      effect_results$best_groups_2 <- results$best_groups[2]
      
      effect_results$no_groups <- no_groups
      
      effect_results$n_per_group <- n_per_group

      return(as.data.frame(effect_results))
    }
    
    
    # returns dataframe with 7 columns: 
    # 1. effect_est 2. CI_lwr 3. CI_upr 4. best_pvalue 5. dsg  6. multiple_test 7. effect_est_glht 8. CI_glht_lwr 9. CI_glht_upr 10. best_pvalue_glht
    
    # if multiple_test = NULL:
    # columns 1,2,3,4 are used for the plots 
    # if multiple_test != NULL: 
    # columns 7,8,9,10 are used for the plots
    
  }
  
  ##############################################################################################  
  ###### function to perform the simulations for each design
  
  run_simulation_for_each_dsg <- function(no_groups, parameters_new, n_per_group, distribution, type_of_comparison,
                                          custom_comparisons, stanmodel, choose_bayes, multiple_test){
    
    if (is.null(multiple_test) == TRUE){
      
      # standard design
      results_standard <- run_simulation(n_per_group = n_per_group, no_groups = no_groups, distribution = distribution, parameters_new = parameters_new, 
                                         custom_comparisons = custom_comparisons, dsg = "standard", type_of_comparison = type_of_comparison, 
                                         stanmodel = stanmodel, multiple_test = NULL) %>% as.data.frame
      
      results_standard_adj <- NULL
      
    } else {
    
    # results when adjusting for multiple comparisons (-> standard + adjustment)
    results_standard_adj <- run_simulation(n_per_group = n_per_group, no_groups = no_groups, distribution = distribution, parameters_new = parameters_new, 
                                       custom_comparisons = custom_comparisons, dsg = "standard", type_of_comparison = type_of_comparison, 
                                       stanmodel = stanmodel, multiple_test = multiple_test) %>% as.data.frame
    
    # standard design without multiple correction 
    results_standard <- run_simulation(n_per_group = n_per_group, no_groups = no_groups, distribution = distribution, parameters_new = parameters_new, 
                                       custom_comparisons = custom_comparisons, dsg = "standard", type_of_comparison = type_of_comparison, 
                                       stanmodel = stanmodel, multiple_test = NULL) %>% as.data.frame
    }
    
    results_repeat_and_pool <- run_simulation(n_per_group = n_per_group, no_groups = no_groups, distribution = distribution, parameters_new = parameters_new, 
                                              custom_comparisons = custom_comparisons, dsg = "repeat_and_pool", type_of_comparison = type_of_comparison, 
                                              stanmodel = stanmodel, multiple_test = NULL) %>% as.data.frame
    
    results_repeat_and_replace <- run_simulation(n_per_group = n_per_group, no_groups = no_groups, distribution = distribution, parameters_new = parameters_new, 
                                                 custom_comparisons = custom_comparisons, dsg = "repeat_and_replace", type_of_comparison = type_of_comparison, 
                                                 stanmodel = stanmodel, multiple_test = NULL) %>% as.data.frame
    
    if(choose_bayes == "Yes"){
      results_bayes <- run_simulation(n_per_group = n_per_group, no_groups = no_groups, distribution = distribution, parameters_new = parameters_new, 
                                      custom_comparisons = custom_comparisons, dsg = "bayes", type_of_comparison = type_of_comparison, 
                                      stanmodel = stanmodel, multiple_test = NULL) %>% as.data.frame
    } else{results_bayes <- data.frame(NULL)}
    
    return(bind_rows(results_standard_adj,results_standard, results_repeat_and_pool, results_repeat_and_replace, results_bayes))
  }
  

  #########################################################################################################   
  # Reactive expression to perform the simulation -------------------
  gen_data <- eventReactive(input$run, {
    
    if(input$choose_bayes == "Yes"){
      stanmodel <- stan_model("mixture_model_ttest.stan")
    } else{stanmodel <- ""}
    
    # the distribution can be provided by two different inputs, so we need a variable containing this information to pass it further into functions
    if(any(c(input$choose_distribution_custom == "normal", input$choose_distribution == "normal"))){
      distribution = "normal" 
    } else if(any(c(input$choose_distribution_custom == "lognormal", input$choose_distribution == "lognormal"))){
      distribution = "lognormal" 
    } else if(any(c(input$choose_distribution_custom == "cauchy", input$choose_distribution == "cauchy"))){
      distribution = "cauchy" 
    } 
    
    # the parameters specified need to be transformed to efficiently sample the data
    # specifically, a list containing two vectors of length k like this: list(c(0,1,2,3,4), c(1,2,3,4,5)) 
    # should become a list which contains k vectors of length 2, like this: list(c(0,1), c(1,2), c(2,3), c(3,4), c(4,5))
    # if no parameters were specified but instead an effect size for all comparisons to control was defined,
    # these need also be defined in the correct form, i.e. standard normal for control and means = effect_size for treatments
    if(input$define_comparisons == "all_vs_control" & !is.null(input$define_effect) && input$define_effect == "all_equal"){
      parameters_new <- mapply(c, c(0, rep(input$choose_effect, input$no_groups-1)), rep(1, input$no_groups), SIMPLIFY = FALSE)
      parameters_new_tls <- list()
      
      # if traffic light system was chosen, parameters for no_groups, tls_no_groups, tls_reference_green, tls_reference_yellow have to be saved
       if (input$traffic_light_system == "Yes"){
         
         # collect all input parameters 
         inputs_no_groups_tls <- c(as.numeric(input$no_groups), as.numeric(input$tls_no_groups),as.numeric(input$tls_reference_green), as.numeric(input$tls_reference_yellow))
         inputs_no_groups_tls <- unique(inputs_no_groups_tls) # keep only unique values
         
         # generate parameters for other scenarios to appear in traffic light system 
         parameters_new_tls <- lapply(inputs_no_groups_tls, function(i) {
           mapply(c, c(rep(input$choose_effect, i)), rep(1, i), SIMPLIFY = FALSE)})
         }
                                                
      
    } else {
      parameters_new <- mapply(c, 
                               as.numeric(unlist(strsplit(x = input$user_def_par1, split = ","))),
                               as.numeric(unlist(strsplit(x = input$user_def_par2, split = ","))), 
                               SIMPLIFY = FALSE)
      parameters_new_tls <- list()
      
      # if traffic light system was chosen, parameters for no_groups, no_groups_tls2 & no_groups_tls3 have to be computed & saved
        if (input$traffic_light_system == "Yes"){
          
          # collect all input parameters
          inputs_no_groups_tls <- c(as.numeric(input$no_groups), as.numeric(input$tls_no_groups),as.numeric(input$tls_reference_green), as.numeric(input$tls_reference_yellow))
          inputs_no_groups_tls <- unique(inputs_no_groups_tls) # keep only unique values
          
          parameters_all <- mapply(c, 
                               as.numeric(unlist(strsplit(x = input$user_def_par1, split = ","))),
                               as.numeric(unlist(strsplit(x = input$user_def_par2, split = ","))), 
                               SIMPLIFY = FALSE)
          
          parameters <- parameters_all[[1]] # since tls is only valid under null hypothesis: take parameters of group 1 and replicate for no_groups inputs
          
          # generate parameters for other scenarios to appear in traffic light system 
          parameters_new_tls <- lapply(inputs_no_groups_tls, function(i) {
            mapply(c, c(rep(parameters[[1]], i)), rep(parameters[[2]], i), SIMPLIFY = FALSE)})
        }
    }
          
    # custom comparisons instead of all comparisons (input for custom_comparisons in run_simulation_for_each_dsg)
    if(input$define_comparisons == "custom") {cust_combs <- input$custom_comparisons} else {cust_combs <- NULL}  
    
    
    ### collect all inputs together and sample the data
    
    # input <- list("no_groups" = 10, "n_per_group" = 10, "choose_bayes" = "Yes", "define_comparisons" = "all_pairwise", "choose_distribution" = "normal",
    #               "user_def_par1" = rep("0", 10), "user_def_par2" = rep("1", 10), "nsim" = 10)
    # only required for test runs

   
    multiple_test <- if(input$multiple_test == "No") NULL else selected_methods()

    
    
    ### results: either with or without data for traffic light system
    
    simulation_results <- lapply(1:input$nsim, function(i) {
      
      run_simulation_for_each_dsg(no_groups = input$no_groups, parameters_new = parameters_new, n_per_group = input$n_per_group,
                                  distribution = distribution, type_of_comparison = input$define_comparisons,
                                  custom_comparisons = cust_combs, stanmodel = stanmodel,
                                  choose_bayes = input$choose_bayes, multiple_test = multiple_test)
      
    }) %>% do.call(rbind,.) %>% as_tibble()
    
    
    if (input$traffic_light_system == "Yes"){
      
      if (input$define_comparisons != "custom"){
        
      n_per_group_tls <- unique(c(as.numeric(input$n_per_group), as.numeric(input$tls_n_per_group)))
      no_groups_tls <- unique(c(as.numeric(input$no_groups), as.numeric(input$tls_no_groups), as.numeric(input$tls_reference_green), as.numeric(input$tls_reference_yellow)))
      
      tls_grid <- expand.grid(
        n_per_group = n_per_group_tls,
        no_groups = no_groups_tls
      )
      
      tls_simulation_results <- foreach(i = 1:input$nsim, .combine = 'rbind') %do% {
        foreach(t = 1:nrow(tls_grid), .combine = 'rbind') %do% {
          run_simulation_for_each_dsg(no_groups = tls_grid$no_groups[t], parameters_new = parameters_new_tls[[which(lengths(parameters_new_tls) == tls_grid$no_groups[[t]])]], 
                                      n_per_group = tls_grid$n_per_group[t],
                                      distribution = distribution, type_of_comparison = input$define_comparisons,
                                      custom_comparisons = cust_combs, stanmodel = stanmodel,
                                      choose_bayes = input$choose_bayes, multiple_test = multiple_test)
        }
      } %>% as_tibble()
      
      
      
      # for custom comparisons: all pairwise comparisons will be displayed in the traffic light system
      # -> two simulations to get results for custom comparisons and for pairwise comparisons
      } else if (input$define_comparisons == "custom") {
        
        n_per_group_tls_pw <- unique(c(as.numeric(input$tls_n_per_group), as.numeric(input$n_per_group)))
        no_groups_tls_pw <- unique(c(as.numeric(input$tls_no_groups), as.numeric(input$tls_reference_green), as.numeric(input$tls_reference_yellow)))
        
        tls_grid_pw <- expand.grid(
          n_per_group = n_per_group_tls_pw,
          no_groups = no_groups_tls_pw
        )
        
        tls_simulation_results_pw <- foreach(i = 1:input$nsim, .combine = 'rbind') %do% {
          foreach(t = 1:nrow(tls_grid_pw), .combine = 'rbind') %do% {
            run_simulation_for_each_dsg(no_groups = tls_grid_pw$no_groups[t], parameters_new = parameters_new_tls[[which(lengths(parameters_new_tls) == tls_grid_pw$no_groups[[t]])]], 
                                        n_per_group = tls_grid_pw$n_per_group[t],
                                        distribution = distribution, type_of_comparison = "all_pairwise",
                                        custom_comparisons = cust_combs, stanmodel = stanmodel,
                                        choose_bayes = input$choose_bayes, multiple_test = multiple_test)
          }
        } %>% as_tibble()
        
        n_per_group_cust <- as.numeric(input$n_per_group)#unique(c(as.numeric(input$n_per_group), as.numeric(input$tls_n_per_group)))
        no_groups_cust <- as.numeric(input$no_groups)
        
        tls_grid_cust <- expand.grid(
          n_per_group = n_per_group_cust,
          no_groups = no_groups_cust
        )
        
        tls_simulation_results_cust <- foreach(i = 1:input$nsim, .combine = 'rbind') %do% {
          foreach(t = 1:nrow(tls_grid_cust), .combine = 'rbind') %do% {
            run_simulation_for_each_dsg(no_groups = tls_grid_cust$no_groups[t], parameters_new = parameters_new_tls[[which(lengths(parameters_new_tls) == tls_grid_cust$no_groups[[t]])]], 
                                        n_per_group = tls_grid_cust$n_per_group[t],
                                        distribution = distribution, type_of_comparison = "custom",
                                        custom_comparisons = cust_combs, stanmodel = stanmodel,
                                        choose_bayes = input$choose_bayes, multiple_test = multiple_test)
          }
      } %>% as_tibble()
        
        tls_simulation_results <- bind_rows(tls_simulation_results_cust, tls_simulation_results_pw)
      }
    }
      
    
    ## preparation for traffic light plots 
    
    if (input$traffic_light_system == "Yes") {
    
    traffic_light_system_mse <- tls_simulation_results %>%
      filter(multiple_test == "none") %>%
      rowwise() %>%
      mutate(
        true_effect = parameters_new_tls[[which(lengths(parameters_new_tls) == no_groups)]][[best_groups_1]][1] - parameters_new_tls[[which(lengths(parameters_new_tls) == no_groups)]][[best_groups_2]][1],
        squared_error = diff(c(true_effect, effect_est))^2,
        bias = diff(c(true_effect, effect_est))) %>%
      ungroup() %>%
      group_by(dsg, no_groups, n_per_group) %>%
      summarise(MSE = mean(squared_error),
                bias = mean(bias),
                Var_est = var(effect_est), .groups = "keep") %>%
      group_by(n_per_group) %>%
      mutate(MSE_ref = MSE / MSE[dsg == "standard" & no_groups == as.numeric(input$tls_reference_green)]) %>%
      mutate(colors = case_when(MSE_ref <= MSE_ref[dsg == "standard" & no_groups == as.numeric(input$tls_reference_green)]*1.05 ~ "green",
                                MSE_ref <= MSE_ref[dsg == "standard" & no_groups == as.numeric(input$tls_reference_yellow)]*1.05 ~ "yellow",
                                TRUE ~ "red"))  
    } else { 
      traffic_light_system_mse = NULL}

    
    # compute rejection rate
    
    rejection_rates <- simulation_results %>%
     mutate(pvalue_used = ifelse(multiple_test == "none", best_pvalue, best_pvalue_glht),
            method = ifelse(multiple_test == "none", "unadjusted", multiple_test)) %>%
      group_by(dsg, method) %>% 
      summarise(rejection_rate = mean(pvalue_used < 0.05), .groups = "drop") %>% 
      ungroup()
    
    # prepare data for plots per adjustment method
    # adjustments only for standard design
    results_standard_adj <- simulation_results %>% filter(dsg == "standard" & multiple_test != "none") 
    
    
    return(list(simulation_results = simulation_results, rejection_rates = rejection_rates, results_standard_adj = results_standard_adj, traffic_light_system_mse = traffic_light_system_mse))
  })
  ###########################################################
    

  # Plots -------------------
  output$plot <- renderPlot({

  data <- gen_data()$simulation_results
  data_tls_mse <- gen_data()$traffic_light_system_mse
  results_standard_adj <- gen_data()$results_standard_adj
  rejection_rates <- gen_data()$rejection_rates
  
  
  

  p1.1 <- data %>% filter(dsg == "standard" & multiple_test == "none") %>% ggplot(.) + geom_histogram(aes(x=effect_est), binwidth = 0.16) + ggtitle("Standard design") +
    xlab("Observed effect estimate with smallest p-value") + ylab("Number of simulations") +
    coord_cartesian(xlim= c(-4,4)) + #scale_x_continuous(limits = c(-4,4)) + 
    theme_bw()
  
  p2.1 <- data %>% filter(dsg == "standard" & multiple_test == "none") %>% arrange(CI_lwr) %>% rowid_to_column("sim") %>%
    mutate(CI_excludes_0 = case_when(CI_lwr < 0 & CI_upr > 0 ~ FALSE, TRUE ~ TRUE)) %>%
    mutate(CI_excludes_0 = factor(CI_excludes_0, levels = c(TRUE, FALSE))) %>%
    ggplot(.) + geom_errorbar(aes(y = sim, xmin = CI_lwr, xmax = CI_upr, colour = CI_excludes_0)) +
    scale_colour_manual(values = c("red", "grey50"), breaks = c(TRUE, FALSE), guide = "none") +
    ggtitle("Standard design")  +
    xlab("Effect size") + theme_bw() +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  
  p3.1 <- data %>% filter(dsg == "standard" & multiple_test == "none") %>% ggplot(.) + geom_histogram(aes(x=best_pvalue), bins = 50, boundary = 0) + ggtitle("Standard design") +
    xlab("Smallest observed p-value") + ylab("Number of simulations") + theme_bw() + scale_x_continuous(limits = c(0,1))
  
  p4.1 <- rejection_rates %>% filter(dsg == "standard" & method == "unadjusted") %>% ggplot(., aes(x = method, y = rejection_rate, fill = method)) +
    geom_col() + 
    geom_text(aes(label = round(rejection_rate, 3)), vjust = -0.5) +
    xlab("Method") + ylab("Rejection Rate") + ggtitle("Rejection Rate") + theme_bw() + coord_cartesian(ylim = c(0,1))
  
  p1.2 <- data %>% filter(dsg == "repeat_and_pool") %>% ggplot() + 
    geom_histogram(aes(x=effect_est), binwidth = 0.16 ) + ggtitle("Add follow-up") +
    xlab("Observed effect estimate with smallest p-value") + ylab("Number of simulations") + 
    coord_cartesian(xlim = c(-4,4)) +#scale_x_continuous(limits = c(-4,4)) + 
    theme_bw()
  
  p2.2 <- data %>% filter(dsg == "repeat_and_pool") %>% arrange(CI_lwr) %>% rowid_to_column("sim") %>% 
    mutate(CI_excludes_0 = case_when(CI_lwr < 0 & CI_upr > 0 ~ FALSE, TRUE ~ TRUE)) %>% 
    mutate(CI_excludes_0 = factor(CI_excludes_0, levels = c(TRUE, FALSE))) %>% 
    ggplot(.) + geom_errorbar(aes(y = sim, xmin = CI_lwr, xmax = CI_upr, colour = CI_excludes_0)) + 
    scale_colour_manual(values = c("red", "grey50"), breaks = c(TRUE, FALSE), guide = "none") + 
    ggtitle("Add follow-up")  +
    xlab("Effect size") + theme_bw() + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  
  p3.2 <- data %>% filter(dsg == "repeat_and_pool") %>% ggplot(.) + geom_histogram(aes(x=best_pvalue), bins = 50, boundary = 0) + ggtitle("Add follow-up") +
    xlab("Smallest observed p-value") + ylab("Number of simulations") + theme_bw() + scale_x_continuous(limits = c(0,1))
  
  p4.2 <- rejection_rates %>% filter(dsg == "repeat_and_pool") %>% ggplot(., aes(x = method, y = rejection_rate, fill = method)) +
    geom_col() + 
    geom_text(aes(label = round(rejection_rate, 3)), vjust = -0.5) +
    xlab("Method") + ylab("Rejection Rate") + ggtitle("Rejection Rate") + theme_bw() + coord_cartesian(ylim = c(0,1))
  
  p1.3 <- data %>% filter(dsg == "repeat_and_replace") %>% ggplot() + 
    geom_histogram(aes(x=effect_est), binwidth = 0.16) + ggtitle("Repeat and replace") +
    xlab("Observed effect estimate with smallest p-value") + ylab("Number of simulations") + 
    coord_cartesian(xlim = c(-4,4)) + #scale_x_continuous(limits = c(-4,4)) + 
    theme_bw()
  
  p2.3 <- data %>% filter(dsg == "repeat_and_replace") %>% arrange(CI_lwr) %>% rowid_to_column("sim") %>% 
    mutate(CI_excludes_0 = case_when(CI_lwr < 0 & CI_upr > 0 ~ FALSE, TRUE ~ TRUE)) %>% 
    mutate(CI_excludes_0 = factor(CI_excludes_0, levels = c(TRUE, FALSE))) %>% 
    ggplot(.) + geom_errorbar(aes(y = sim, xmin = CI_lwr, xmax = CI_upr, colour = CI_excludes_0)) + 
    scale_colour_manual(values = c("red", "grey50"), breaks = c(TRUE, FALSE), guide = "none") + 
    ggtitle("Repeat and replace")  +
    xlab("Effect size") + theme_bw() + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  
  p3.3 <- data %>% filter(dsg == "repeat_and_replace") %>% ggplot(.) + geom_histogram(aes(x=best_pvalue), bins = 50, boundary = 0) + ggtitle("Repeat and replace") +
    xlab("Smallest observed p-value") + ylab("Number of simulations") + theme_bw() + scale_x_continuous(limits = c(0,1))
  
  p4.3 <- rejection_rates %>% filter(dsg == "repeat_and_replace") %>% ggplot(., aes(x = method, y = rejection_rate, fill = method)) +
    geom_col() + 
    geom_text(aes(label = round(rejection_rate, 3)), vjust = -0.5) +
    xlab("Method") + ylab("Rejection Rate") + ggtitle("Rejection Rate") + theme_bw() + coord_cartesian(ylim = c(0,1))
  
  if(input$choose_bayes == "Yes"){
    p1.4 <- data %>% filter(dsg == "bayes") %>% ggplot() +
      geom_histogram(aes(x=effect_est), binwidth = 0.16) + ggtitle("Robust mixture prior") +
      xlab("Observed effect estimate with smallest p-value") + ylab("Number of simulations") +
      coord_cartesian(xlim = c(-4,4)) + #scale_x_continuous(limits = c(-4,4)) + 
      theme_bw()
    
    p2.4 <- data %>% filter(dsg == "bayes") %>% arrange(CI_lwr) %>% rowid_to_column("sim") %>%
      mutate(CI_excludes_0 = case_when(CI_lwr < 0 & CI_upr > 0 ~ FALSE, TRUE ~ TRUE)) %>%
      mutate(CI_excludes_0 = factor(CI_excludes_0, levels = c(TRUE, FALSE))) %>%
      ggplot(.) + geom_errorbar(aes(y = sim, xmin = CI_lwr, xmax = CI_upr, colour = CI_excludes_0)) +
      scale_colour_manual(values = c("red", "grey50"), breaks = c(TRUE, FALSE), guide = "none") +
      ggtitle("Robust mixture prior")  +
      xlab("Effect size") + theme_bw() +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
    
    p3.4 <- data %>% filter(dsg == "bayes") %>% ggplot(.) + geom_histogram(aes(x=best_pvalue), bins = 50, boundary = 0) + ggtitle("Robust mixture prior") +
      xlab("Smallest observed p-value") + ylab("Number of simulations") + theme_bw() + scale_x_continuous(limits = c(0,1))
 
    
    p4.4 <- rejection_rates %>% filter(dsg == "bayes") %>% ggplot(., aes(x = method, y = rejection_rate, fill = method)) +
      geom_col() + 
      geom_text(aes(label = round(rejection_rate, 3)), vjust = -0.5) +
      xlab("Method") + ylab("Rejection Rate") + ggtitle("Rejection Rate") + theme_bw() + coord_cartesian(ylim = c(0,1))
    
  }
  
  # first three plots for each correction method
  plots_mt <- results_standard_adj %>% group_split(multiple_test) %>%
    map(function(df) {
      method = unique(df$multiple_test)
      
      p1.m <- ggplot(df) + geom_histogram(aes(x = effect_est_glht), binwidth = 0.16) + ggtitle(paste(method)) +
        xlab("Observed effect estimate with smallest p-value") + ylab("Number of simulations") + 
        coord_cartesian(xlim = c(-4,4)) + #scale_x_continuous(limits = c(-4,4)) + 
        theme_bw()
      
      p2.m <- df %>% arrange(CI_glht_lwr) %>% rowid_to_column("sim") %>%
        mutate(CI_excludes_0 = case_when(CI_glht_lwr < 0 & CI_glht_upr > 0 ~ FALSE, TRUE ~ TRUE)) %>%
        mutate(CI_excludes_0 = factor(CI_excludes_0, levels = c(TRUE, FALSE))) %>%
        ggplot(.) + geom_errorbar(aes(y = sim, xmin = CI_glht_lwr, xmax = CI_glht_upr, colour = CI_excludes_0)) +
        scale_colour_manual(values = c("red", "grey50"), breaks = c(TRUE, FALSE), guide = "none") +
        ggtitle(paste(method))  +
        xlab("Effect size") + theme_bw() +
        theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
      
      p3.m <- ggplot(df) + geom_histogram(aes(x=best_pvalue_glht), bins = 50, boundary = 0) + ggtitle(paste(method)) +
        xlab("Smallest observed p-value") + ylab("Number of simulations") + theme_bw() + scale_x_continuous(limits = c(0,1))
      
      list(method = method, plots = list(p1.m,p2.m,p3.m))
    })
  
 
  # all plots (including rejection rate) for each correction method 
  plots_mt_full <- map(plots_mt, function(x) {
    method <- x$method
    plots <- x$plots
    
    # find rejection rate
    rejection_rate <- rejection_rates %>% 
      filter(dsg == "standard", method == !!method) %>%
      pull(rejection_rate)
    
    p4.m <- ggplot(tibble(method = method, rejection_rate = rejection_rate),
                   aes(x = method, y = rejection_rate, fill = method)) + 
      geom_col() +
      geom_text(aes(label = round(rejection_rate,3)), vjust = -0.5) +
      xlab("Method") + ylab ("Rejection Rate") + ggtitle("Rejection Rate") + theme_bw() + coord_cartesian(ylim = c(0,1))
    
    all_plots <- c(plots, list(p4.m))
    all_plots
      
  })
  
  # plots traffic light system
  if(!(is.null(data_tls_mse))){
  scenario_labeller <- labeller(
    dsg = c(
      repeat_and_pool = "Repeat and pool",
      repeat_and_replace = "Repeat and replace",
      standard = "Standard"
    ),
    n_per_group = function(x) {
      paste0("N per group:\n", x)
    }
  )
  
  
  p.tls <- data_tls_mse %>%
        ggplot(.) + geom_point(aes(x = no_groups, y = MSE)) +
    geom_rect(aes(xmin = no_groups-0.499, xmax = no_groups+0.501, ymin=0, ymax=Inf, fill = colors), alpha = 0.3) +
    facet_grid(dsg ~ n_per_group, labeller = scenario_labeller) +
    scale_fill_manual(values = c("#009e73", "#f0e442", "#d55e00"), breaks = c("green", "yellow", "red")) +
    scale_x_continuous(breaks = 2:30) + xlab("Number of groups") +
    theme_bw() + theme(panel.grid.minor = element_blank()) 
  

  
  traffic_light_plot <- arrangeGrob(p.tls, ncol = 1)
  }
  
  if(input$choose_bayes == "Yes"){
    if(input$traffic_light_system == "Yes") {
    small_plots_bayes <- arrangeGrob(
      grobs = c(list(p1.1, p2.1, p3.1, p4.1,
                     p1.2, p2.2, p3.2, p4.2,
                     p1.3, p2.3, p3.3, p4.3,
                     p1.4, p2.4, p3.4, p4.4),
                flatten(plots_mt_full)
      ), ncol = 4
    )
    grid.arrange(
      small_plots_bayes,
      traffic_light_plot, 
      ncol = 1, 
      heights = c(6,4)
    )} else {
      grid.arrange(grobs = c(
        list(p1.1, p2.1, p3.1, p4.1, 
             p1.2, p2.2, p3.2, p4.2, 
             p1.3, p2.3, p3.3, p4.3, 
             p1.4, p2.4, p3.4, p4.4), 
        flatten(plots_mt_full)), ncol = 4) 
    }
    
    

    } else if  (input$choose_bayes == "No"){
    if (input$traffic_light_system == "Yes") {
    small_plots <- arrangeGrob(
      grobs = c(list(p1.1, p2.1, p3.1, p4.1, 
                     p1.2, p2.2, p3.2, p4.2, 
                     p1.3, p2.3, p3.3, p4.3), 
                flatten(plots_mt_full)
      ), ncol = 4
    )
    grid.arrange(
      small_plots,
      traffic_light_plot,
      ncol = 1,
      heights = c(6,4)  
    )
    } else {
      grid.arrange(grobs = c(
        list(p1.1, p2.1, p3.1, p4.1, 
             p1.2, p2.2, p3.2, p4.2, 
             p1.3, p2.3, p3.3, p4.3), 
        flatten(plots_mt_full)), ncol = 4) 
  }
      

  }}, height = 1100)

  output$text <- renderText({
    all_combs <- gen_all_combs()
    paste0("all combinations:\n",paste0(all_combs, collapse=", "),"\nselected combinations:\n", paste0(input$custom_comparisons, collapse=", "))
    })
  
}

# Create Shiny app ----
shinyApp(ui, server)