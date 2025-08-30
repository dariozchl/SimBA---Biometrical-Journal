This is the simulation code accompanying the manuscript
"A novel approach to the design and sample size planning of animal experiments based on treatment effect estimation"
by Dario Zocholl, Henrike Solveen, and Matthias Schmid

Authors of this simulation code: Dario Zocholl & Henrike Solveen
Corresponding author: Dario Zocholl, dario.zocholl@imbie.uni-bonn.de
Last modifications: Aug 29, 2025


All results and figures in the manuscript can be reproduced by running "Master File.R".
Since the simulations require significant computation time (approx. 14 hours distributed to 15 cores),
it is also possible to just load the data or run it with fewer simulations by manually editing the 
variable 'nsim' in "simulation_setup.R" and "simulation_setup_supplementary_material.R.

The following files are included:
File name                         Purpose of the file
"Master File.R"                   Master File
"sample_data.R"                   Function to sample data
"find_best_comparison.R"          Function to find the best comparison
"estimate_effect.R"               Function to estimate the effect
"get_data.R"                      Function to generate some data
"simulation_wrapper.R"            Function to combine all previously defined functions
"simulation_setup.R"              Script to set up the parameters for the simulation study
"run_simulation.R"                Function to perform the simulations
"simulated_data_2025-05-27.Rds"   Simulation results 
"data_manipulation.R"             Script to prepare the data for visualization
"Figure_2.R"                      Script to create Figure 2
"Figure_3.R"                      Script to create Figure 3
"Figure_7.R"                      Script to create Figure 7
"Figure_8.R"                      Script to create Figure 8
"Figure_9.R"                      Script to create Figure 9
"Figure_10.R"                     Script to create Figure 10
"simulation_setup_supplementary_material.R" Set up the parameters for for supplementary material simulation. 
"run_simulation_supplementary_material.R." Script to run simulations for supplementary material.
"figures_supplementary_material.R" Make figures for supplementary material. 
"app.R"                           Script to run the shiny app version which is available online.
                                  Note that this is a computationally limited version. 
"app_bayesian_model.R"            Script to run the shiny app with full functionality.
"mixture_model_ttest.stan"        Stan-model for robust mixture priors


The subfolder "figures" contains figures produced by running the corresponding script 
in the Master file.

The subfolder "data" contains data produced by running the scripts 
run_simulation.R and run_simulation_supplementary_material.R. 



The version information is provided in the following: 

> sessionInfo()
R version 4.3.3 (2024-02-29 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 26100)

Matrix products: default


Random number generation:
 RNG:     L'Ecuyer-CMRG 
 Normal:  Inversion 
 Sample:  Rejection 
 
locale:
[1] LC_COLLATE=German_Germany.utf8  LC_CTYPE=German_Germany.utf8    LC_MONETARY=German_Germany.utf8
[4] LC_NUMERIC=C                    LC_TIME=German_Germany.utf8    

time zone: Europe/Berlin
tzcode source: internal

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] future.apply_1.11.2 furrr_0.3.1         future_1.34.0       doParallel_1.0.17   iterators_1.0.14   
 [6] snow_0.4-4          ggh4x_0.2.8         rlist_0.4.6.2       foreach_1.5.2       rstan_2.32.5       
[11] StanHeaders_2.32.5  emmeans_1.10.3      multcomp_1.4-25     TH.data_1.1-2       MASS_7.3-60.0.1    
[16] survival_3.5-8      mvtnorm_1.2-4       gridExtra_2.3       ggpubr_0.6.0        lubridate_1.9.3    
[21] forcats_1.0.0       stringr_1.5.1       dplyr_1.1.4         purrr_1.1.0         readr_2.1.5        
[26] tidyr_1.3.1         tibble_3.2.1        ggplot2_3.5.0       tidyverse_2.0.0     shinyjs_2.1.0      
[31] shiny_1.8.1.1      

loaded via a namespace (and not attached):
 [1] inline_0.3.19      sandwich_3.1-0     rlang_1.1.3        magrittr_2.0.3     matrixStats_1.2.0 
 [6] compiler_4.3.3     loo_2.7.0          systemfonts_1.0.5  vctrs_0.6.5        pkgconfig_2.0.3   
[11] fastmap_1.1.1      backports_1.4.1    ellipsis_0.3.2     labeling_0.4.3     utf8_1.2.4        
[16] promises_1.2.1     tzdb_0.4.0         ragg_1.2.7         cachem_1.0.8       jsonlite_1.8.8    
[21] later_1.3.2        rlecuyer_0.3-8     broom_1.0.5        R6_2.5.1           bslib_0.6.1       
[26] stringi_1.8.3      parallelly_1.38.0  car_3.1-2          jquerylib_0.1.4    estimability_1.5.1
[31] Rcpp_1.0.12        zoo_1.8-12         httpuv_1.6.14      Matrix_1.6-5       splines_4.3.3     
[36] timechange_0.3.0   tidyselect_1.2.0   rstudioapi_0.15.0  abind_1.4-5        codetools_0.2-19  
[41] curl_5.2.0         listenv_0.9.1      pkgbuild_1.4.3     lattice_0.22-5     withr_3.0.0       
[46] coda_0.19-4.1      RcppParallel_5.1.7 pillar_1.9.0       carData_3.0-5      rsconnect_1.4.1   
[51] stats4_4.3.3       generics_0.1.3     hms_1.1.3          munsell_0.5.0      scales_1.3.0      
[56] globals_0.16.3     xtable_1.8-4       glue_1.7.0         tools_4.3.3        data.table_1.15.0 
[61] ggsignif_0.6.4     cowplot_1.1.3      grid_4.3.3         QuickJSR_1.1.3     colorspace_2.1-0  
[66] cli_3.6.2          textshaping_0.3.7  fansi_1.0.6        V8_4.4.2           gtable_0.3.4      
[71] rstatix_0.7.2      sass_0.4.8         digest_0.6.34      farver_2.1.1       memoise_2.0.1     
[76] htmltools_0.5.7    lifecycle_1.0.4    mime_0.12         