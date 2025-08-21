##############################################################################################
# this is the simulation code accompanying the manuscript
# "A novel approach to the design and sample size planning of animal experiments 
# based on treatment effect estimation"
# by Dario Zocholl, Henrike Solveen, and Matthias Schmid

# authors of this simulation code: Dario Zocholl & Henrike Solveen
# corresponding author: Dario Zocholl, dario.zocholl@imbie.uni-bonn.de
# last modifications: Aug 19, 2025
##############################################################################################


# load libraries
library(tidyverse)
library(future)
library(future.apply)
library(furrr)
library(rstan)
library(emmeans)
library(multcomp)
library(ggpubr)
library(ggh4x)


# load functions
source("sample_data.R")
source("find_best_comparison.R")
source("estimate_effect.R")
source("get_data.R")
source("simulation_wrapper.R")


# set up simulation
source("simulation_setup.R")


# run simulation
source("run_simulation.R")


# load data (optional instead of running simulation)
data <- readRDS("C:/Users/zocholl/Documents/GitHub/SimBA/simulated_data_2025-05-27.Rds")

# prepare data for visualization
source("data_manipulation.R")

# Figures
source("Figure_2.R")
source("Figure_3.R")
source("Figure_7.R")
source("Figure_8.R")
source("Figure_9.R")
source("Figure_10.R")



# set up simulation for supplementary material
source("simulation_setup_supplementary_material.R")

# run simulation for supplementary material
source("run_simulation.R")
