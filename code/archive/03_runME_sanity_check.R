## Matt Ryan
## 23/01/2024
## Sanity check on the CTMC. Make sure steady states are maintained and don't do anything
## silly.  No wPip introduced, wAlbAB only
## Requires: pacman, tidyverse, here, glue, parallel, progress

# libraries ---------------------------------------------------------------

pacman::p_load(tidyverse)

source(here::here("code/model_functions.R"))
source(here::here("code/run_simulations_sanity_check.R"))


# Flags -------------------------------------------------------------------
## Allows you to choose either expected value (both false) or min/max values of parameters
MIN_VALUES <- FALSE
MAX_VALUES <- FALSE

## Allows you to choose whether to transform the birth rate or the death rate to allow for 
## biologically feasible conditions
TRANSFORM_DEATHRATE <- TRUE
TRANSFORM_BIRTHRATE <- FALSE

param_load_name <- glue::glue("data/param_list_MIN_VALUES_{MIN_VALUES}_MAX_VALUES_{MAX_VALUES}_TRANSFORM_DEATHRATE_{TRANSFORM_DEATHRATE}_TRANSFORM_BIRTHRATE_{TRANSFORM_BIRTHRATE}.Rds")


GENERATE_MODEL_PARAMETERS <- FALSE
GENERATE_SCENARIO_PARAMETERS <- FALSE

# variable definitions ----------------------------------------------------
# Local
seed <- 230124
num_runs <- 1
start_time <- 0 
end_time <- 100
mc.cores <- parallel::detectCores() - 1 # CHANGE WHEN ON BOWEN

# Bowen
# seed <- 230124
# num_runs <- 100
# start_time <- 0 
# end_time <- 500
# mc.cores <- 8# parallel::detectCores() - 1 # CHANGE WHEN ON BOWEN


# load parameters ---------------------------------------------------------

# Model parameters
if(GENERATE_MODEL_PARAMETERS){ # fixme: Due to the min/max flags and transformation flags, running this may not align with what you want to do.
  source(here::here("code/create_model_parameters.R"))
}
param_list <- read_rds(here::here(param_load_name))

if(!dir.exists(here::here("outputs/sanityCheck"))){
  dir.create(here::here("outputs/sanityCheck"))
}

if(!dir.exists(here::here(glue::glue("outputs/sanityCheck/{param_list$parameter_set}")))){
  dir.create(here::here(glue::glue("outputs/sanityCheck/{param_list$parameter_set}")))
}


SanityCheckDate <- today()

# Run each simulation type ------------------------------------------------


print(paste("Doing scenario on ", SanityCheckDate))

# No immigration in cage scenario
female_immigration_rate <- 0
male_immigration_rate <- 0

# Define intervention function for the scenario.  Easier to do here than elsewhere.
interventionFunction <- function(state_vector, startTime, endTime, interventionTime){
  # Don't run intervention
  return(list(state_vector = state_vector, numReleased = 0))
}

model_parameters <- param_list

# calculate immigration rate vector for scenario
immigration_and_emigration_rates <- create_immigration_emmigration_vectors(state_vector = model_parameters$state_vector, overallMaleImmigrationRate = male_immigration_rate, 
                                                                           overallFemaleImmigrationRate = female_immigration_rate, maleDeathRate = model_parameters$maleDeathRate, 
                                                                           numMaleClasses = model_parameters$numMaleClasses, phi = model_parameters$phi, 
                                                                           theta = model_parameters$theta, M1 = model_parameters$M1, carryingCapacityByPopulation_vector=model_parameters$carryingCapacityByPopulation_vector)
immigrationRate_vector <- immigration_and_emigration_rates$immigrationRate_vector
emigrationRate_vector <- immigration_and_emigration_rates$emigrationRate_vector

model_parameters$immigrationRate_vector <- immigrationRate_vector
model_parameters$emigrationRate_vector <- emigrationRate_vector

res <- run_simulations_sanity_check(SanityCheckDate = SanityCheckDate,
                                    interventionFunction = interventionFunction,
                                    model_parameters = model_parameters, 
                                    seed = seed, num_runs = num_runs,
                                    start_time = start_time, end_time = end_time, mc.cores = mc.cores)

