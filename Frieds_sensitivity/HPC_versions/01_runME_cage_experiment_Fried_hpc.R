## Matt Ryan
## 14/12/2023
## Run all scenario simulations under consideration
## Requires: pacman, tidyverse, here, glue, parallel, progress

# libraries ---------------------------------------------------------------

pacman::p_load(tidyverse)

source(here::here("code/model_functions.R"))
source(here::here("Frieds_sensitivity/run_simulations_cage_Fried.R"))


print(paste("start time ", Sys.time()))

# Load system variables ---------------------------------------------------

print(date())
## Detect cores from phoenix
slurm_ntasks <- as.numeric(Sys.getenv("SLURM_NTASKS")) # Obtain environment variable SLURM_NTASKS
if (!is.na(slurm_ntasks)) {
  mc.cores = slurm_ntasks # if slurm_ntasks is numerical, then assign it to cores
}else {
  mc.cores = 1 # Figure out how many cores there are
}

idx <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

print(str_c("Doing job ", idx))

# Flags -------------------------------------------------------------------
## Allows you to choose either expected value (both false) or min/max values of parameters
Fried_type <- "low"

## Allows you to choose whether to transform the birth rate or the death rate to allow for 
## biologically feasible conditions
TRANSFORM_DEATHRATE <- TRUE
TRANSFORM_BIRTHRATE <- FALSE

param_load_name <- glue::glue("data/param_list_Fried_{Fried_type}.Rds")


GENERATE_MODEL_PARAMETERS <- FALSE
GENERATE_SCENARIO_PARAMETERS <- FALSE

# variable definitions ----------------------------------------------------
# Local
# seed <- 141223
# num_runs <- 5
# start_time <- 0 
# end_time <- 300
# mc.cores <- parallel::detectCores() - 1 # CHANGE WHEN ON BOWEN

# Bowen
seed <- 141223
num_runs <- 500
start_time <- 0
end_time <- 500
# mc.cores <- 10# parallel::detectCores() - 1 # CHANGE WHEN ON BOWEN


# load parameters ---------------------------------------------------------

# Model parameters
if(GENERATE_MODEL_PARAMETERS){ # fixme: Due to the min/max flags and transformation flags, running this may not align with what you want to do.
  source(here::here("code/create_model_parameters.R"))
}
param_list <- read_rds(here::here(param_load_name))
# Scenario parameters
if(GENERATE_SCENARIO_PARAMETERS){
  source(here::here("code/create_scenario_parameters_cage.R"))
}
scenario_parameters <- read_csv(here::here("data/scenario_parameters_cage.csv"), 
                                col_types = cols()) %>% 
  filter(direction=="bi")



if(!dir.exists(here::here(glue::glue("outputs/Fried/cage_simulations/{param_list$parameter_set}")))){
  dir.create(here::here(glue::glue("outputs/Fried/cage_simulations/{param_list$parameter_set}")), recursive = TRUE)
}

# Run each simulation type ------------------------------------------------

# for(idx in 1:nrow(scenario_parameters)){
  
  print(paste("Doing scenario ", idx))
  
  # Pull out scenario parameters
  initialPopulation <- 420#scenario_parameters$initialPopulation[idx] #fixme: parameter choices in model_parameters depend on this carrying capacity
  propPip <- scenario_parameters$propPip[idx]
  femaleProportion <- scenario_parameters$femaleProportion[idx]
  direction <- scenario_parameters$direction[idx]
  # No immigration in cage scenario
  female_immigration_rate <- 0
  male_immigration_rate <- 0
  
  # Define intervention function for the scenario.  Easier to do here than elsewhere.
  interventionFunction <- function(state_vector, startTime, endTime, interventionTime, ...){
    # Don't run intervention
    return(list(state_vector = state_vector, numReleased = 0, stop_intervention=TRUE))
  }
  
  model_parameters <- param_list
  
  state_vector <- create_cage_state_vector(model_parameters$state_vector, propPip = propPip, 
                                           total_capacity = initialPopulation, female_split = femaleProportion)
  
  model_parameters$state_vector <- state_vector
  
  # calculate immigration rate vector for scenario
  immigration_and_emigration_rates <- create_immigration_emmigration_vectors(state_vector = model_parameters$state_vector, overallMaleImmigrationRate = male_immigration_rate, 
                                                                             overallFemaleImmigrationRate = female_immigration_rate, maleDeathRate = model_parameters$maleDeathRate, 
                                                                             numMaleClasses = model_parameters$numMaleClasses, phi = model_parameters$phi, 
                                                                             theta = model_parameters$theta, M1 = model_parameters$M1, carryingCapacityByPopulation_vector=model_parameters$carryingCapacityByPopulation_vector)
  immigrationRate_vector <- immigration_and_emigration_rates$immigrationRate_vector
  emigrationRate_vector <- immigration_and_emigration_rates$emigrationRate_vector
  
  model_parameters$immigrationRate_vector <- immigrationRate_vector
  model_parameters$emigrationRate_vector <- emigrationRate_vector
  
  if(isTRUE(direction=="uni")){
    #Creating the CI parameter vector.  Start by giving every type of individual 0 CI
    ## If a wAlbAB female is mated with a wPip, they cannot give birth (CI=1)
    ## If a wPip female is mated with a wAlbAB male, they can give birth (CI=0)
    ## This means that CI is uni-directional, it is in favour of the wPip
    ## Taking into account that we assume wolbachia strain is passed on from the females (wPip mother means wPip baby),
    populationNames <- param_list$populationNames
    stateNames <- names(model_parameters$state_vector)
    ciByName_vector <- rep(0, length(stateNames)) #don't change
    names(ciByName_vector) <- stateNames #don't change
    ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wPip", sex = "f", ageClass = NA, mate = "wAlbAB") #find the indices for state names for wPip females that were mated by wAlbAB males
    ciByName_vector[ind] <- 0 # give those females a cytoplasmic incompatibility of 1 (no offspring)
    ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wAlbAB", sex = "f", ageClass = NA, mate = "wPip") #find the indices for state names for wAlbAB females that were mated by wPip males
    ciByName_vector[ind] <- 1 # give those females a cytoplasmic incompatibility of 1 (no offspring)
    ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wPip", sex = "f", ageClass = NA, mate = "wAlbAB", mateStage = 15:19) #find the indices for state names for wPip females that were mated by wAlbAB males who were aged 15 to 19 days old
    ciByName_vector[ind] <- 0. # give those females a cytoplasmic incompatibility of 0.68 (diminished probability of offspring)
    ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wPip", sex = "f", ageClass = NA, mate = "wAlbAB", mateStage = 20) #find the indices for state names for wPip females that were mated by wAlbAB males who were aged 20 days old (males can't get any older, from 20 days onwards they are all acting the same way)
    ciByName_vector[ind] <- 0 # give those females a cytoplasmic incompatibility of 0 can produce offspring
    
    model_parameters$ciByName_vector <- ciByName_vector 
  }
  
  if(isTRUE(direction=="bi_no_decay")){
    #Creating the CI parameter vector.  Start by giving every type of individual 0 CI
    ## If a wAlbAB female is mated with a wPip, they cannot give birth (CI=1)
    ## If a wPip female is mated with a wAlbAB male, they can give birth (CI=0)
    ## This means that CI is uni-directional, it is in favour of the wPip
    ## Taking into account that we assume wolbachia strain is passed on from the females (wPip mother means wPip baby),
    populationNames <- param_list$populationNames
    stateNames <- names(model_parameters$state_vector)
    ciByName_vector <- rep(0, length(stateNames)) #don't change
    names(ciByName_vector) <- stateNames #don't change
    ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wPip", sex = "f", ageClass = NA, mate = "wAlbAB") #find the indices for state names for wPip females that were mated by wAlbAB males
    ciByName_vector[ind] <- 1 # give those females a cytoplasmic incompatibility of 1 (no offspring)
    ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wAlbAB", sex = "f", ageClass = NA, mate = "wPip") #find the indices for state names for wAlbAB females that were mated by wPip males
    ciByName_vector[ind] <- 1 # give those females a cytoplasmic incompatibility of 1 (no offspring)
    # ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wPip", sex = "f", ageClass = NA, mate = "wAlbAB", mateStage = 15:19) #find the indices for state names for wPip females that were mated by wAlbAB males who were aged 15 to 19 days old
    # ciByName_vector[ind] <- 0. # give those females a cytoplasmic incompatibility of 0.68 (diminished probability of offspring)
    # ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wPip", sex = "f", ageClass = NA, mate = "wAlbAB", mateStage = 20) #find the indices for state names for wPip females that were mated by wAlbAB males who were aged 20 days old (males can't get any older, from 20 days onwards they are all acting the same way)
    # ciByName_vector[ind] <- 0 # give those females a cytoplasmic incompatibility of 0 can produce offspring
    
    model_parameters$ciByName_vector <- ciByName_vector 
  }
  
  res <- run_simulations_cage(initialPopulation = initialPopulation, propPip = propPip, direction = direction,
                              femaleProportion = femaleProportion, interventionFunction = interventionFunction,
                              model_parameters = model_parameters, seed = seed, num_runs = num_runs,
                              start_time = start_time, end_time = end_time, mc.cores = mc.cores)
# }
print(paste("End time ", Sys.time()))

