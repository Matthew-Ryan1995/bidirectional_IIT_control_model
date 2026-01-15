## Matt Ryan
## 14/12/2023
## Run all scenario simulations under consideration
## Requires: pacman, tidyverse, here, glue, parallel, progress

# libraries ---------------------------------------------------------------

pacman::p_load(tidyverse)

source(("code/model_functions.R"))
source(here::here("Frieds_sensitivity/run_simulations_Fried.R"))

print(paste("start time ", Sys.time()))

# Load system variables ---------------------------------------------------

print(date())
## Detect cores from phoenix
slurm_ntasks <- as.numeric(Sys.getenv("SLURM_NTASKS")) # Obtain environment variable SLURM_NTASKS
if (!is.na(slurm_ntasks)) {
  mc.cores = slurm_ntasks # if slurm_ntasks is numerical, then assign it to cores
}else {
  mc.cores = 10 # Figure out how many cores there are
}

idx <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if(is.na(idx)){
  idx <- 1
}

print(str_c("Doing job ", idx))

# MIN_VALUES <- as.logical(Sys.getenv("MIN_VALUES"))
# MAX_VALUES <- as.logical(Sys.getenv("MAX_VALUES"))


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
# seed <- 230124
# num_runs <- 1
# start_time <- 0 
# end_time <- 200
# mc.cores <- parallel::detectCores() - 1 # CHANGE WHEN ON BOWEN


# Bowen
seed <- 230124
num_runs <- 500
start_time <- 0
end_time <- 2*365 + 180 # MM: added extra 6 months (+ 180)
# mc.cores <- 10# parallel::detectCores() - 1 # CHANGE WHEN ON BOWEN


# load parameters ---------------------------------------------------------

# Model parameters
if(GENERATE_MODEL_PARAMETERS){ # fixme: Due to the min/max flags and transformation flags, running this may not align with what you want to do.
  source(("code/create_model_parameters.R"))
}
param_list <- read_rds((param_load_name))
# Scenario parameters
if(GENERATE_SCENARIO_PARAMETERS){
  source(("code/create_scenario_parameters.R"))
}
scenario_parameters <- read_csv(("data/scenario_parameters.csv"), 
                                col_types = cols())

if(!dir.exists((glue::glue("outputs/Fried/intervention_simulations/{param_list$parameter_set}")))){
  dir.create((glue::glue("outputs/Fried/intervention_simulations/{param_list$parameter_set}")), recursive = TRUE)
}


# Run each simulation type ------------------------------------------------

# for(idx in 1:nrow(scenario_parameters)){

print(paste("Doing scenario ", idx))

# Pull out scenario parameters
OF <- scenario_parameters$OF[idx]
FCP <- scenario_parameters$FCP[idx]
female_immigration_rate <- scenario_parameters$female_immigration_rate[idx]
male_immigration_rate <- scenario_parameters$male_immigration_rate[idx]
control_proportion <- scenario_parameters$control_proportion[idx]
intervention_type <- scenario_parameters$intervention_type[idx]

# Set global parameter, turned off for some interventions 
# so we can define intervention function
# stop_intervention <- FALSE


# Pull out initial population for intervention function
## This is the carrying capacity of wild type adults in the population
initial_population <- sum(param_list$carryingCapacityByPopulation_vector[1]) # Hardcode 1 population

# Define intervention function for the scenario.  Easier to do here than elsewhere.
interventionFunction <- function(state_vector, startTime, endTime, interventionTime, stop_intervention,  # MM: need edits
                                 overfloodingDump = OF, femaleContamination = FCP,
                                 wPip_prop = control_proportion, intervention = intervention_type,
                                 init_prop = initial_population, ...){
  
  ARGS <- list(...)
  if(isTRUE(ARGS$timeout)){
    return(list(state_vector = state_vector, numReleased = 0, stop_intervention = stop_intervention))
  }
  
  # If maintain, start release again if proportion of wPip to carrying capacity drops below
  # a certain threshold (80% of wPip_prop)
  if(isTRUE(intervention=="maintain")){
    if(isTRUE(stop_intervention)){
      indAlb <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, 
                                        genotype = "wAlbAB", sex=c("f", "m"), ageClass = NA, mate = NA,
                                        mateStage = NA)
      indPip <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, 
                                        genotype = "wPip", sex=c("f", "m"), ageClass = NA, mate = NA,
                                        mateStage = NA)
      
      # Calculate proportion of adult wPip compared to carrying capacity in adult population
      propPip <- sum(state_vector[indPip])/init_prop
      
      propAlb <- sum(state_vector[indAlb])/init_prop
      
      # Start interventions again if the proportion of wPip in adults drops below
      # 80% of the establishment threshold proportion (wPip_prop) and
      # wAlb are not below suppression level.
      if(isTRUE((propPip < 0.8 * wPip_prop) & (propAlb > 0.1))){
        stop_intervention <- FALSE
      }
    }
  }
  
  # Don't run intervention
  if(isTRUE(stop_intervention)){
    return(list(state_vector = state_vector, numReleased = 0, stop_intervention = stop_intervention))
  }
  
  ## Do we need to stop our interventions
  ## Conditions:
  ### 1. only stop interventions if more than 100 days
  ### 2. If elimination, stop if the proportion of wAlb drops below 10% of carrying capacity
  ### 3. If complete_stop, stop when proportion of wPip adults to carrying capacity drops below wPip_prop (critical threshold)
  ### 4. If maintain, stop release if proportion of wPip adults to carrying capacity drops below wPip_prop (critical threshold)
  if(isTRUE(interventionTime > 100)){# MM: Edit
    # Elimination condition
    # This might be redundant code as I think this is hardcoded in runModel
    # if(isTRUE(intervention=="elimination")){ # MM: this is the naive intervention 
    indAlb <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, 
                                      genotype = "wAlbAB", sex=c("f", "m"), ageClass = NA, mate = NA,
                                      mateStage = NA)
    propAlb <- sum(state_vector[indAlb])/init_prop
    
    # Stop when wAlbAB drops below 10% of original population
    if(isTRUE(propAlb < 0.1)){ 
      stop_intervention <- TRUE
      return(list(state_vector = state_vector, numReleased = 0, stop_intervention = stop_intervention))
    }
    # }else{
    if(isFALSE(intervention=="elimination")){
      # Complete stop  and maintain condition condition
      
      indPip <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, 
                                        genotype = "wPip", sex=c("f", "m"), ageClass = NA, mate = NA,
                                        mateStage = NA)
      
      # Calculate proportion of adult wPip in adult population
      propPip <- sum(state_vector[indPip])/(init_prop)
      
      if(isTRUE(propPip > wPip_prop)){
        stop_intervention <- TRUE
        
        # # Start our ticking time bomb if we are doing the complete stop intervention
        # if(isTRUE(intervention=="complete_stop")){
        #   stop_sims <<- TRUE
        #   stop_counter <<- 0
        # }
        
        return(list(state_vector = state_vector, numReleased = 0, stop_intervention = stop_intervention))
      }
    }
  }
  
  
  # See how many adult wAlbAB, the "bad guys"
  badguyInds <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, 
                                        genotype = "wAlbAB", sex = "m", ageClass = NA, mate = NA,
                                        mateStage = NA)
  numBadGuys <- sum(state_vector[badguyInds])
  
  ## commented code counts how many young wPip there are, currently not in use
  goodguyInds <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA,
                                         genotype = "wPip", sex = "m", ageClass = 1, mate = NA,
                                         mateStage = NA)
  # numGoodGuys <- sum(state_vector[goodguyInds])
  
  # Find indices for where to put released females
  contaminatedFemaleInd <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, 
                                                   genotype = "wPip", sex = "f", ageClass = 1, 
                                                   mate = "Unmated", mateStage = NA) 
  
  
  # How many to release based on overflooding ratio
  numReleased <- overfloodingDump*numBadGuys #??- numGoodGuys
  # Randomise release number
  numMalesReleased <- rbinom(1, numReleased, (1 - femaleContamination))
  numFemalesReleased <- numReleased - numMalesReleased
  
  # Update state vector
  state_vector[goodguyInds] <- state_vector[goodguyInds] + numMalesReleased
  state_vector[contaminatedFemaleInd] <- state_vector[contaminatedFemaleInd] + numFemalesReleased
  
  
  
  return(list(state_vector = state_vector, numReleased = numReleased))
}



model_parameters <- param_list

# calculate immigration rate vector for scenario
immigration_and_emigration_rates <- create_immigration_emmigration_vectors(state_vector = model_parameters$state_vector, overallMaleImmigrationRate = male_immigration_rate, 
                                                                           overallFemaleImmigrationRate = female_immigration_rate, maleDeathRate = model_parameters$maleDeathRate, 
                                                                           numMaleClasses = model_parameters$numMaleClasses, phi = model_parameters$phi, theta = model_parameters$theta, M1 = model_parameters$M1,
                                                                           carryingCapacityByPopulation_vector = model_parameters$carryingCapacityByPopulation_vector)
immigrationRate_vector <- immigration_and_emigration_rates$immigrationRate_vector
emigrationRate_vector <- immigration_and_emigration_rates$emigrationRate_vector

model_parameters$immigrationRate_vector <- immigrationRate_vector
model_parameters$emigrationRate_vector <- emigrationRate_vector

res <- run_simulations(OF = OF, FCP = FCP, female_immigration_rate = female_immigration_rate,
                       male_immigration_rate = male_immigration_rate, control_proportion = control_proportion,
                       intervention_type = intervention_type, interventionFunction = interventionFunction,
                       model_parameters = model_parameters, seed = seed, num_runs = num_runs,
                       start_time = start_time, end_time = end_time, mc.cores = mc.cores)
# }
print(paste("End time ", Sys.time()))
