## Matt Ryan
## 19/12/2023
## Code to take in a set of parameters and run multiple simulations under those conditions
## Require: pacman, tidyverse, glue, here, parallel


#' run_simulations_cage
#' Run simulations to emulate cage experiment.  Looking to identify proportion of wPip leading to 
#' establishment.
#'
#' @param initialPopulation - Integer.  Initial population being considered for file naming purposes.
#' @param propPip - Numeric.  Proportion of wPip in initial state.
#' @param femaleProportion - Numeric.  Proportion of female in initial state.
#' @param interventionFunction - Function.  The intervention function controlling mosquito release
#' @param model_parameters - List.  The model parameters from `create_model_parameters.R`.  Need to modify immigration rates outside function
#' @param seed - Integer.  Seed value, to be altered by OF and FCP for reproducibility
#' @param num_runs - Integer.  Number of simulation to run
#' @param start_time - Integer.  Day 0 of simulation.
#' @param end_time - Integer.  End time of simulation.
#' @param mc.cores - Integer.  Number of cores to run code on, defaults to 1
#'
#' @return saves outputs in ~/outputs
#' @export
#'
#' @examples
run_simulations_cage <- function(initialPopulation, propPip, femaleProportion, direction,
                                 interventionFunction, model_parameters, 
                                 seed=1234, num_runs = 100, start_time = 0, end_time = 600,
                                 mc.cores = 1){
  
  # expand list to accessible parameters, 
  # probably makes them global but that should be fine as long as we don't double parameterise
  list2env(model_parameters, envir = globalenv())
  
  # Name of scenario
  scenario <- glue::glue("initialPopulation_{initialPopulation}_propPip_{propPip}_femaleProportion_{femaleProportion}_direction_{direction}_start_{start_time}_end_{end_time}_seed_{seed}")
  if(!dir.exists(here::here(glue::glue("outputs/cage_simulations/{parameter_set}/scenario_{scenario}")))){
    dir.create(here::here(glue::glue("outputs/cage_simulations/{parameter_set}/scenario_{scenario}")))
  }
  
  # Set seed of reproducibility
  # RNGkind("L'Ecuyer-CMRG")
  # set.seed(round(seed * initialPopulation / propPip))# + num_runs))
  # 
  # Run the simulations
  res <- parallel::mclapply(
    1:num_runs,
    function(j){
      if(all(transitionRate_imm*I_bar/(femaleBirthRate*F_bar_mated) < 1)) # Make sure the parameters are biologically feasible
      {
        set.seed(round(seed * initialPopulation / propPip + j))
        file_name <- glue::glue("outputs/cage_simulations/{parameter_set}/scenario_{scenario}/outputs-run{j}.csv")
        if(file.exists(here::here(file_name))){ # Don't repeat if run has been done
          return(NULL)
        }
        print(glue::glue("Starting job {scenario} run {j}")) # Tracker
        modelOutputs <- runModel(state_vector, start_time, end_time, interventionFunction, 
                                 # Run interventions every 7 days, starting on day 1
                                 interventionTimes = seq(start_time+1, end_time, 7), 
                                 recordTimes = (start_time+1):end_time, numClassesBySexAndGenotype_matrix, 
                                 matingRateByName_vector, ciByName_vector, offspringGenotypeByName_vector, 
                                 friedsByName_vector, deathRatesByName_vector, birthRateByName_vector, 
                                 transitionRateByName_vector, emmigrationRate_matrix, 
                                 immigrationRate_vector, emigrationRate_vector, carryingCapacityByPopulation_vector, 
                                 ImaxByPopulation_vector, progress = FALSE, cage = TRUE)
        
        # Should I save modelOutputs as well?
        # Maybe for future work
        ts <- aggregateStates(names(state_vector), modelOutputs, populationNames, genotypeNames)
        ts$released_to_date <- modelOutputs[nrow(modelOutputs), ]
        # save
        write.csv(ts, here::here(file_name))
        print(glue::glue("Completed job {scenario} run {j}"))
      }
    },
    mc.cores = mc.cores # Use most available cores
  )
  
  return("Sims complete")
  
}