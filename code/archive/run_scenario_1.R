## Matt Ryan
## Run results for scenario one
## Required: here, parallel, glue


# load functions ----------------------------------------------------------

source(here::here("code/model_functions.R"))


# load parameters ---------------------------------------------------------

param_list <- readRDS(here::here("data/param_list.Rds"))
list2env(param_list, envir = globalenv())


# Run simulations ---------------------------------------------------------

#Check that the parameter set meets parametric constraints that are necessary to ensure model behaves properly
scenario <-1
if(!dir.exists(here::here(glue::glue("outputs/scenario_{scenario}")))){
  dir.create(here::here(glue::glue("outputs/scenario_{scenario}")))
}
num_runs <- 2
start_time <- 0
end_time <- 10

RNGkind("L'Ecuyer-CMRG")
set.seed(13112023)

res <- parallel::mclapply(
  1:num_runs,
  function(j){
    if(all(transitionRate_imm*I_bar/(femaleBirthRate*F_bar_mated) < 1))
    {
      print(paste0("Starting job ", j))
      modelOutputs <- runModel(state_vector, start_time, end_time, interventionFunction, 
                               interventionTimes = seq(start_time+1, end_time, 7), 
                               recordTimes = (start_time+1):end_time, numClassesBySexAndGenotype_matrix, 
                               matingRateByName_vector, ciByName_vector, offspringGenotypeByName_vector, 
                               friedsByName_vector, deathRatesByName_vector, birthRateByName_vector, 
                               transitionRateByName_vector, emmigrationRate_matrix, 
                               immigrationRate_vector, carryingCapacityByPopulation_vector, 
                               ImaxByPopulation_vector)
      
      
      ts <- aggregateStates(names(state_vector), modelOutputs, populationNames, genotypeNames)
      ts$released_to_date <- modelOutputs[nrow(modelOutputs), ]
      # save
      write.csv(ts, here::here(glue::glue("outputs/scenario_{scenario}"),paste("Test-outputs-run",j,".csv")))
      print(paste0("Completed job ", j))
      
      
    }
  },
  mc.cores = parallel::detectCores() - 1
)


