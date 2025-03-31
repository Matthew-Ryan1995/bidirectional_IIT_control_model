## Matt Ryan
## Run results for scenario one
## Required: here, parallel, glue


# load functions ----------------------------------------------------------

source(here::here("code/model_functions.R"))


# load parameters ---------------------------------------------------------

param_list <- readRDS(here::here("data/param_list.Rds"))
list2env(param_list, envir = globalenv())

OF <- 1
FCP <- 1e-2


# Intervention function ---------------------------------------------------

interventionFunction <- function(state_vector, startTime, endTime, interventionTime,
                                 overfloodingDump = OF, femaleContamination = FCP)
{
  # return(state_vector)
  badguyInds <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, genotype = "wAlbAB", sex = "m", ageClass = NA, mate = NA, mateStage = NA)
  goodguyInds <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, genotype = "wPip", sex = "m", ageClass = 1, mate = NA, mateStage = NA) #query: Why age class 1 here but not above?
  #In the following line, you need to decide what state the females are in, mated or unmated?
  contaminatedFemaleInd <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, genotype = "wPip", sex = "f", ageClass = 1, mate = "Unmated", mateStage = NA) #query: Why only wPip?
  
  numBadGuys <- sum(state_vector[badguyInds])
  numGoodGuys <- sum(state_vector[goodguyInds])
  # I have no idea where these numbers come from
  # Only release eggs if there are more bad guys than good guys?
  # Does this need a stop if negative?
  numReleased <- overfloodingDump*numBadGuys #??- numGoodGuys
  # Probabilisticly assign how many males released.
  numMalesReleased <- rbinom(1, numReleased, (1 - femaleContamination))
  numFemalesReleased <- numReleased - numMalesReleased
  state_vector[goodguyInds] <- state_vector[goodguyInds] + numMalesReleased
  state_vector[contaminatedFemaleInd] <- state_vector[contaminatedFemaleInd] + numFemalesReleased
  
  indPipFemale <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, genotype = "wPip", sex = "f", ageClass = NA, mate = NA, mateStage = NA)
  indAlbABFemale <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, genotype = "wAlbAB", sex = "f", ageClass = NA, mate = NA, mateStage = NA)
  
  propPip <- sum(state_vector[indPipFemale])/(sum(state_vector[indAlbABFemale]) + sum(state_vector[indPipFemale]))
  
  # Full on stop if more than 100 interventions or more than 50% pip females?
  if(isTRUE(interventionTime > 100 & propPip > 0.5))
  {
    state_vector <- state_vector*0  
    numReleased <- 0
  }
  
  
  
  return(list(state_vector = state_vector, numReleased = numReleased))
}



# Run simulations ---------------------------------------------------------

#Check that the parameter set meets parametric constraints that are necessary to ensure model behaves properly
scenario <- glue::glue("OF_{OF}_FCP_{FCP}")
if(!dir.exists(here::here(glue::glue("outputs/scenario_{scenario}")))){
  dir.create(here::here(glue::glue("outputs/scenario_{scenario}")))
}
num_runs <- 102
start_time <- 0
end_time <- 600

# Set seed of reproducibility
RNGkind("L'Ecuyer-CMRG")
set.seed(13112023 * OF / FCP)

res <- parallel::mclapply(
  1:num_runs,
  function(j){
    if(all(transitionRate_imm*I_bar/(femaleBirthRate*F_bar_mated) < 1))
    {
      if(file.exists(here::here(glue::glue("outputs/scenario_{scenario}"),paste("outputs-run",j,".csv")))){
        return(NULL)
      }
      print(paste0("Starting job ", j))
      
      # A hack fix, these results used no immigration or emigration
      immigrationRate_vector <- rep(0, length(state_vector))
      names(immigrationRate_vector) <- names(state_vector)
      emigrationRate_vector <- immigrationRate_vector
      
      modelOutputs <- runModel(state_vector = state_vector,
                               startTime = start_time, 
                               endTime = end_time,
                               interventionFunction = interventionFunction, 
                               interventionTimes = seq(start_time+1, end_time, 7), 
                               recordTimes = (start_time+1):end_time, 
                               numClassesBySexAndGenotype_matrix = numClassesBySexAndGenotype_matrix, 
                               matingRateByName_vector = matingRateByName_vector,
                               ciByName_vector = ciByName_vector,
                               offspringGenotypeByName_vector = offspringGenotypeByName_vector, 
                               friedsByName_vector = friedsByName_vector, 
                               deathRatesByName_vector = deathRatesByName_vector,
                               birthRateByName_vector = birthRateByName_vector, 
                               transitionRateByName_vector = transitionRateByName_vector,
                               emmigrationRate_matrix = emmigrationRate_matrix, 
                               immigrationRate_vector = immigrationRate_vector,
                               emigrationRate_vector = emigrationRate_vector,
                               carryingCapacityByPopulation_vector = carryingCapacityByPopulation_vector, 
                               ImaxByPopulation_vector = ImaxByPopulation_vector)
      
      
      ts <- aggregateStates(names(state_vector), modelOutputs, populationNames, genotypeNames)
      ts$released_to_date <- modelOutputs[nrow(modelOutputs), ]
      # save
      write.csv(ts, here::here(glue::glue("outputs/scenario_{scenario}"),paste("outputs-run",j,".csv")))
      print(paste0("Completed job ", j))
    }
  },
  mc.cores = parallel::detectCores() - 1
)


