## Requires: glue, here, readxl, progress
## Modified by Matt Ryan
##
library(progress)
#' getIndicesForStateQuery
#' Finds which indices in the state names match the given queries of interest, 
#' i.e. Wolbachia strain
#'
#' @param stateNames 
#' @param populationName 
#' @param genotype 
#' @param sex 
#' @param ageClass 
#' @param mate 
#' @param mateStage 
#'
#' @return
#' @export
#'
#' @examples
getIndicesForStateQuery <- function(stateNames, populationName = NA, genotype = NA, 
                                    sex = NA, ageClass = NA, mate = NA, mateStage = NA)
{
  match_indices <- c()
  for(i in 1:length(stateNames))
  {
    bits <- unlist(strsplit(stateNames[i], "_"))
    
    this_populationName <- bits[1]
    this_genotype <- bits[2]
    this_sex <- bits[3]
    this_ageClass <- bits[4]
    this_mate <- bits[5]
    this_mateStage <- bits[6]
    
    test1 <- FALSE
    if(all(is.na(populationName)) | this_populationName %in% populationName)
    {
      test1 <- TRUE
    }
    
    test2 <- FALSE
    if(all(is.na(genotype)) | this_genotype %in% genotype)
    {
      test2 <- TRUE
    }
    
    test3 <- FALSE
    if(all(is.na(sex)) | this_sex %in% sex)
    {
      test3 <- TRUE
    }
    
    test4 <- FALSE
    if(all(is.na(ageClass)) | this_ageClass %in% ageClass)
    {
      test4 <- TRUE
    }
    
    test5 <- FALSE
    if(all(is.na(mate)) | this_mate %in% mate)
    {
      test5 <- TRUE
    }
    
    test6 <- FALSE
    if(all(is.na(mateStage)) | this_mateStage %in% mateStage)
    {
      test6 <- TRUE
    }
    
    if(test1 & test2 & test3 & test4 & test5 & test6)
    {
      match_indices <- c(match_indices, i)
    }
  }
  return(match_indices)
}

aggregateStates <- function(stateNames, modelOutputs, populationNames, genotypeNames)
{
  timeSeries <- list()
  # For each population and genotpye
  for(p in 1:length(populationNames))
  {
    for(g in 1:length(genotypeNames))
    {
      # find males, females, and immatures
      theseInds_m <- getIndicesForStateQuery(stateNames, populationName = populationNames[p], 
                                             genotype = genotypeNames[g], sex = "m", 
                                             ageClass = NA, mate = NA)
      theseInds_f <- getIndicesForStateQuery(stateNames, populationName = populationNames[p], 
                                             genotype = genotypeNames[g], sex = "f", 
                                             ageClass = NA, mate = NA)
      theseInds_imm <- getIndicesForStateQuery(stateNames, populationName = populationNames[p], 
                                               genotype = genotypeNames[g], sex = "imm", 
                                               ageClass = NA, mate = NA)
      
      # Add along all model outputs for particular genotype population, and age
      if(length(theseInds_m) == 1){
        timeSeries[[paste0(populationNames[p], "_", genotypeNames[g], "_m")]] <- modelOutputs[theseInds_m, ]
      }else{
        timeSeries[[paste0(populationNames[p], "_", genotypeNames[g], "_m")]] <- apply(modelOutputs[theseInds_m, ], 2, sum) 
      }
      timeSeries[[paste0(populationNames[p], "_", genotypeNames[g], "_f")]] <- apply(modelOutputs[theseInds_f, ], 2, sum)
      timeSeries[[paste0(populationNames[p], "_", genotypeNames[g], "_imm")]] <- apply(modelOutputs[theseInds_imm, ], 2, sum)
    }
  }
  return(timeSeries)
}


getStateNames <- function(n_m, n_f, n_i, genotypeNames, populationNames, numClassesBySexAndGenotype_matrix)
{
  #Each population contains:
  #males of different genotypes and with n_m classes for each
  #females of different genotypes and with n_f classes for each
  #immatures of different genotypes and with n_i classes for each (equal numbers for male and female immatures)
  
  stateNames <- c()
  for(p in 1:length(populationNames))
  {
    for(g in 1:length(genotypeNames))
    {
      stateNames <- c(stateNames, paste0(populationNames[p], "_", genotypeNames[g], "_m_", 1:n_m[g]), paste0(populationNames[p], "_", genotypeNames[g], "_f_", 1:n_f[g], "_Unmated_None"), paste0(populationNames[p], "_", genotypeNames[g], "_imm_", 1:n_i[g]))
      for(m in 1:length(genotypeNames))
      {
        for(mc in 1:numClassesBySexAndGenotype_matrix[genotypeNames[m], "m"])
        {
          stateNames <- c(stateNames, paste0(populationNames[p], "_", genotypeNames[g], "_f_", 1:n_f[g], "_", genotypeNames[m], "_", mc))
        }
      }
    }
  }
  return(stateNames)
}


getTransitionRates <- function(state_vector, numClassesBySexAndGenotype_matrix, matingRateByName_vector, 
                               ciByName_vector, offspringGenotypeByName_vector, friedsByName_vector, 
                               deathRatesByName_vector, birthRateByName_vector, transitionRateByName_vector, 
                               emmigrationRate_matrix, immigrationRate_vector, emigrationRate_vector, carryingCapacityByPopulation_vector, 
                               ImaxByPopulation_vector, propFemale=1/2, propMale=1-propFemale)
{
  #Name Convention For States: populationName_genotype_sex_stage_mate_mateStage (where _mate only exists for the female sex and should be set to _Unmated for unmated females)
  #Populations: these are individual subpopulations in a metapopulation
  #Genotype: this can be used to describe wildtype and wolbachia infected types of individuals
  #Sex: this should be one of m, f or imm corresponding to male, female of immature individuals.
  #Stage: this is a numeric value that corresponds to one or more "stages" that individuals can progress through.  For immatures, this corresponds to development time
  #       before emerging as an adult, and for male and female types, this corresponds to adult life stages that can be given different parameters (CI, frieds index, mating rate, birth rate, death rate etc)
  #Mate: the type by which a female has been mated (use "Unmated" for virgin females).  This suffix is not present for male and immature states.
  
  #Inputs:
  #state_vector: a named vector (named using the state naming convention) containing the number of individuals in each class
  #numClassesBySexAndGenotype_matrix: a matrix with named rows (genotypes) and columns (sex - one of "m", "f" or "imm")
  #matingRateByName_vector: a named vector containing the mating rate for all female states
  #ciByName_vector: a named vector containing the cytoplasmic incompatibility for all female states (1 indicates full incompatibility, 0 indicates full compatibility, numbers between can be used for partial CI)
  #offspringGenotypeByName_vector: a named vector that gives the genotype of immature that can be birthed by a mated female of this state name.
  #friedsByName_vector: a named vector that gives the mating competitiveness coefficient (frieds index) for males of different state names
  #deathRatesByName_vector: a named vector that gives the per capita death rates of all male and female state names
  #birthRateByName_vector: a named vector that gives the per capita birth rates of all male and female state names
  #transitionRateByName_vector: a named vector that gives the per capita rate of transitions through stages for all male, female and immature state names
  #emmigrationRate_matrix: a matrix (named rows and columns) containing the per capita rate at which males and females emmigrate from one population to another.
  #                        row names are the "from" populations and column names are the "to" populations.
  #
  #immigrationRate_vector: a named vector containing the actual rates (not per capita rates) at which a new individual should appear in this state from "outside" (i.e. not one of the other populations int he model)
  #carryingCapacityByPopulation_vector: a vector containing the carrying capacity of a wildtype population in each of the populations
  
  #Return two items: one is a vector of transition rates; and two is a matrix whose columns correspond to the state-change
  
  ## Symbol definition
  # g - Number of genus
  # n_{m, i} - Number of males classes for genus i
  # n_{f, i} - Number of female classes for genus i
  # n_{I, i} - Number of immature classes for genus i
  # w - wildtype genus
  # () - transition rate
  # [] - further explanation
  #
  ## Male transitions
  # Each male can die ($\sum_{i = 1}^g n_{m, i}$)
  # Each male except the oldest can transition ($\sum_{i = 1}^g (n_{m, i} - 1)$)
  # Each male can experience emigration ($\sum_{i = 1}^g n_{m, i}$)
  # Wild type males can experience immigration ($n_{m, w}$)
  
  ## Female transitions
  # Each female can die ($\sum_{j = 1}^g n_{f, j} (1 + \sum_{i = 1}^g n_{m, i})$) [1 : unmated, n_{m, i} : for each mating option]
  # Each female can transition except the oldest ($\sum_{j = 1}^g (n_{f, j} - 1) (1 + \sum_{i = 1}^g n_{m, i})$)
  # Each female can experience emigration ($\sum_{j = 1}^g n_{f, j} (1 + \sum_{i = 1}^g n_{m, i})$)
  # Wild type female can experience immigration ($n_{f, w} (1 + n_{m, w})$)
  # Unmated females can mate ($\sum_{i, j = 1}^g n_{f, j} n_{m, i}$)
  
  ## Immature transitions
  # Immatures can be birthed ($\sum_{i, j = 1}^g n_{f, j} n_{m, i}$) [Assuming worst case of no CI]
  # Immatures can transition, except for oldest class ($\sum_{k = 1}^g (n_{I, k} - 1)$)
  # Immatures can transition to adulthood ($2g$) [can become male or female for each genus]
  
  
  # The following is a bunch of counting to see how many possible transitions can occur.
  stateNames_vector <- names(state_vector)
  
  populationNames <- names(carryingCapacityByPopulation_vector)
  numPopulations <- length(carryingCapacityByPopulation_vector) 
  
  genotypeNames <- rownames(numClassesBySexAndGenotype_matrix)
  
  totalMaleDeathTransitions <- sum(numClassesBySexAndGenotype_matrix[, "m"])
  totalMaleClassTransitions <- sum(numClassesBySexAndGenotype_matrix[, "m"] - 1)
  totalPossibleMaleEmmigrations <- sum(numClassesBySexAndGenotype_matrix[, "m"])
  totalMaleImmigrationTransitions <- numClassesBySexAndGenotype_matrix["wAlbAB", "m"] # wild type hardcoded
  
  totalFemaleDeathTransitions <- sum(numClassesBySexAndGenotype_matrix[, "f"] * (1 + sum(numClassesBySexAndGenotype_matrix[, "m"])))
  totalFemaleClassTransitions <- sum((numClassesBySexAndGenotype_matrix[, "f"] - 1) * (1 + sum(numClassesBySexAndGenotype_matrix[, "m"])))
  totalPossibleFemaleEmmigrations <- sum(numClassesBySexAndGenotype_matrix[, "f"] * (1 + sum(numClassesBySexAndGenotype_matrix[, "m"])))
  totalFemaleImmigrationTransitions <- numClassesBySexAndGenotype_matrix["wAlbAB", "f"] * (1 + numClassesBySexAndGenotype_matrix["wAlbAB", "m"])  # wild type hardcoded
  totalFemaleMatingTransitions <- sum(numClassesBySexAndGenotype_matrix[, "f"] * sum(numClassesBySexAndGenotype_matrix[, "m"]))
  
  totalImmatureBirthTransitions <- sum(numClassesBySexAndGenotype_matrix[, "f"] * (sum(numClassesBySexAndGenotype_matrix[, "m"])))
  totalImmatureClassTransitions <- sum(numClassesBySexAndGenotype_matrix[, "imm"] - 1)
  totalAdultBirthTransitions <- 2 * length(genotypeNames)
  n_transitions <- numPopulations*(totalMaleImmigrationTransitions + totalFemaleImmigrationTransitions + totalPossibleMaleEmmigrations + totalPossibleFemaleEmmigrations + totalMaleDeathTransitions + totalMaleClassTransitions + totalFemaleDeathTransitions + totalFemaleClassTransitions + totalFemaleMatingTransitions +  totalImmatureBirthTransitions + totalImmatureClassTransitions + totalAdultBirthTransitions)
  
  transitionRates_vector <- rep(0, n_transitions)
  transitions <- matrix(0, length(state_vector), n_transitions)
  
  maleIndices <- which(grepl( "_m_", stateNames_vector, fixed = TRUE) & state_vector > 0)
  femaleIndices <- which(grepl( "_f_", stateNames_vector, fixed = TRUE) & state_vector > 0)
  immatureIndices <- which(grepl( "_imm_", stateNames_vector, fixed = TRUE) & state_vector > 0)
  
  # Hardcoding that emigration = immigration
  # Emigration coded wrong
  ## CURRENTLY OVER-INFLATING MATED FEMALE MOVEMENT
  ## Emigration is currently wrong.  Does not track emigration of mates alb x pip
  # Overthinking this
  # emigrationRate_vector <- immigrationRate_vector
  # wPip_ind <- getIndicesForStateQuery(stateNames, genotype = "wPip", sex = c("m", "f"))
  # wAlbAB_ind <- getIndicesForStateQuery(stateNames, genotype = "wAlbAB", sex = c("m", "f"))
  # f_wAlbAB_wAlbAB_ind <- getIndicesForStateQuery(stateNames, genotype = "wAlbAB", sex = "f", 
  #                                                mate = "wAlbAB")
  # f_wAlbAB_wPip_ind <- getIndicesForStateQuery(stateNames, genotype = "wAlbAB", sex = "f", 
  #                                              mate = "wPip")
  # emigrationRate_vector[f_wAlbAB_wPip_ind] <- emigrationRate_vector[f_wAlbAB_wAlbAB_ind]
  # emigrationRate_vector[wPip_ind] <- emigrationRate_vector[wAlbAB_ind]
  
  #Start a counter for which transition we are up to in the matrix
  counter <- 1
  
  #All transitions for male states with at least one individual
  #We are also going to store each of the relevant male names by population
  relevantMalesByPopulation <- list()
  if(length(maleIndices) > 0)
  {
    for(i in 1:length(maleIndices))
    {
      #all male death transitions
      transitionRates_vector[counter] <- state_vector[maleIndices[i]]*deathRatesByName_vector[stateNames_vector[maleIndices[i]]]
      transitions[maleIndices[i], counter] <- -1
      names(transitionRates_vector)[counter] <- paste0("death_", stateNames_vector[maleIndices[i]])
      counter <- counter + 1
      
      #all male age transitions
      ind <- maleIndices[i]
      bits <- unlist(strsplit(stateNames_vector[ind], "_"))
      populationName <- bits[1]
      genotype <- bits[2]
      ageClass <- as.numeric(bits[4])
      relevantMalesByPopulation[[populationName]] <- c(relevantMalesByPopulation[[populationName]], stateNames_vector[maleIndices[i]])
      #print(numClassesBySexAndGenotype_matrix[genotype, "m"])
      if(ageClass < numClassesBySexAndGenotype_matrix[genotype, "m"])
      {
        nextInd <- which(stateNames_vector == paste0(populationName, "_", genotype, "_m_", ageClass + 1))
        transitionRates_vector[counter] <- state_vector[ind]*transitionRateByName_vector[stateNames_vector[ind]]
        transitions[ind, counter] <- -1
        transitions[nextInd, counter] <- 1
        counter <- counter + 1
      }
      
      #all male emmigration
      connections <- colnames(emmigrationRate_matrix)[which(emmigrationRate_matrix[populationName, ] > 0)]
      if(length(connections) > 0)
      {
        for(c in 1:length(connections))
        {
          connectedInd <- which(stateNames_vector == paste0(connections[c], "_", genotype, "_m_", ageClass))
          transitionRates_vector[counter] <- state_vector[ind]*emmigrationRate_matrix[populationName, connections[c]] 
          transitions[ind, counter] <- -1
          transitions[connectedInd, counter] <- 1
          counter <- counter + 1
        }
      }
    }
  }
  
  
  
  #All transitions for immature states with at least one individual
  #Need to count the number of immatures in each population
  ItotalByPopulations <- rep(0, length(populationNames))
  names(ItotalByPopulations) <- populationNames
  if(length(immatureIndices) > 0)
  {
    for(i in 1:length(immatureIndices))
    {
      ind <- immatureIndices[i]
      bits <- unlist(strsplit(stateNames_vector[ind], "_"))
      populationName <- bits[1]
      genotype <- bits[2]
      ageClass <- as.numeric(bits[4])
      
      ItotalByPopulations[[populationName]] <- ItotalByPopulations[[populationName]] + state_vector[ind]
      
      if(ageClass < numClassesBySexAndGenotype_matrix[genotype, "imm"])
      {
        nextInd <- which(stateNames_vector == paste0(populationName, "_", genotype, "_imm_", ageClass + 1))
        transitionRates_vector[counter] <- state_vector[ind]*transitionRateByName_vector[stateNames_vector[ind]]
        transitions[ind, counter] <- -1
        transitions[nextInd, counter] <- 1
        counter <- counter + 1
      }
      if(ageClass == numClassesBySexAndGenotype_matrix[genotype, "imm"])
      { # Note the 50% chance of male vs female appearing through 0.5*state_vector
        nextInd_m <- which(stateNames_vector == paste0(populationName, "_", genotype, "_m_", 1))
        transitionRates_vector[counter] <- propMale*state_vector[ind]*transitionRateByName_vector[stateNames_vector[ind]]
        transitions[ind, counter] <- -1
        transitions[nextInd_m, counter] <- 1
        counter <- counter + 1
        
        nextInd_f <- which(stateNames_vector == paste0(populationName, "_", genotype, "_f_", 1, "_Unmated_None"))
        transitionRates_vector[counter] <- propFemale*state_vector[ind]*transitionRateByName_vector[stateNames_vector[ind]]
        transitions[ind, counter] <- -1
        transitions[nextInd_f, counter] <- 1
        counter <- counter + 1
      }
    }
  }
  
  
  
  #All transitions for female states with at least one individual
  
  if(length(femaleIndices) > 0)
  {
    for(i in 1:length(femaleIndices))
    {
      #all female death transitions
      transitionRates_vector[counter] <- state_vector[femaleIndices[i]]*deathRatesByName_vector[stateNames_vector[femaleIndices[i]]]
      transitions[femaleIndices[i], counter] <- -1
      names(transitionRates_vector)[counter] <- paste0("death_", stateNames_vector[femaleIndices[i]])
      counter <- counter + 1
      
      #all female age transitions
      ind <- femaleIndices[i]
      bits <- unlist(strsplit(stateNames_vector[ind], "_"))
      populationName <- bits[1]
      genotype <- bits[2]
      ageClass <- as.numeric(bits[4])
      mate <- bits[5]
      mateClass <- bits[6]
      
      
      if(ageClass < numClassesBySexAndGenotype_matrix[genotype, "f"])
      {
        nextInd <- which(stateNames_vector == paste0(populationName, "_", genotype, "_f_", ageClass + 1, "_", mate, "_", mateClass))
        transitionRates_vector[counter] <- state_vector[ind]*transitionRateByName_vector[stateNames_vector[ind]]
        transitions[ind, counter] <- -1
        transitions[nextInd, counter] <- 1
        counter <- counter + 1
      }
      
      #all female emmigration
      connections <- colnames(emmigrationRate_matrix)[which(emmigrationRate_matrix[populationName, ] > 0)]
      if(length(connections) > 0)
      {
        for(c in 1:length(connections))
        {
          connectedInd <- which(stateNames_vector == paste0(connections[c], "_", genotype, "_f_", ageClass, "_", mate, "_", mateClass))
          transitionRates_vector[counter] <- state_vector[ind]*emmigrationRate_matrix[populationName, connections[c]]
          transitions[ind, counter] <- -1
          transitions[connectedInd, counter] <- 1
          counter <- counter + 1
        }
      }
      
      #female matings
      if(mate == "Unmated")
      {
        relevantMaleStates <- relevantMalesByPopulation[[populationName]]
        
        if(!is.null(relevantMaleStates))
        {
          maleGenotypeAndStageMatingPressure <- matrix(0, length(genotypeNames), max(numClassesBySexAndGenotype_matrix[, "m"]))
          rownames(maleGenotypeAndStageMatingPressure) <- genotypeNames
          for(j in 1:length(relevantMaleStates))
          {
            m_bits <- unlist(strsplit(relevantMaleStates[j], "_"))
            m_genotype <- m_bits[2]
            m_ageClass <- as.numeric(m_bits[4])
            maleGenotypeAndStageMatingPressure[m_genotype, m_ageClass] <- maleGenotypeAndStageMatingPressure[m_genotype, m_ageClass] + state_vector[relevantMaleStates[j]]*friedsByName_vector[relevantMaleStates[j]]
          }
          
          for(j in 1:length(genotypeNames))
          {
            for(k in 1:numClassesBySexAndGenotype_matrix[genotypeNames[j], "m"])
            {
              m_genotype <- genotypeNames[j]
              m_ageClass <- k
              # sets C_{geno} = 1?
              transitionRates_vector[counter] <- state_vector[ind]*matingRateByName_vector[stateNames_vector[ind]]*maleGenotypeAndStageMatingPressure[m_genotype, m_ageClass]/sum(maleGenotypeAndStageMatingPressure)
              matedInd <- which(stateNames_vector == paste0(populationName, "_", genotype, "_f_", ageClass, "_", genotypeNames[j], "_", k))
              #print("in8")
              #print(dim(transitions))
              #print(ind)
              #print(counter)
              transitions[ind, counter] <- -1
              #print("out8")
              transitions[matedInd, counter] <- 1
              counter <- counter + 1
            }
          }
        }
      }
      else
      {
        #Mated and can birth new individuals
        transitionRates_vector[counter] <- state_vector[ind]*(birthRateByName_vector[stateNames_vector[ind]])*(1 - ciByName_vector[stateNames_vector[ind]])*(ImaxByPopulation_vector[populationName] - ItotalByPopulations[[populationName]])/ImaxByPopulation_vector[populationName]
        if(transitionRates_vector[counter] < 0){
          transitionRates_vector[counter] <- 0
        }
        birthInd <- which(stateNames_vector == paste0(populationName, "_", offspringGenotypeByName_vector[stateNames_vector[ind]], "_imm_1")) 
        transitions[birthInd, counter] <- 1  
        names(transitionRates_vector)[counter] <- paste0("birth_", stateNames_vector[ind])
        counter <- counter + 1
        
      }
    }
  }
  
  
  # Immigration
  # the rate of immigration is not proportional to state size
  # This means there is no steady state possible. Not true: steady state depends on relative immigration to emigration ratio. Equal means steady state as expected. Immigration>emigration means exponential growth to infinity. Immigration<emigration means tends to zero
  immigrationInds <- which(immigrationRate_vector > 0)
  names(immigrationRate_vector) <- stateNames_vector
  if(length(immigrationInds) > 0)
  {
    for(i in 1:length(immigrationInds))
    {
      immStateInd <- which(stateNames_vector == names(immigrationRate_vector)[immigrationInds[i]])
      transitionRates_vector[counter] <- immigrationRate_vector[immigrationInds[i]]
      transitions[immStateInd, counter] <- 1
      counter <- counter + 1
    }
  }
  
  # Emigration
  # the rate of emigration is  proportional to state size
  emigrationInds <- which(emigrationRate_vector > 0)
  names(emigrationRate_vector) <- stateNames_vector
  if(length(emigrationInds) > 0)
  {
    for(i in 1:length(emigrationInds))
    {
      emmStateInd <- which(stateNames_vector == names(emigrationRate_vector)[emigrationInds[i]])
      transitionRates_vector[counter] <- emigrationRate_vector[emigrationInds[i]] * state_vector[emmStateInd]
      transitions[emmStateInd, counter] <- -1
      counter <- counter + 1
    }
  }
  
  transitionRates_vector <- transitionRates_vector[1:(counter - 1)]
  if((counter - 1) == 1)
  {
    transitions <- matrix(transitions[, 1:(counter - 1)], length(state_vector), (counter - 1))
  }
  else
  {
    transitions <- transitions[, 1:(counter - 1)]
  }
  # stop(glue::glue("Counter is {counter}"))
  
  return(list(transitionRates = transitionRates_vector, transitions = transitions))
}



runModel <- function(state_vector, startTime, endTime, interventionFunction, interventionTimes, recordTimes,
                     numClassesBySexAndGenotype_matrix, matingRateByName_vector, ciByName_vector, offspringGenotypeByName_vector, 
                     friedsByName_vector, deathRatesByName_vector, birthRateByName_vector, transitionRateByName_vector, emmigrationRate_matrix, 
                     immigrationRate_vector, emigrationRate_vector, carryingCapacityByPopulation_vector, ImaxByPopulation_vector,
                     progress=FALSE, cage = FALSE, stop_intervention = FALSE)
{
  currentTime <- startTime
  interventionTimeIndex <- 1
  recordTimeIndex <- 1
  
  
  # These indices won't work for generalised code
  m_f_ind <- getIndicesForStateQuery(stateNames = names(state_vector),
                                     sex = c("m", "f"), genotype = "wAlbAB")
  init_pop <- sum(state_vector[m_f_ind])
  
  
  ## Create stopping conditions for the simulation
  ## Currently consider:
  ### Cage simulations: stop if either population completely takes over
  ### IIT interventions: stop if suppression reached
  #### There will be other stopping conditions that can be met depending on the intervention, i.e.
  ##### Complete stop: if wPip go above threshold, stop intervention and simulation
  stop_sims <- FALSE
  stop_counter <- -1
  if(isTRUE(cage)){
    # Flags to navigate simulations
    first_run <- TRUE
    switch_props <- FALSE
    stopping_condition <- function(...){
      ARGS <- list(...)
      state_vector <- ARGS$state_vector
      init_pop <- ARGS$init_pop[1] # Hardcoded for one population
      
      # Stop if make up adult population
      alb_ind <- getIndicesForStateQuery(stateNames = names(state_vector),
                                         genotype = "wAlbAB", sex = c("m", "f"))
      # pip_ind <- getIndicesForStateQuery(stateNames = names(state_vector),
      #                                    genotype = "wPip", sex = c("m", "f"))
      # 
      # Stop if make up full population
      # alb_ind <- getIndicesForStateQuery(stateNames = names(state_vector),
      #                                    genotype = "wAlbAB")
      pip_ind <- getIndicesForStateQuery(stateNames = names(state_vector),
                                         genotype = "wPip")
      
      
      # Stop if either 
      ## wAlbAB is suppressed (Adult population less than 10% of population)
      ## Total wPip including immatures is 0 (wPip extinct)
      prop_wAlb <- sum(state_vector[alb_ind])/init_pop
      total_wPip <- sum(state_vector[pip_ind])
      
      if(isTRUE(first_run)){
        if(prop_wAlb < 0.5){
          switch_props <<- TRUE
        }
        first_run <<- FALSE
      }
      
      if(isTRUE(switch_props)){ # If start with more wPip than wAlb
        # Stop if make up adult population
        # alb_ind <- getIndicesForStateQuery(stateNames = names(state_vector),
        #                                    genotype = "wAlbAB", sex = c("m", "f"))
        pip_ind <- getIndicesForStateQuery(stateNames = names(state_vector),
                                           genotype = "wPip", sex = c("m", "f"))
        # 
        # Stop if make up full population
        alb_ind <- getIndicesForStateQuery(stateNames = names(state_vector),
                                           genotype = "wAlbAB")
        # pip_ind <- getIndicesForStateQuery(stateNames = names(state_vector),
        #                                    genotype = "wPip")
        
        
        # Stop if either 
        ## wAlbAB is suppressed (Adult population less than 10% of population)
        ## Total wPip including immatures is 0 (wPip extinct)
        prop_wPip <- sum(state_vector[pip_ind])/init_pop
        total_wAlb <- sum(state_vector[alb_ind])
        
        if((total_wAlb==0) | (prop_wPip<0.1)){
          return(list(stop_sims=TRUE, stop_counter=0, stop_intervention=TRUE))
        }else{
          return(list(stop_sims=FALSE, stop_counter=-1))
        }
      }
      
      if((total_wPip==0) | (prop_wAlb<0.1)){
        return(list(stop_sims=TRUE, stop_counter=0, stop_intervention=TRUE))
      }else{
        return(list(stop_sims=FALSE, stop_counter=-1))
      }
    }
  }else{
    stopping_condition <- function(...){ # MM: need edits
      ARGS <- list(...)
      
      if(isTRUE(ARGS$currentTime > 365*2)){
        return(list(timeout = TRUE))
      }else{
        return(list(timeout = FALSE))
      }
      
      #   state_vector <- ARGS$state_vector
      #   init_pop <- ARGS$init_pop
      #   alb_ind <- getIndicesForStateQuery(stateNames = names(state_vector),
      #                                      genotype = "wAlbAB", sex = c("m", "f"))
      #   
      #   prop_wAlb <- sum(state_vector[alb_ind])/init_pop
      # 
      #   
      #   # Stop if population drops below 10% of original
      #   # This is the "suppression" level
      #   if(prop_wAlb < 0.1){ 
      #     return(list(stop_sims=TRUE, stop_counter=0, stop_intervention=TRUE))
      #   }else{
      #     return(list(stop_sims=FALSE, stop_counter=-1))
      #   }
      # }
    }
  }
  
  storage <- matrix(0, length(state_vector) + 1, length(recordTimes) + 1)
  # Record initial state vector
  save_vector <- c(state_vector, "released_to_date"=0)
  storage[, recordTimeIndex] <- save_vector
  recordTimeIndex <- recordTimeIndex + 1
  
  # Create a progress bar
  tick <- 0
  if(progress){
    pd <- progress_bar$new(format = "  Run time [:bar] :current/:total (:percent), ellapsed: :elapsed eta: :eta",
                           total = end_time, clear = FALSE)
    pd$tick(0)
    update <- function(currentTime, tick){
      if(floor(currentTime) >= tick){ # Increase progress bar every time step, not intermediate time
        # print(currentTime)
        tick <<- tick + 1 # advance tick counter
        pd$tick(1)
        # return(pd)
      }
    }
  }else{
    update <- function(currentTime, tick){
      NULL
    }
  }
  
  
  # if(is.null(stoppingCondition_function)){
  #   stoppingCondition_function <- function(...){
  #     NULL
  #   }
  # }else{
  #   ARGS <- purrr::map(stopping_args,
  #                      function(s){
  #                        sym <- as.symbol(s)
  #                        return(sym)
  #                      })
  #   names(ARGS) <- as.character(stopping_args)
  # }
  
  total_released <- 0
  # num_events <- 0
  timeout <- FALSE
  while(currentTime < endTime)
  {
    update(currentTime, tick)
    
    # if(isTRUE(stoppingCondition_function(purrr::map(ARGS, eval)))){
    #   break
    # }
    if(all(state_vector == 0))
    {
      break  
    }
    
    
    if(isFALSE(stop_sims)){ # MM: need edits
      tmp <- stopping_condition(state_vector = state_vector, init_pop = init_pop, currentTime=currentTime)
      if(!is.null(tmp$stop_sims)){
        stop_sims <- tmp$stop_sims
      }
      if(!is.null(tmp$stop_counter)){
        stop_counter <- tmp$stop_counter 
      }
      if(!is.null(tmp$stop_intervention)){
        stop_intervention <- tmp$stop_intervention
      }
      if(!is.null(tmp$timeout)){
        timeout <- tmp$timeout
      }
      # print(stop_sims)
    }
    # if(isTRUE(sum(state_vector[m_f_ind])/init_pop < 0.1)){
    #   break
    # }
    transitions <- getTransitionRates(state_vector, numClassesBySexAndGenotype_matrix, matingRateByName_vector, ciByName_vector, offspringGenotypeByName_vector, 
                                      friedsByName_vector, deathRatesByName_vector, birthRateByName_vector, transitionRateByName_vector, emmigrationRate_matrix, 
                                      immigrationRate_vector, emigrationRate_vector, carryingCapacityByPopulation_vector, ImaxByPopulation_vector)
    
    rates <- transitions$transitionRates
    
    rand_transitionIndex <- tryCatch(sample(1:length(rates), 1, prob = rates/sum(rates)),
                                     error = function(e) print(transitions$transitionRates))
    rand_time <- rexp(1, rate = sum(rates))
    proposedTime <- currentTime + rand_time
    
    if(proposedTime > interventionTimes[interventionTimeIndex] & interventionTimeIndex <= length(interventionTimes))
    {
      intervention <- interventionFunction(state_vector, startTime, endTime, interventionTimes[interventionTimeIndex], 
                                           stop_intervention=stop_intervention, timeout=timeout)
      state_vector <- intervention$state_vector
      if(!is.null(intervention$stop_intervention)){
        stop_intervention <- intervention$stop_intervention
      }
      currentTime <- interventionTimes[interventionTimeIndex]
      proposedTime <- currentTime
      interventionTimeIndex <- interventionTimeIndex + 1
      total_released <- total_released + intervention$numReleased
    }
    if(proposedTime > recordTimes[recordTimeIndex] & recordTimeIndex  <= length(recordTimes))
    {
      save_vector <- c(state_vector, "released_to_date"=total_released)
      storage[, recordTimeIndex] <- save_vector
      recordTimeIndex <- recordTimeIndex + 1
      if(isTRUE(stop_counter>=0)){
        stop_counter <- stop_counter + 1
      }
    }
    if(proposedTime > currentTime)
    {
      currentTime <- proposedTime
      transitionVector <- transitions$transitions[, rand_transitionIndex]
      state_vector_new <- round(state_vector + transitionVector) #unclear why i am forced to round here.  R seems to mess up integer addition here.
      state_vector <- state_vector_new
    }
    # Turn on to stop once wAlbAB eliminated
    # wAlb_ind <- getIndicesForStateQuery(names(state_vector), genotype = "wAlbAB")
    # wAlb_mates_ind <- getIndicesForStateQuery(names(state_vector), genotype = "wPip", mate = "wAlbAB")
    # if(sum(state_vector[c(wAlb_mates_ind, wAlb_ind)]) == 0){
    #   break
    # }
    # num_events <- num_events + 1
    
    # stop 6 months after sufficient stopping condition reached.
    if(isTRUE(stop_counter >= 180)){ 
      break
    }
  }
  
  stop_intervention <<- FALSE
  # print(num_events)
  return(storage)
}



# Intervention function ---------------------------------------------------

#HERE YOU CAN WRITE A FUNCTION THAT RELEASES MOSQUITOES AT SPECIFIED TIMES TO TEST A SCENARIO
#Currently set to not release any mosquitoes
# This seems to be a function designed to dump in wPip mosquitos?
# Also, this isn't releasing eggs?  This is releasing age 1 adult males?
# count number of adult male AblAB
# Release immature wPip: NO
# Should not actually impact the breeding limits
# The should be released at age 1
# Boxes released, so "release date" is on the maturation date of mozzies
interventionFunction <- function(state_vector, startTime, endTime, interventionTime)
{
  # return(state_vector)
  overfloodingDump <- 1
  femaleContamination <- 0.01 # Change this for the contamination rate
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




## State vector for cage experiments
create_cage_state_vector <- function(state_vector, propPip, total_capacity = 400, female_split = 0.5){
  
  stateNames_vector <- names(state_vector)
  cage_state <- rep(0, length(stateNames_vector))
  names(cage_state) <- stateNames_vector
  
  num_males <- round((1-female_split) * total_capacity)
  num_females <- round(female_split * total_capacity)
  
  m_wAlbAB_ind <- getIndicesForStateQuery(stateNames = stateNames_vector, genotype = "wAlbAB", 
                                          sex = "m", ageClass = 1)
  f_wAlbAB_ind <- getIndicesForStateQuery(stateNames = stateNames_vector, genotype = "wAlbAB", 
                                          sex = "f", ageClass = 1, mate = "Unmated")
  m_wPip_ind <- getIndicesForStateQuery(stateNames = stateNames_vector, genotype = "wPip", 
                                        sex = "m", ageClass = 1)
  f_wPip_ind <- getIndicesForStateQuery(stateNames = stateNames_vector, genotype = "wPip", 
                                        sex = "f", ageClass = 1, mate = "Unmated")
  
  cage_state[m_wAlbAB_ind] <- round(num_males * (1-propPip))
  cage_state[f_wAlbAB_ind] <- round(num_females * (1-propPip))
  cage_state[m_wPip_ind] <- round(num_males * (propPip))
  cage_state[f_wPip_ind] <- round(num_females * (propPip))
  
  return(cage_state)
  
}

## Helper function for immigration values
get_age_probability <- function(i, K = numMaleClasses[1],
                                muD = maleDeathRate[1], tau = transitionRate_m[1]){
  # Calculate the probability of getting to a certain age before dying
  if(i==K){
    1 - pexp((1/tau)*(i-1), rate = muD)
  }else{
    pexp((1/tau)*(i), rate = muD) - pexp((1/tau)*(i-1), rate=muD) 
  }
}

## Immigration vector
create_immigration_emmigration_vectors <- function(state_vector, overallMaleImmigrationRate,
                                                   overallFemaleImmigrationRate, maleDeathRate,
                                                   numMaleClasses, phi, theta, M1, carryingCapacityByPopulation_vector,
                                                   transitionRate_m = 1){
  # Initialise vectors
  immigrationRate_vector <- rep(0, length(state_vector))
  names(immigrationRate_vector) <- names(state_vector)
  emigrationRate_vector <- rep(0, length(state_vector))
  names(emigrationRate_vector) <- names(state_vector)
  
  # Define quantities
  # immature_ind <- getIndicesForStateQuery(stateNames = names(state_vector), sex = "imm")
  C <- carryingCapacityByPopulation_vector[1] # Assuming a single block population#sum(state_vector[-immature_ind]) # Carrying capacity, assumes state vector is steady state of wild type
  K <- numMaleClasses[1]
  
  m_ind <- getIndicesForStateQuery(stateNames = names(state_vector), genotype = "wAlbAB", sex = "m")
  m_ind_pip <- getIndicesForStateQuery(stateNames = names(state_vector), genotype = "wPip", sex = "m")
  
  f_ind <- getIndicesForStateQuery(stateNames = names(state_vector), genotype = "wAlbAB", sex = "f",
                                   mate = c("Unmated", "wAlbAB"))
  f_ind_total <- getIndicesForStateQuery(stateNames = names(state_vector), sex = "f")
  
  
  # Indices are hardcoded to the wild type class to get steady state values
  # Wild-type is assumed to be the first genotype.
  
  j <- 1:numMaleClasses[1]
  Pj <- lapply(j, get_age_probability,  K = numMaleClasses[1],
               muD = maleDeathRate[1], tau = transitionRate_m[1]) |>
    unlist()
  
  # Male immigration and emigration
  
  M_denominator <- M1[1] * (sum(Pj[1:(K-1)] * phi[1]^((1:(K-1)) - 1)) + transitionRate_m*(phi[1]^(K-2)) * Pj[K] / maleDeathRate[1] )
  
  xi_j_M <- (overallMaleImmigrationRate * Pj) / M_denominator
  
  immigrationRate_vector[m_ind] <- xi_j_M * state_vector[m_ind]
  emigrationRate_vector[m_ind] <- emigrationRate_vector[m_ind_pip] <-  xi_j_M 
  
  # Female immigration and emigration
  
  xi_F <- ((1 + theta[1]) * overallFemaleImmigrationRate) / (theta[1] * C)
  
  emigrationRate_vector[f_ind_total] <- xi_F
  immigrationRate_vector[f_ind] <- xi_F * state_vector[f_ind]
  
  
  return(list(immigrationRate_vector = immigrationRate_vector, emigrationRate_vector = emigrationRate_vector))
}