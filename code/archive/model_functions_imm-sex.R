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
  # For each popluation and genotpye
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
      theseInds_immM <- getIndicesForStateQuery(stateNames, populationName = populationNames[p], 
                                                genotype = genotypeNames[g], sex = "immM", 
                                                ageClass = NA, mate = NA)
      theseInds_immF <- getIndicesForStateQuery(stateNames, populationName = populationNames[p],
                                                genotype = genotypeNames[g], sex = "immF", 
                                                ageClass = NA, mate = NA)
      # todo: ensure there are no "_imm" whole words left -- must have either "_m" or "_f" appended
      # Add along all model outputs for particular genotype population, and age
      if(length(theseInds_m) == 1){
        timeSeries[[paste0(populationNames[p], "_", genotypeNames[g], "_m")]] <- modelOutputs[theseInds_m, ]
      }else{
        timeSeries[[paste0(populationNames[p], "_", genotypeNames[g], "_m")]] <- apply(modelOutputs[theseInds_m, ], 2, sum) 
      }
      timeSeries[[paste0(populationNames[p], "_", genotypeNames[g], "_f")]] <- apply(modelOutputs[theseInds_f, ], 2, sum)
      timeSeries[[paste0(populationNames[p], "_", genotypeNames[g], "_immM")]] <- apply(modelOutputs[theseInds_immM, ], 2, sum)
      timeSeries[[paste0(populationNames[p], "_", genotypeNames[g], "_immF")]] <- apply(modelOutputs[theseInds_immF, ], 2, sum)
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
      stateNames <- c(stateNames, paste0(populationNames[p], "_", genotypeNames[g], "_m_", 1:n_m[g]), paste0(populationNames[p], "_", genotypeNames[g], "_f_", 1:n_f[g], "_Unmated_None"), paste0(populationNames[p], "_", genotypeNames[g], "_immM_", 1:n_i[g]), paste0(populationNames[p], "_", genotypeNames[g], "_immF_", 1:n_i[g]))
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
                               emmigrationRate_matrix, immigrationRate_vector, carryingCapacityByPopulation_vector, 
                               ImaxByPopulation_vector)
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
  
  #Possible transition rates within each population are:
  #Each male class can die ($\sum_{i = 1}^g n_{m, i}$)
  #Each male class except the last one can transition ($\sum_{i = 1}^g (n_{m, i} - 1)$))
  #Each female class can die ($\sum_{i = 1}^g n_{f, i}$)
  #Each female class except the last one can transition ($\sum_{i = 1}^g (n_{f, i} - 1)$))
  #Each immature genotype can have a birth into the first class ($n_g$)
  #Each immature class except the last one can transition into another one ($\sum_{j = 1}^g (n_{i, j} - 1)$))
  #The last class in each immature genotype can birth either a male or a female (2 x n_g)
  
  # The following is a bunch of counting to see how many possible transitions can occur.
  stateNames_vector <- names(state_vector)
  
  populationNames <- names(carryingCapacityByPopulation_vector)
  numPopulations <- length(carryingCapacityByPopulation_vector) #query: This is the only place carring capacity is brough in, and that is only used to get the number of populations
  
  genotypeNames <- rownames(numClassesBySexAndGenotype_matrix)
  
  totalMaleDeathTransitions <- sum(numClassesBySexAndGenotype_matrix[, "m"])
  totalMaleClassTransitions <- sum(numClassesBySexAndGenotype_matrix[, "m"] - 1)
  totalFemaleDeathTransitions <- sum(numClassesBySexAndGenotype_matrix[, "f"])
  totalFemaleClassTransitions <- sum(numClassesBySexAndGenotype_matrix[, "f"] - 1)*sum(numClassesBySexAndGenotype_matrix[, "m"]) #query: Why the number of males multiplied?
  totalImmatureBirthTransitions <- sum(numClassesBySexAndGenotype_matrix[, "f"])*length(genotypeNames) #query: Why is this multplying by length(genotypeNames)
  totalImmatureMaleClassTransitions <- sum(numClassesBySexAndGenotype_matrix[, "immM"] - 1)
  totalImmatureFemaleClassTransitions <- sum(numClassesBySexAndGenotype_matrix[, "immF"] - 1)
  totalAdultBirthTransitions <- length(genotypeNames) #query: Should this be multiplied by 2?
  
  totalMaleImmigrationTransitions <- sum(numClassesBySexAndGenotype_matrix[, "m"])
  totalFemaleImmigrationTransitions <- sum(numClassesBySexAndGenotype_matrix[, "f"])
  totalPossibleMaleEmmigrations <- sum(numClassesBySexAndGenotype_matrix[, "m"])*length(which(emmigrationRate_matrix > 0))
  totalPossibleFemaleEmmigrations <- sum(numClassesBySexAndGenotype_matrix[, "f"])*length(which(emmigrationRate_matrix > 0))
  
  #query: Why is this multiplied by 5?
  n_transitions <- 5*numPopulations*(totalMaleImmigrationTransitions + totalFemaleImmigrationTransitions + totalPossibleMaleEmmigrations + totalPossibleFemaleEmmigrations + totalMaleDeathTransitions + totalMaleClassTransitions + totalFemaleDeathTransitions + totalFemaleClassTransitions +  totalImmatureBirthTransitions + totalImmatureMaleClassTransitions + totalImmatureFemaleClassTransitions + totalAdultBirthTransitions)
  
  transitionRates_vector <- rep(0, n_transitions)
  transitions <- matrix(0, length(state_vector), n_transitions)
  
  maleIndices <- which(grepl( "_m_", stateNames_vector, fixed = TRUE) & state_vector > 0)
  femaleIndices <- which(grepl( "_f_", stateNames_vector, fixed = TRUE) & state_vector > 0)
  immatureMaleIndices <- which(grepl( "_immM_", stateNames_vector, fixed = TRUE) & state_vector > 0)
  immatureFemaleIndices <- which(grepl( "_immF_", stateNames_vector, fixed = TRUE) & state_vector > 0)
  
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
  
  
  
  #All transitions for male immature states with at least one individual
  #Need to count the number of immatures in each population
  ItotalByPopulations <- rep(0, length(populationNames))
  names(ItotalByPopulations) <- populationNames
  if(length(immatureMaleIndices) > 0)
  {
    for(i in 1:length(immatureMaleIndices))
    {
      ind <- immatureMaleIndices[i]
      bits <- unlist(strsplit(stateNames_vector[ind], "_"))
      populationName <- bits[1]
      genotype <- bits[2]
      ageClass <- as.numeric(bits[4])
      
      ItotalByPopulations[[populationName]] <- ItotalByPopulations[[populationName]] + state_vector[ind]
      
      if(ageClass < numClassesBySexAndGenotype_matrix[genotype, "immM"])
      {
        nextInd <- which(stateNames_vector == paste0(populationName, "_", genotype, "_immM_", ageClass + 1))
        transitionRates_vector[counter] <- state_vector[ind]*transitionRateByName_vector[stateNames_vector[ind]]
        transitions[ind, counter] <- -1
        transitions[nextInd, counter] <- 1
        counter <- counter + 1
      }
      if(ageClass == numClassesBySexAndGenotype_matrix[genotype, "immM"])
      {
        nextInd_m <- which(stateNames_vector == paste0(populationName, "_", genotype, "_m_", 1))
        transitionRates_vector[counter] <- state_vector[ind]*transitionRateByName_vector[stateNames_vector[ind]]
        transitions[ind, counter] <- -1
        transitions[nextInd_m, counter] <- 1
        counter <- counter + 1
      }
    }
  }
  #All transitions for female immature states with at least one individual
  #Need to count the number of immatures in each population
  # edit:
  # ItotalByPopulations <- rep(0, length(populationNames))
  # names(ItotalByPopulations) <- populationNames
  if(length(immatureFemaleIndices) > 0)
  {
    for(i in 1:length(immatureFemaleIndices))
    {
      ind <- immatureFemaleIndices[i]
      bits <- unlist(strsplit(stateNames_vector[ind], "_"))
      populationName <- bits[1]
      genotype <- bits[2]
      ageClass <- as.numeric(bits[4])
      
      ItotalByPopulations[[populationName]] <- ItotalByPopulations[[populationName]] + state_vector[ind]
      
      if(ageClass < numClassesBySexAndGenotype_matrix[genotype, "immF"])
      {
        nextInd <- which(stateNames_vector == paste0(populationName, "_", genotype, "_immF_", ageClass + 1))
        transitionRates_vector[counter] <- state_vector[ind]*transitionRateByName_vector[stateNames_vector[ind]]
        transitions[ind, counter] <- -1
        transitions[nextInd, counter] <- 1
        counter <- counter + 1
      }
      if(ageClass == numClassesBySexAndGenotype_matrix[genotype, "immF"])
      {
        nextInd_f <- which(stateNames_vector == paste0(populationName, "_", genotype, "_f_", 1, "_Unmated_None"))
        transitionRates_vector[counter] <- state_vector[ind]*transitionRateByName_vector[stateNames_vector[ind]]
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
        # Male births
        transitionRates_vector[counter] <- state_vector[ind]*(birthRateByName_vector[stateNames_vector[ind]]/2)*(1 - ciByName_vector[stateNames_vector[ind]])*(ImaxByPopulation_vector[populationName] - ItotalByPopulations[[populationName]])/ImaxByPopulation_vector[populationName]
        if(transitionRates_vector[counter] < 0){
          transitionRates_vector[counter] <- 0
        }
        birthInd <- which(stateNames_vector == paste0(populationName, "_", offspringGenotypeByName_vector[stateNames_vector[ind]], "_immM_1"))
        transitions[birthInd, counter] <- 1  
        names(transitionRates_vector)[counter] <- paste0("birth_male_", stateNames_vector[ind])
        counter <- counter + 1
        
        # female births
        transitionRates_vector[counter] <- state_vector[ind]*(birthRateByName_vector[stateNames_vector[ind]]/2)*(1 - ciByName_vector[stateNames_vector[ind]])*(ImaxByPopulation_vector[populationName] - ItotalByPopulations[[populationName]])/ImaxByPopulation_vector[populationName]
        if(transitionRates_vector[counter] < 0){
          transitionRates_vector[counter] <- 0
        }
        birthInd <- which(stateNames_vector == paste0(populationName, "_", offspringGenotypeByName_vector[stateNames_vector[ind]], "_immF_1"))
        transitions[birthInd, counter] <- 1  
        names(transitionRates_vector)[counter] <- paste0("birth_female_", stateNames_vector[ind])
        counter <- counter + 1
        # #Mated and can birth new individuals
        # # print("Imax")
        # # print(ImaxByPopulation_vector[populationName])
        # # print("Itotal")
        # # print(ItotalByPopulations[[populationName]])
        # transitionRates_vector[counter] <- state_vector[ind]*birthRateByName_vector[stateNames_vector[ind]]*(1 - ciByName_vector[stateNames_vector[ind]])*(ImaxByPopulation_vector[populationName] - ItotalByPopulations[[populationName]])/ImaxByPopulation_vector[populationName]
        # if (runif(1)<0.5) {
        #   birthInd <- which(stateNames_vector == paste0(populationName, "_", offspringGenotypeByName_vector[stateNames_vector[ind]], "_immM_1"))
        # } else {
        #   birthInd <- which(stateNames_vector == paste0(populationName, "_", offspringGenotypeByName_vector[stateNames_vector[ind]], "_immF_1"))
        # }
        # # todo: double check this is for single birth process at a time, otherwise need a 0.5 somewhere
        # # ans: I think this is for a single process at a time.  From what I gather, this is asking for 
        # # a particular mated female whether they can birth and male or a female offspring.
        # # The transition matrix is generated every loop of the model I believe, so this will randomly
        # # switch between generating males and females.
        # transitions[birthInd, counter] <- 1  
        # counter <- counter + 1
      }
    }
  }
  
  
  
  
  
  # Immigration
  # Immigration is coded wrong
  immigrationInds <- which(immigrationRate_vector > 0)
  names(immigrationRate_vector) <- stateNames_vector
  if(length(immigrationInds) > 0)
  {
    for(i in 1:length(immigrationInds))
    {
      emmStateInd <- which(stateNames_vector == names(immigrationRate_vector)[immigrationInds[i]])
      transitionRates_vector[counter] <- immigrationRate_vector[immigrationInds[i]]
      transitions[emmStateInd, counter] <- 1
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
  
  
  return(list(transitionRates = transitionRates_vector, transitions = transitions))
}



runModel <- function(state_vector, startTime, endTime, interventionFunction, interventionTimes, recordTimes,
                     numClassesBySexAndGenotype_matrix, matingRateByName_vector, ciByName_vector, offspringGenotypeByName_vector, 
                     friedsByName_vector, deathRatesByName_vector, birthRateByName_vector, transitionRateByName_vector, emmigrationRate_matrix, 
                     immigrationRate_vector, carryingCapacityByPopulation_vector, ImaxByPopulation_vector)
{
  currentTime <- startTime
  interventionTimeIndex <- 1
  recordTimeIndex <- 1
  
  # These indices won't work for generalised code
  m_f_ind <- getIndicesForStateQuery(stateNames = names(state_vector),sex = c("m", "f"), genotype = "wAlbAB")
  init_pop <- sum(state_vector[m_f_ind])
  
  storage <- matrix(0, length(state_vector) + 1, length(recordTimes))
  # Create a progress bar
  pd <- progress_bar$new(format = "  Run time [:bar] :current/:total (:percent), ellapsed: :elapsed eta: :eta",
                         total = end_time, clear = FALSE)
  pd$tick(0)
  tick <- 0
  total_released <- 0
  while(currentTime < endTime)
  {
    if(floor(currentTime) >= tick){ # Increase progress bar every time step, not intermediate time
      # print(currentTime)
      pd$tick()
      tick <- tick + 1 # advance tick counter
    }
    if(all(state_vector == 0))
    {
      break  
    }
    # if(isTRUE(sum(state_vector[m_f_ind])/init_pop < 0.1)){
    #   break
    # }
    transitions <- getTransitionRates(state_vector, numClassesBySexAndGenotype_matrix, matingRateByName_vector, ciByName_vector, offspringGenotypeByName_vector, 
                                      friedsByName_vector, deathRatesByName_vector, birthRateByName_vector, transitionRateByName_vector, emmigrationRate_matrix, 
                                      immigrationRate_vector, carryingCapacityByPopulation_vector, ImaxByPopulation_vector)
    
    rates <- transitions$transitionRates
    
    rand_transitionIndex <- tryCatch(sample(1:length(rates), 1, prob = rates/sum(rates)),
                                     error = function(e) print(transitions$transitionRates))
    rand_time <- rexp(1, rate = sum(rates))
    proposedTime <- currentTime + rand_time
    
    if(proposedTime > interventionTimes[interventionTimeIndex] & interventionTimeIndex <= length(interventionTimes))
    {
      intervention <- interventionFunction(state_vector, startTime, endTime, interventionTimes[interventionTimeIndex])
      state_vector <- intervention$state_vector
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
  }
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

