
getIndicesForStateQuery <- function(stateNames, populationName = NA, genotype = NA, sex = NA, ageClass = NA, mate = NA, mateStage = NA)
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
    for(p in 1:length(populationNames))
    {
      for(g in 1:length(genotypeNames))
      {
        theseInds_m <- getIndicesForStateQuery(stateNames, populationName = populationNames[p], genotype = genotypeNames[g], sex = "m", ageClass = NA, mate = NA)
        theseInds_f <- getIndicesForStateQuery(stateNames, populationName = populationNames[p], genotype = genotypeNames[g], sex = "f", ageClass = NA, mate = NA)
        theseInds_immM <- getIndicesForStateQuery(stateNames, populationName = populationNames[p], genotype = genotypeNames[g], sex = "immM", ageClass = NA, mate = NA)
        theseInds_immF <- getIndicesForStateQuery(stateNames, populationName = populationNames[p], genotype = genotypeNames[g], sex = "immF", ageClass = NA, mate = NA)
        # todo: ensure there are no "_imm" whole words left -- must have either "_m" or "_f" appended
        timeSeries[[paste0(populationNames[p], "_", genotypeNames[g], "_m")]] <- apply(modelOutputs[theseInds_m, ], 2, sum)
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


getTransitionRates <- function(state_vector, numClassesBySexAndGenotype_matrix, matingRateByName_vector, ciByName_vector, offspringGenotypeByName_vector, friedsByName_vector, deathRatesByName_vector, birthRateByName_vector, transitionRateByName_vector, emmigrationRate_matrix, immigrationRate_vector, carryingCapacityByPopulation_vector, ImaxByPopulation_vector)
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

  
  stateNames_vector <- names(state_vector)
  
  populationNames <- names(carryingCapacityByPopulation_vector)
  numPopulations <- length(carryingCapacityByPopulation_vector)
  
  genotypeNames <- rownames(numClassesBySexAndGenotype_matrix)
  
  totalMaleDeathTransitions <- sum(numClassesBySexAndGenotype_matrix[, "m"])
  totalMaleClassTransitions <- sum(numClassesBySexAndGenotype_matrix[, "m"] - 1)
  totalFemaleDeathTransitions <- sum(numClassesBySexAndGenotype_matrix[, "f"])
  totalFemaleClassTransitions <- sum(numClassesBySexAndGenotype_matrix[, "f"] - 1)*sum(numClassesBySexAndGenotype_matrix[, "m"])
  totalImmatureBirthTransitions <- sum(numClassesBySexAndGenotype_matrix[, "f"])*length(genotypeNames)
  totalImmatureMaleClassTransitions <- sum(numClassesBySexAndGenotype_matrix[, "immM"] - 1)
  totalImmatureFemaleClassTransitions <- sum(numClassesBySexAndGenotype_matrix[, "immF"] - 1)
  totalAdultBirthTransitions <- length(genotypeNames)
  
  totalMaleImmigrationTransitions <- sum(numClassesBySexAndGenotype_matrix[, "m"])
  totalFemaleImmigrationTransitions <- sum(numClassesBySexAndGenotype_matrix[, "f"])
  totalPossibleMaleEmmigrations <- sum(numClassesBySexAndGenotype_matrix[, "m"])*length(which(emmigrationRate_matrix > 0))
  totalPossibleFemaleEmmigrations <- sum(numClassesBySexAndGenotype_matrix[, "f"])*length(which(emmigrationRate_matrix > 0))
  
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
  ItotalByPopulations <- rep(0, length(populationNames))
  names(ItotalByPopulations) <- populationNames
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
        transitionRates_vector[counter] <- state_vector[ind]*birthRateByName_vector[stateNames_vector[ind]]*(1 - ciByName_vector[stateNames_vector[ind]])*(ImaxByPopulation_vector[populationName] - ItotalByPopulations[[populationName]])/ImaxByPopulation_vector[populationName]
        if (runif(1)<0.5) {
          birthInd <- which(stateNames_vector == paste0(populationName, "_", offspringGenotypeByName_vector[stateNames_vector[ind]], "_immM_1"))
        } else {
          birthInd <- which(stateNames_vector == paste0(populationName, "_", offspringGenotypeByName_vector[stateNames_vector[ind]], "_immF_1"))
        }
        transitions[birthInd, counter] <- 1  # todo: double check this is for single birth process at a time, otherwise need a 0.5 somewhere
        counter <- counter + 1
      }
    }
  }
  
  

  
  
  #Immigration
  immigrationInds <- which(immigrationRate_vector > 0)
  if(length(immigrationInds) > 0)
  {
    for(i in 1:length(immigrationInds))
    {
      emmStateInd <- which(stateNames_vector == names(immigrationRate_vector)[immigrationInds])
      transitionRates_vector[counter] <- immigrationRate_vector[immigrationInds[i]]
      transitions[emmStateInd, counter] <- 1
      counter <- counter + 1
    }
  }

  transitionRates_vector <- transitionRates_vector[1:(counter - 1)]
  if(counter - 1 == 1)
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
  
  storage <- matrix(0, length(state_vector), length(recordTimes))
  
  while(currentTime < endTime)
  {
    print(currentTime)
    if(all(state_vector == 0))
    {
      break  
    }
    transitions <- getTransitionRates(state_vector, numClassesBySexAndGenotype_matrix, matingRateByName_vector, ciByName_vector, offspringGenotypeByName_vector, 
                                      friedsByName_vector, deathRatesByName_vector, birthRateByName_vector, transitionRateByName_vector, emmigrationRate_matrix, 
                                      immigrationRate_vector, carryingCapacityByPopulation_vector, ImaxByPopulation_vector)
    
    rates <- transitions$transitionRates
    
    rand_transitionIndex <- sample(1:length(rates), 1, prob = rates/sum(rates))
    rand_time <- rexp(1, rate = sum(rates))
    proposedTime <- currentTime + rand_time

    if(proposedTime > interventionTimes[interventionTimeIndex] & interventionTimeIndex <= length(interventionTimes))
    {
      state_vector <- interventionFunction(state_vector, startTime, endTime, interventionTimes[interventionTimeIndex])
      currentTime <- interventionTimes[interventionTimeIndex]
      proposedTime <- currentTime
      interventionTimeIndex <- interventionTimeIndex + 1
    }
    if(proposedTime > recordTimes[recordTimeIndex] & recordTimeIndex  <= length(recordTimes))
    {
      storage[, recordTimeIndex] <- state_vector
      recordTimeIndex <- recordTimeIndex + 1
    }
    if(proposedTime > currentTime)
    {
      currentTime <- proposedTime
      transitionVector <- transitions$transitions[, rand_transitionIndex]
      state_vector_new <- round(state_vector + transitionVector) #unclear why i am forced to round here.  R seems to mess up integer addition here.
      state_vector <- state_vector_new
    }
    #print(state_vector[state_vector > 0])
    #Sys.sleep(1)
  }
  
  return(storage)
}


#Example with a single population 
library(xlsx)
library(here)
parameters<-parameters <- read.xlsx2(here("code","parameters.xlsx"),sheetIndex=1, row.names=1)
parameters$Description<-NULL
parameters$Source<-NULL

genotypeNames <- unlist(parameters["genotypeNames",]) #The different genotypes or wolbachia types
populationNames <- parameters["populationNames",1] #fixme: don't hardcode assumption wAlbAB contains all # c("block01") #subpopulations in the metapopulation
numMaleClasses <- as.numeric(unlist(parameters["numMaleClasses",])) # rep(20, length(genotypeNames)) #The number of age classes (days) to use in the model for males of the different genotypes (don't have to be the same for each genotype)
numFemaleClasses <- as.numeric(unlist(parameters["numFemaleClasses",])) # rep(1, length(genotypeNames)) #The number of age classes (days) to use in the model for females of the different genotypes (don't have to be the same for each genotype)
numImmatureClasses <- as.numeric(unlist(parameters["numImmatureClasses",])) #rep(10, length(genotypeNames)) #The number of immature classes (days) to use in the model.  These give rise to a lag between birth and emergence as an adult
names(numMaleClasses) <- genotypeNames
names(numFemaleClasses) <- genotypeNames
names(numImmatureClasses) <- genotypeNames
names(numImmatureClasses) <- genotypeNames
numClassesBySexAndGenotype_matrix <- cbind(numMaleClasses, numFemaleClasses, numImmatureClasses, numImmatureClasses) #Don't edit this line
colnames(numClassesBySexAndGenotype_matrix) <- c("m", "f", "immM", "immF") #don't edit this line
rownames(numClassesBySexAndGenotype_matrix) <- genotypeNames #don't edit this line
carryingCapacityByPopulation_vector <- rep(as.numeric(unlist(parameters["carryingCapacityByPopulation_vector",1])), length(populationNames)) #todo: not hardcode assumptions single population in file #The number of adult mosquitoes in each of the populations at equilibrium
names(carryingCapacityByPopulation_vector) <- populationNames #Don't edit this line
transitionRate_imm <- 1 #don't edit this line: 1 as assuming average time spent in each is 1 day (hence rate is 1/1)
transitionRate_m <- 1 #don't edit this line
transitionRate_f <- 1 #don't edit this line
maleDeathRate <- as.numeric(unlist(parameters["maleDeathRate",])) # 0.2 #the death rate for males (can be different for different genotypes if needed).  Note, this isn't the % mortality per day, it is the rate.
femaleDeathRate <- as.numeric(unlist(parameters["femaleDeathRate",])) # 0.1 #the death rate for females (can be different for different genotypes if needed). Note, this isn't the % mortality per day, it is the rate.
propMated <- as.numeric(unlist(parameters["propMated",])) #0.8 #the proportion of females in the wildtype population that you expect to be mated at equilibrium
femaleBirthRate <- as.numeric(unlist(parameters["femaleBirthRate",])) # 0.4 #the birth rate of females

#Initial population.  Starting out as a vector of zeros (i.e. no mosquitoes)
stateNames <- getStateNames(numMaleClasses, numFemaleClasses, numImmatureClasses, genotypeNames, populationNames, numClassesBySexAndGenotype_matrix ) #don't edit this line
state_vector <- rep(0, length(stateNames)) #don't edit this line
names(state_vector) <- stateNames #don't edit this line
#print(length(state_vector))
#print(length(stateNames))


theta <- maleDeathRate/femaleDeathRate #don't edit this line
F_bar <- carryingCapacityByPopulation_vector*theta*(1 - propMated)/(1 + theta) #don't edit this line
M_bar <- carryingCapacityByPopulation_vector/(1 + theta) #don't edit this line
F_bar_mated <- propMated*carryingCapacityByPopulation_vector*theta/(1 + theta) #don't edit this line
I_bar <- carryingCapacityByPopulation_vector*(femaleDeathRate*theta + maleDeathRate)/(transitionRate_imm*(1 + theta)) #don't edit this line
I_total <- I_bar*numImmatureClasses #don't edit this line #todo: double check don't need to multiply numImmatureClasses by 2
I_max <- I_total/(1 - (transitionRate_imm*I_bar)/(femaleBirthRate*F_bar_mated)) #don't edit this line
femaleMatingRateByPopulation <- transitionRate_imm*I_bar*(1 + theta)/(2*(1 - propMated)*carryingCapacityByPopulation_vector*theta) - femaleDeathRate #don't edit this line
names(femaleMatingRateByPopulation) <- names(carryingCapacityByPopulation_vector) #don't edit this line

ImaxByPopulation_vector <- rep(I_max[1], length(populationNames)) # fixme:this has assumed single genotype per population or something
names(ImaxByPopulation_vector) <- populationNames

#Blank vector to fill in with mating rate parameters
matingRateByName_vector <- rep(NA, length(stateNames)) #don't edit this line
for(p in 1:length(populationNames))
{
  ind <- getIndicesForStateQuery(stateNames, populationName = populationNames[p], genotype = NA, sex = "f", ageClass = NA, mate = NA) #don't edit this line
  matingRateByName_vector[ind] <- femaleMatingRateByPopulation[p] #don't edit this line
}
names(matingRateByName_vector) <- stateNames #don't edit this line


# Populate the state-vector here
print(length(state_vector))
state_vector["block01_wPip_m_1"] <- 0 #make this whatever you want the population of wPip males (first age class) to be in the population named "block01"
state_vector["block01_wPip_f_1_Unmated_None"  ] <- 0 #make this whatever you want the population of unmated wPip females (first age class) to be in the population named "block01"
state_vector["block01_wAlbAB_m_1"  ] <-  200#as.numeric(round(M_bar)) #this is set to be the wildtype population of males in the first age class
state_vector["block01_wAlbAB_f_1_Unmated_None"  ] <- 40#as.numeric(round(F_bar)) #this is set to be the wildtype population of unmated females in the first age class
state_vector["block01_wAlbAB_f_1_wAlbAB_1"  ] <- 160#as.numeric(round(F_bar_mated)) #this is set to be the wildtype population of mated females (by wildtype males) in the first age class
for (imm_ind in seq(numImmatureClasses)) {
  #for (geno_ind in seq(length(genotypeNames))){
    #state_vector[sprintf("%s_%s_immM_%d", populationNames, genotypeNames[geno_ind], imm_ind)] <- round(I_bar[geno_ind]/2)
    #state_vector[sprintf("%s_%s_immF_%d", populationNames, genotypeNames[geno_ind], imm_ind)] <- round(I_bar[geno_ind]/2)
  #}
    state_vector[sprintf("%s_wAlbAB_immM_%d", populationNames, imm_ind)] <- round(I_bar[1]/2)
    state_vector[sprintf("%s_wAlbAB_immF_%d", populationNames, imm_ind)] <- round(I_bar[1]/2)
}
print(length(state_vector))

#Creating the CI parameter vector.  Start by giving every type of individual 0 CI
ciByName_vector <- rep(0, length(stateNames)) #don't change
names(ciByName_vector) <- stateNames #don't change
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wPip", sex = "f", ageClass = NA, mate = "wAlbAB") #find the indices for state names for wPip females that were mated by wAlbAB males
ciByName_vector[ind] <- 1 # give those females a cytoplasmic incompatibility of 1 (no offspring)
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wAlbAB", sex = "f", ageClass = NA, mate = "wPip") #find the indices for state names for wAlbAB females that were mated by wPip males
ciByName_vector[ind] <- 1 # give those females a cytoplasmic incompatibility of 1 (no offspring)
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wPip", sex = "f", ageClass = NA, mate = "wAlbAB", mateStage = 15:19) #find the indices for state names for wPip females that were mated by wAlbAB males who were aged 15 to 19 days old
ciByName_vector[ind] <- 0.68 # give those females a cytoplasmic incompatibility of 0.68 (diminished probability of offspring)
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wPip", sex = "f", ageClass = NA, mate = "wAlbAB", mateStage = 20) #find the indices for state names for wPip females that were mated by wAlbAB males who were aged 20 days old (males can't get any older, from 20 days onwards they are all acting the same way)
ciByName_vector[ind] <- 0 # give those females a cytoplasmic incompatibility of 0 can produce offspring

#create a parameter vector that gives the offspring genotype for each type of mated female
offspringGenotypeByName_vector <- rep("", length(stateNames)) #don't edit
names(offspringGenotypeByName_vector) <- stateNames #don't edit
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wAlbAB", sex = "f", ageClass = NA, mate = "wAlbAB") #find the indices of all wAlbAB females that have been mated by wAlbAB males
offspringGenotypeByName_vector[ind] <- "wAlbAB" #these females will produce wAlbAB offspring
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wAlbAB", sex = "f", ageClass = NA, mate = "wPip") #find the indices of all wAlbAB females that have been mated by wPip males
offspringGenotypeByName_vector[ind] <- "wAlbAB" #these females will produce wAlbAB offspring
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wPip", sex = "f", ageClass = NA, mate = "wPip") #find the indices of all wPip females that have been mated by wPip males
offspringGenotypeByName_vector[ind] <- "wPip" #these females will produce wPip offspring
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wPip", sex = "f", ageClass = NA, mate = "wAlbAB") #find the indices of all wPip females that have been mated by wAlbAB males
offspringGenotypeByName_vector[ind] <- "wPip" #these females will produce wPip offspring

#Define the frieds index for each type of individuals
friedsByName_vector <- rep(0, length(stateNames)) #don't edit
names(friedsByName_vector) <- stateNames #don't edit
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wAlbAB", sex = "m", ageClass = NA, mate = NA) #find the indices of all wAlbAB males
friedsByName_vector[ind] <- 1 #wAlbAB have mating competitiveness 1
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wPip", sex = "m", ageClass = NA, mate = NA) #find the indices of all wPip males
friedsByName_vector[ind] <- 1.0 #wPip have the same competitiveness

#Define the death rates for each type of individuals.  
deathRatesByName_vector <- rep(0, length(stateNames)) #don't edit
names(deathRatesByName_vector) <- stateNames #don't edit
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = NA, sex = "m", ageClass = NA, mate = NA) #don't edit
deathRatesByName_vector[ind] <- maleDeathRate #don't edit
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = NA, sex = "f", ageClass = NA, mate = NA) #don't edit
deathRatesByName_vector[ind] <- femaleDeathRate #don't edit

birthRateByName_vector <- rep(0, length(stateNames)) #don't edit
names(birthRateByName_vector) <- stateNames #don't edit
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = NA, sex = "f", ageClass = NA, mate = genotypeNames) #don't edit
birthRateByName_vector[ind] <- femaleBirthRate #don't edit


transitionRateByName_vector <- rep(0, length(stateNames)) #don't edit
names(transitionRateByName_vector) <- stateNames #don't edit
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = NA, sex = "m", ageClass = NA, mate = NA) #don't edit
transitionRateByName_vector[ind] <- 1 #don't edit
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = NA, sex = "f", ageClass = NA, mate = NA) #don't edit
transitionRateByName_vector[ind] <- 1 #don't edit
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = NA, sex = "immM", ageClass = NA, mate = NA) #don't edit
transitionRateByName_vector[ind] <- 1 #don't edits
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = NA, sex = "immF", ageClass = NA, mate = NA) #don't edit
transitionRateByName_vector[ind] <- 1 #don't edits


#This matrix defines the rate at which each individual in a block moves to other blocks.
# e.g. emmigrationRate_matrix <- matrix(c(0, 1, 1, 0), 2, 2) for when you need to include other blocks 
emmigrationRate_matrix <- matrix(0.0, length(populationNames), length(populationNames)) #This is where you can specify the per capita rates of emmigration from one subpopulation to another.  Applied to males and females in the same way.
colnames(emmigrationRate_matrix) <- populationNames #don't edit
rownames(emmigrationRate_matrix) <- populationNames #don't edit
immigrationRate_vector <- rep(0.0, length(stateNames)) #This is where you can specify the rate at which new individuals enter each subpopulation from outside (i.e. not from one of the other subpopulations)
overallMaleImmigrationRate = 2/7 #two per week
overallFemaleImmigrationRate = 2/7 #two per week
for(i in 1:numMaleClasses["wAlbAB"])
{
  ind_m_wAlb_i <- getIndicesForStateQuery(stateNames, populationName = NA, genotype = "wAlbAB", sex = "m", ageClass = i, mate = NA, mateStage = NA)
  if(i < numMaleClasses["wAlbAB"])
  {
      immigrationRate_vector[ind_m_wAlb_i] <- overallMaleImmigrationRate*(pexp((1/transitionRate_m)*i, rate = maleDeathRate[1]) - pexp((1/transitionRate_m)*(i - 1), rate = maleDeathRate))
  }
  else
  {
    immigrationRate_vector[ind_m_wAlb_i] <- overallMaleImmigrationRate*(1 - pexp((1/transitionRate_m)*i, rate = maleDeathRate[1]))
  }
      
}

for(i in 1:numFemaleClasses["wAlbAB"])
{
  ind_f_wAlb_i <- getIndicesForStateQuery(stateNames, populationName = NA, genotype = "wAlbAB", sex = "f", ageClass = i, mate = c("Unmated", "wAlbAB"), mateStage = NA)
  if(numFemaleClasses["wAlbAB"] == 1)
  {
    immigrationRate_vector[ind_f_wAlb_i] <- overallFemaleImmigrationRate/length(ind_f_wAlb_i)
  }
  else
  {
    if(i < numFemaleClasses["wAlbAB"])
    {
      immigrationRate_vector[ind_f_wAlb_i] <- overallFemaleImmigrationRate*(pexp((1/transitionRate_f)*i, rate = femaleDeathRate[1]) - pexp((1/transitionRate_f)*(i - 1), rate = femaleDeathRate))/length(ind_f_wAlb_i)
    }
    else
    {
      immigrationRate_vector[ind_f_wAlb_i] <- overallFemaleImmigrationRate*(1 - pexp((1/transitionRate_f)*i, rate = femaleDeathRate[1]))/length(ind_f_wAlb_i)
    }
  }
}


names(immigrationRate_vector) <- stateNames #don't edit

#HERE YOU CAN WRITE A FUNCTION THAT RELEASES MOSQUITOES AT SPECIFIED TIMES TO TEST A SCENARIO
#Currently set to not release any mosquitoes
interventionFunction <- function(state_vector, startTime, endTime, interventionTime)
{
  overfloodingDump <- 1
  femaleContamination <- 0.01 # Change this for the contamination rate
  badguyInds <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, genotype = "wAlbAB", sex = "m", ageClass = NA, mate = NA, mateStage = NA)
  goodguyInds <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, genotype = "wPip", sex = "m", ageClass = 1, mate = NA, mateStage = NA)
  #In the following line, you need to decide what state the females are in, mated or unmated?
  contaminatedFemaleInd <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, genotype = "wPip", sex = "f", ageClass = 1, mate = "Unmated", mateStage = NA)
  
  numBadGuys <- sum(state_vector[badguyInds])
  numGoodGuys <- sum(state_vector[goodguyInds])
  numReleased <- overfloodingDump*numBadGuys - numGoodGuys
  numMalesReleased <- rbinom(1, numReleased, (1 - femaleContamination))
  numFemalesReleased <- numReleased - numMalesReleased
  state_vector[goodguyInds] <- state_vector[goodguyInds] + numMalesReleased
  state_vector[contaminatedFemaleInd] <- state_vector[contaminatedFemaleInd] + numFemalesReleased
  
  indPipFemale <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, genotype = "wPip", sex = "f", ageClass = NA, mate = NA, mateStage = NA)
  indAlbABFemale <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, genotype = "wAlbAB", sex = "f", ageClass = NA, mate = NA, mateStage = NA)
  
  propPip <- sum(state_vector[indPipFemale])/(sum(state_vector[indAlbABFemale]) + sum(state_vector[indPipFemale]))
  
  if(interventionTime > 100 & propPip > 0.5)
  {
    state_vector <- state_vector*0  
  }
  
  return(state_vector)
}


#Check that the parameter set meets parametric constraints that are necessary to ensure model behaves properly
runNumber=1
if(all(transitionRate_imm*I_bar/(femaleBirthRate*F_bar_mated) < 1))
{
  modelOutputs <- runModel(state_vector, 0, 150, interventionFunction, interventionTimes = seq(1, 150, 7), recordTimes = 1:150,
                           numClassesBySexAndGenotype_matrix, matingRateByName_vector, ciByName_vector, offspringGenotypeByName_vector, 
                           friedsByName_vector, deathRatesByName_vector, birthRateByName_vector, transitionRateByName_vector, emmigrationRate_matrix, 
                           immigrationRate_vector, carryingCapacityByPopulation_vector, ImaxByPopulation_vector)
  
  
  ts <- aggregateStates(names(state_vector), modelOutputs, populationNames, genotypeNames)
  # save
  write.csv(ts,here("code/outputs/",paste("outputs-run",runNumber,".csv")))
  
  
  for(i in 1:length(ts))
  {
    plot(ts[[i]], xlab = "Days", ylab = "Population", main = names(ts)[[i]])  
  }
  
}else
{
  print("Parameters do not meet constraints")
}









#Fix for broken print
#ID = seq(1,150,1)
#ts$ID <- ID
#plot(ts$ID, ts$block01_wAlbAB_m, xlab = "Days", ylab = "Population", main = "wAlbAB_m")
#plot(ts$ID, ts$block01_wAlbAB_f, xlab = "Days", ylab = "Population", main = "wAlbAB_f")
#plot(ts$ID, ts$block01_wAlbAB_imm, xlab = "Days", ylab = "Population", main = "wAlbAB_Imm")
#plot(ts$ID, ts$block01_wPip_m, xlab = "Days", ylab = "Population", main = "wPip_m")
#plot(ts$ID, ts$block01_wPip_f, xlab = "Days", ylab = "Population", main = "wPip_f")
#plot(ts$ID, ts$block01_wPip_imm, xlab = "Days", ylab = "Population", main = "wPip_imm")

