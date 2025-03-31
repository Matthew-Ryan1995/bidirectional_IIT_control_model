## Matt Ryan
## Run results for scenario one
## Required: here, parallel, glue


# load functions ----------------------------------------------------------

source(here::here("code/model_functions.R"))


# load parameters ---------------------------------------------------------

param_list <- readRDS(here::here("data/param_list.Rds"))
list2env(param_list, envir = globalenv())

OF <- 5
FCP <- 1e-2
supression_proportion <- 0.22

stop_intervention <- TRUE

# Intervention function ---------------------------------------------------

interventionFunction <- function(state_vector, startTime, endTime, interventionTime,
                                 overfloodingDump = OF, femaleContamination = FCP, 
                                 f_prop = supression_proportion)
{
  if(stop_intervention){
    return(list(state_vector = state_vector, numReleased = 0))
  }
  
  indPipFemale <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, genotype = "wPip", sex = "f", ageClass = NA, mate = NA, mateStage = NA)
  indAlbABFemale <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, genotype = "wAlbAB", sex = "f", ageClass = NA, mate = NA, mateStage = NA)
  
  propPip <- sum(state_vector[indPipFemale])/(sum(state_vector[indAlbABFemale]) + sum(state_vector[indPipFemale]))
  
  if(isTRUE(interventionTime > 100 & (propPip > f_prop)))
  {
    state_vector <- state_vector  
    numReleased <- 0
    stop_intervention <<- TRUE
  }else{
    # return((list(state_vector = state_vector, numReleased = 0)))
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
    
  }
  
  # Turn off intervention
  # # return((list(state_vector = state_vector, numReleased = 0)))
  # badguyInds <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, genotype = "wAlbAB", sex = "m", ageClass = NA, mate = NA, mateStage = NA)
  # goodguyInds <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, genotype = "wPip", sex = "m", ageClass = 1, mate = NA, mateStage = NA) #query: Why age class 1 here but not above?
  # #In the following line, you need to decide what state the females are in, mated or unmated?
  # contaminatedFemaleInd <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, genotype = "wPip", sex = "f", ageClass = 1, mate = "Unmated", mateStage = NA) #query: Why only wPip?
  # 
  # numBadGuys <- sum(state_vector[badguyInds])
  # numGoodGuys <- sum(state_vector[goodguyInds])
  # # I have no idea where these numbers come from
  # # Only release eggs if there are more bad guys than good guys?
  # # Does this need a stop if negative?
  # numReleased <- overfloodingDump*numBadGuys #??- numGoodGuys
  # # Probabilisticly assign how many males released.
  # numMalesReleased <- rbinom(1, numReleased, (1 - femaleContamination))
  # numFemalesReleased <- numReleased - numMalesReleased
  # state_vector[goodguyInds] <- state_vector[goodguyInds] + numMalesReleased
  # state_vector[contaminatedFemaleInd] <- state_vector[contaminatedFemaleInd] + numFemalesReleased
  # 
  # indPipFemale <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, genotype = "wPip", sex = "f", ageClass = NA, mate = NA, mateStage = NA)
  # indAlbABFemale <- getIndicesForStateQuery(stateNames = names(state_vector), populationName = NA, genotype = "wAlbAB", sex = "f", ageClass = NA, mate = NA, mateStage = NA)
  # 
  # propPip <- sum(state_vector[indPipFemale])/(sum(state_vector[indAlbABFemale]) + sum(state_vector[indPipFemale]))
  # 
  # # Full on stop if more than 100 interventions or more than 50% pip females?
  # if(isTRUE(interventionTime > 100 & propPip > 0.5))
  # {
  #   state_vector <- state_vector*0  
  #   numReleased <- 0
  # }
  
  
  
  return(list(state_vector = state_vector, numReleased = numReleased))
}



# Create initial state vector -----------------------------------------------


parameters <- readxl::read_excel(here::here("code","parameters.xlsx"),sheet=1, 
                                 col_names =TRUE, n_max = 11)
parameters$Description<-NULL
parameters$Source<-NULL

## Code wants dataframe with row names
parameters <- as.data.frame(parameters)
rownames(parameters) <- parameters$Parameter

parameters$Parameter <- NULL


# Wrangling ---------------------------------------------------------------
# 
# genotypeNames <- unlist(parameters["genotypeNames",]) #The different genotypes or wolbachia types
# populationNames <- parameters["populationNames",1] #fixme: don't hardcode assumption wAlbAB contains all # c("block01") #subpopulations in the metapopulation
# numMaleClasses <- as.numeric(unlist(parameters["numMaleClasses",])) # rep(20, length(genotypeNames)) #The number of age classes (days) to use in the model for males of the different genotypes (don't have to be the same for each genotype)
# numFemaleClasses <- as.numeric(unlist(parameters["numFemaleClasses",])) # rep(1, length(genotypeNames)) #The number of age classes (days) to use in the model for females of the different genotypes (don't have to be the same for each genotype)
# numImmatureClasses <- as.numeric(unlist(parameters["numImmatureClasses",])) #rep(10, length(genotypeNames)) #The number of immature classes (days) to use in the model.  These give rise to a lag between birth and emergence as an adult
# names(numMaleClasses) <- genotypeNames
# names(numFemaleClasses) <- genotypeNames
# names(numImmatureClasses) <- genotypeNames #NB: when re-combining the male and female immature classes, I removed what appeared to be a duplicate entry on the now-nonexistent next line and current next line
# numClassesBySexAndGenotype_matrix <- cbind(numMaleClasses, numFemaleClasses, numImmatureClasses) #Don't edit this line
# colnames(numClassesBySexAndGenotype_matrix) <- c("m", "f", "imm") #don't edit this line
# rownames(numClassesBySexAndGenotype_matrix) <- genotypeNames #don't edit this line
# carryingCapacityByPopulation_vector <- rep(as.numeric(unlist(parameters["carryingCapacityByPopulation_vector",1])), length(populationNames)) #todo: not hardcode assumptions single population in file #The number of adult mosquitoes in each of the populations at equilibrium
# names(carryingCapacityByPopulation_vector) <- populationNames #Don't edit this line
# transitionRate_imm <- 1 #don't edit this line: 1 as assuming average time spent in each is 1 day (hence rate is 1/1)
# transitionRate_m <- 1 #don't edit this line
# transitionRate_f <- 1 #don't edit this line
# # maleDeathRate <- as.numeric(unlist(parameters["maleDeathRate",])) # 0.2 #the death rate for males (can be different for different genotypes if needed).  Note, this isn't the % mortality per day, it is the rate.
# maleDeathRate <-  0.2 #the death rate for males (can be different for different genotypes if needed).  Note, this isn't the % mortality per day, it is the rate.
# # femaleDeathRate <- as.numeric(unlist(parameters["femaleDeathRate",])) # 0.1 #the death rate for females (can be different for different genotypes if needed). Note, this isn't the % mortality per day, it is the rate.
# femaleDeathRate <- 0.1 #the death rate for females (can be different for different genotypes if needed). Note, this isn't the % mortality per day, it is the rate.
# # propMated <- as.numeric(unlist(parameters["propMated",])) #0.8 #the proportion of females in the wildtype population that you expect to be mated at equilibrium
# propMated <- 0.8 #the proportion of females in the wildtype population that you expect to be mated at equilibrium
# # femaleBirthRate <- as.numeric(unlist(parameters["femaleBirthRate",])) # 0.4 #the birth rate of females
# femaleBirthRate <- .4#as.numeric(unlist(parameters["femaleBirthRate",])) # 0.4 #the birth rate of females


genotypeNames <- unique(unlist(parameters["genotypeNames",])) #The different unique genotypes of wolbachia types
populationNames <- parameters["populationNames",1] #fixme: don't hardcode assumption wAlbAB contains all # c("block01") #subpopulations in the metapopulation

genotype_selector <- genotypeNames

numMaleClasses <- as.numeric(unlist(parameters["numMaleClasses", genotype_selector])) # rep(20, length(genotypeNames)) #The number of age classes (days) to use in the model for males of the different genotypes (don't have to be the same for each genotype)
numFemaleClasses <- as.numeric(unlist(parameters["numFemaleClasses", genotype_selector])) # rep(1, length(genotypeNames)) #The number of age classes (days) to use in the model for females of the different genotypes (don't have to be the same for each genotype)
numImmatureClasses <- as.numeric(unlist(parameters["numImmatureClasses", genotype_selector])) #rep(10, length(genotypeNames)) #The number of immature classes (days) to use in the model.  These give rise to a lag between birth and emergence as an adult
names(numMaleClasses) <- genotypeNames
names(numFemaleClasses) <- genotypeNames
names(numImmatureClasses) <- genotypeNames #NB: when re-combining the male and female immature classes, I removed what appeared to be a duplicate entry on the now-nonexistent next line and current next line

numClassesBySexAndGenotype_matrix <- cbind(numMaleClasses, numFemaleClasses, numImmatureClasses) #Don't edit this line
colnames(numClassesBySexAndGenotype_matrix) <- c("m", "f", "imm") #don't edit this line
rownames(numClassesBySexAndGenotype_matrix) <- genotypeNames #don't edit this line

carryingCapacityByPopulation_vector <- rep(as.numeric(unlist(parameters["carryingCapacityByPopulation_vector",1])), length(populationNames)) #todo: not hardcode assumptions single population in file #The number of adult mosquitoes in each of the populations at equilibrium
names(carryingCapacityByPopulation_vector) <- populationNames #Don't edit this line
transitionRate_imm <- 1 #don't edit this line: 1 as assuming average time spent in each is 1 day (hence rate is 1/1)
transitionRate_m <- 1 #don't edit this line
transitionRate_f <- 1 #don't edit this line
# maleDeathRate <- as.numeric(unlist(parameters["maleDeathRate",])) # 0.2 #the death rate for males (can be different for different genotypes if needed).  Note, this isn't the % mortality per day, it is the rate.
maleDeathRate <-  0.2 #the death rate for males (can be different for different genotypes if needed).  Note, this isn't the % mortality per day, it is the rate.
# femaleDeathRate <- as.numeric(unlist(parameters["femaleDeathRate",])) # 0.1 #the death rate for females (can be different for different genotypes if needed). Note, this isn't the % mortality per day, it is the rate.
femaleDeathRate <- 0.1 #the death rate for females (can be different for different genotypes if needed). Note, this isn't the % mortality per day, it is the rate.
# propMated <- as.numeric(unlist(parameters["propMated",])) #0.8 #the proportion of females in the wildtype population that you expect to be mated at equilibrium
propMated <- 0.8 #the proportion of females in the wildtype population that you expect to be mated at equilibrium
# femaleBirthRate <- as.numeric(unlist(parameters["femaleBirthRate",])) # 0.4 #the birth rate of females
femaleBirthRate <- .4#as.numeric(unlist(parameters["femaleBirthRate",])) # 0.4 #the birth rate of females

propFemale <- as.numeric(unlist(parameters["propFemale", genotype_selector])) # 0.5 # the proportion of immature mosquitoes that mature into adult females
# names(maleDeathRate) <- genotypeNames # Name columns
# names(femaleDeathRate) <- genotypeNames
# names(propMated) <- genotypeNames
# names(femaleBirthRate) <- genotypeNames
# names(propFemale) <- genotypeNames



# carryingCapacityByPopulation_vector[1] <- 100

#Initial population.  Starting out as a vector of zeros (i.e. no mosquitoes)
stateNames <- getStateNames(numMaleClasses, numFemaleClasses, numImmatureClasses, genotypeNames, populationNames, numClassesBySexAndGenotype_matrix ) #don't edit this line
state_vector <- rep(0, length(stateNames)) #don't edit this line
names(state_vector) <- stateNames #don't edit this line

theta <- maleDeathRate/femaleDeathRate #don't edit this line
phi <- 1/(1+maleDeathRate)
F_bar <- carryingCapacityByPopulation_vector*theta*(1 - propMated)/(1 + theta) #don't edit this line
M_bar <- carryingCapacityByPopulation_vector/((1 + theta)) #don't edit this line
F_bar_mated <- propMated*carryingCapacityByPopulation_vector*theta/(1 + theta) #don't edit this line
M1 <- M_bar / ((1-phi^(numMaleClasses - 1))/(1-phi) + (phi^(numMaleClasses - 2))/maleDeathRate)

I_bar <- 2 * M1[1]/(phi * transitionRate_imm)
I_total <- numImmatureClasses*I_bar
I_max <- I_total/(1 - (transitionRate_imm*I_bar)/(femaleBirthRate*F_bar_mated))

femaleMatingRateByPopulation <- ((transitionRate_imm/2) * I_bar - femaleDeathRate * F_bar)/(F_bar * M_bar)
# femaleMatingRateByPopulation <- 0.4
names(femaleMatingRateByPopulation) <- names(carryingCapacityByPopulation_vector) #don't edit this line


matingRateByName_vector <- rep(NA, length(stateNames)) #don't edit this line
for(p in 1:length(populationNames))
{
  ind <- getIndicesForStateQuery(stateNames, populationName = populationNames[p], genotype = NA, sex = "f", ageClass = NA, mate = "Unmated") #don't edit this line
  matingRateByName_vector[ind] <- femaleMatingRateByPopulation[p] * M_bar[p] #don't edit this line
}
names(matingRateByName_vector) <- stateNames #don't edit this line

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


populate_m_values <- function(k, numMaleClasses, phi, maleDeathRate, M1, M_bar){
  if(numMaleClasses[1]==1){
    return(round(M_bar, 1))
  }
  if(k <= numMaleClasses[1]){
    return(round(phi[1]^(k-1) * M1[1]))
  }else{
    return(round(phi[1]^(k-2)/maleDeathRate[1] * M1[1]))
  }
}
populate_f_m_values <- function(k, femaleMatingRateByPopulation, femaleDeathRate, state_vector, F_bar){
  return(round(femaleMatingRateByPopulation[1] * state_vector[sprintf("%s_wAlbAB_m_%d", populationNames, k)] * F_bar[1] / femaleDeathRate[1] ))
}
# Populate the state-vector here
print(length(state_vector))
state_vector["block01_wPip_m_1"] <- 0 #make this whatever you want the population of wPip males (first age class) to be in the population named "block01"
state_vector["block01_wPip_f_1_Unmated_None"  ] <- 0 #make this whatever you want the population of unmated wPip females (first age class) to be in the population named "block01"
# state_vector["block01_wAlbAB_m_1"  ] <-  as.numeric(round(M_bar)) #this is set to be the wildtype population of males in the first age class
state_vector["block01_wAlbAB_f_1_Unmated_None"  ] <- as.numeric(round(F_bar[1])) #this is set to be the wildtype population of unmated females in the first age class
# state_vector["block01_wAlbAB_f_1_wAlbAB_1"  ] <- as.numeric(round(F_bar_mated[1])) #this is set to be the wildtype population of mated females (by wildtype males) in the first age class
for (k in seq(numMaleClasses[1])){
  state_vector[sprintf("%s_wAlbAB_m_%d", populationNames, k)] <- populate_m_values(k, numMaleClasses, phi, maleDeathRate, M1, M_bar)
  state_vector[sprintf("%s_wAlbAB_f_1_wAlbAB_%d", populationNames, k)] <- populate_f_m_values(k, femaleMatingRateByPopulation, femaleDeathRate, state_vector, F_bar)
}
for (imm_ind in seq(numImmatureClasses[1])) {
  state_vector[sprintf("%s_wAlbAB_imm_%d", populationNames, imm_ind)] <- round(I_bar[1])
}

ImaxByPopulation_vector <- rep(I_max[1], length(populationNames)) # fixme:this has assumed single genotype per population or something
names(ImaxByPopulation_vector) <- populationNames



# immigrationRate_vector <- emigrationRate_vector <- rep(0, length(state_vector))
state_vector <- create_cage_state_vector(state_vector = state_vector, propPip = 0.,
                                         total_capacity = carryingCapacityByPopulation_vector[1], female_split = 0.5)
# tmp <- names(state_vector)
# state_vector <- rep(1, length(state_vector))
# names(state_vector) <- tmp

# Run model ---------------------------------------------------------------

# recordTimes <- c(5)
# recordTimeIndex <- 1

#Check that the parameter set meets parametric constraints that are necessary 
#to ensure model behaves properly
set.seed(1234)
runNumber=1
num_runs <- 1#50
start_time <- 0
end_time <- 500
outputs <- list()
tictoc::tic()
for(j in 1:num_runs){
  if(all(transitionRate_imm*I_bar/(femaleBirthRate*F_bar_mated) < 1))
  {
    
    modelOutputs <- runModel(state_vector, start_time, end_time, interventionFunction, 
                             interventionTimes = seq(start_time+1, end_time, 7), 
                             recordTimes = (start_time+1):end_time, numClassesBySexAndGenotype_matrix, 
                             matingRateByName_vector, ciByName_vector, offspringGenotypeByName_vector, 
                             friedsByName_vector, deathRatesByName_vector, birthRateByName_vector, 
                             transitionRateByName_vector, emmigrationRate_matrix, 
                             immigrationRate_vector, emigrationRate_vector, carryingCapacityByPopulation_vector, 
                             ImaxByPopulation_vector, progress=TRUE, cage=FALSE)
    
    ts <- aggregateStates(names(state_vector), modelOutputs, populationNames, genotypeNames)
    # save
    # write.csv(ts,here("code/outputs/",paste("outputs-run",runNumber,".csv")))
    
    # for(i in 1:length(ts))
    # {
    #   plot(ts[[i]], xlab = "Days", ylab = "Population", main = names(ts)[[i]], type="l")
    #   # lines(ts[[i]])
    # }
    
  }else
  {
    print("Parameters do not meet constraints")
  }
  # outputs[[j]] <- ts[[1]]
}

plot(ts[[1]], xlab = "Days", ylab = "Population", main = names(ts)[[1]], type="l")
abline(M_bar, 0)
plot(ts[[2]], xlab = "Days", ylab = "Population", main = names(ts)[[2]], type="l")
abline(F_bar + F_bar_mated, 0)
plot(ts[[3]], xlab = "Days", ylab = "Population", main = names(ts)[[3]], type="l")
abline(I_total, 0)

tictoc::toc()