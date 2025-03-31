# Matt Ryan
# Create the parameters needed for the simulations

# libraries ---------------------------------------------------------------

library(readxl)
library(here)
library(progress)


# functions ---------------------------------------------------------------

source(here("code/model_functions.R"))


# Load data ---------------------------------------------------------------

parameters <- read_excel(here("code","parameters.xlsx"),sheet=1, col_names =TRUE)
parameters$Description<-NULL
parameters$Source<-NULL

## Code wants dataframe with row names
parameters <- as.data.frame(parameters)
rownames(parameters) <- parameters$Parameter

parameters$Parameter <- NULL


# Wrangling ---------------------------------------------------------------

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



# ::NEED CHANGING::
# carryingCapacityByPopulation_vector[1] <- 100

# Create initial population and parameters -----------------------------------------------

#Initial population.  Starting out as a vector of zeros (i.e. no mosquitoes)
stateNames <- getStateNames(numMaleClasses, numFemaleClasses, numImmatureClasses, genotypeNames, populationNames, numClassesBySexAndGenotype_matrix ) #don't edit this line
state_vector <- rep(0, length(stateNames)) #don't edit this line
names(state_vector) <- stateNames #don't edit this line

# I think most of these all still hold when thinking about total values, not
# Particularly I think the I values are talking about "total immatures (M + F)"
theta <- maleDeathRate/femaleDeathRate #don't edit this line
phi <- 1/(1+maleDeathRate)
F_bar <- carryingCapacityByPopulation_vector*theta*(1 - propMated)/(1 + theta) #don't edit this line
M_bar <- carryingCapacityByPopulation_vector/((1 + theta)) #don't edit this line
F_bar_mated <- propMated*carryingCapacityByPopulation_vector*theta/(1 + theta) #don't edit this line
M1 <- M_bar / ((1-phi^(numMaleClasses - 1))/(1-phi) + (phi^(numMaleClasses - 2))/maleDeathRate)

I_bar <- M1/(phi * transitionRate_imm)
I_total <- 2*numImmatureClasses*I_bar
I_max <- I_total/(1 - (transitionRate_imm*I_bar * 2)/(femaleBirthRate*F_bar_mated))

femaleMatingRateByPopulation <- (transitionRate_imm * I_bar - femaleDeathRate * F_bar)/(F_bar * M_bar)
# femaleMatingRateByPopulation <- 0.4
names(femaleMatingRateByPopulation) <- names(carryingCapacityByPopulation_vector) #don't edit this line

matingRateByName_vector <- rep(NA, length(stateNames)) #don't edit this line
for(p in 1:length(populationNames))
{
  ind <- getIndicesForStateQuery(stateNames, populationName = populationNames[p], genotype = NA, sex = "f", ageClass = NA, mate = "Unmated") #don't edit this line
  matingRateByName_vector[ind] <- femaleMatingRateByPopulation[p] * M_bar[p] #don't edit this line
}
names(matingRateByName_vector) <- stateNames #don't edit this line

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
  state_vector[sprintf("%s_wAlbAB_immM_%d", populationNames, imm_ind)] <- round(I_bar[1])
  state_vector[sprintf("%s_wAlbAB_immF_%d", populationNames, imm_ind)] <- round(I_bar[1])
}
print(length(state_vector))

# state_vector[state_vector>0] <- state_vector[state_vector>0] + 10

ImaxByPopulation_vector <- rep(I_max[1], length(populationNames)) # fixme:this has assumed single genotype per population or something
names(ImaxByPopulation_vector) <- populationNames


#Creating the CI parameter vector.  Start by giving every type of individual 0 CI
# Assumes pip can bread with Alb, but not vice versa
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

# This variable is defined as 
# a named vector that gives the per capita rate of transitions through stages for all male, female and immature state names
# Why is it a vector of ones?
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
    immigrationRate_vector[ind_m_wAlb_i] <- overallMaleImmigrationRate*(pexp((1/transitionRate_m)*i, rate = maleDeathRate[1]) - pexp((1/transitionRate_m)*(i - 1), rate = maleDeathRate[1]))
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
      immigrationRate_vector[ind_f_wAlb_i] <- overallFemaleImmigrationRate*(pexp((1/transitionRate_f)*i, rate = femaleDeathRate[1]) - pexp((1/transitionRate_f)*(i - 1), rate = femaleDeathRate[1]))/length(ind_f_wAlb_i)
    }
    else
    {
      immigrationRate_vector[ind_f_wAlb_i] <- overallFemaleImmigrationRate*(1 - pexp((1/transitionRate_f)*i, rate = femaleDeathRate[1]))/length(ind_f_wAlb_i)
    }
  }
}


# Save parameters ---------------------------------------------------------


param_list <- list(
  transitionRate_imm = transitionRate_imm,
  I_bar = I_bar,
  femaleBirthRate = femaleBirthRate,
  F_bar_mated = F_bar_mated,
  state_vector = state_vector,
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
  carryingCapacityByPopulation_vector = carryingCapacityByPopulation_vector,
  ImaxByPopulation_vector = ImaxByPopulation_vector,
  populationNames = populationNames,
  genotypeNames = genotypeNames
)

write_rds(param_list, "data/param_list.Rds")