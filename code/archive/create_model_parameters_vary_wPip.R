# Matt Ryan
# Create the parameters needed for the simulations
# These parameters sets will keep wAlb at the expected level, but wPip at the High or Low ranges

# libraries ---------------------------------------------------------------

library(readxl)
library(readr)
library(here)
library(progress)

# todo: Large wPip with small wAlb and vice versa

# functions ---------------------------------------------------------------

source(here("code/model_functions.R"))


# Parameter flags ---------------------------------------------------------
## Allows you to choose either expected value (both false) or min/max values of parameters
MIN_VALUES <- TRUE
# MAX_VALUES <- FALSE

wpip_level <- ifelse(MIN_VALUES, "low", "high")

## Allows you to choose whether to transform the birth rate or the death rate to allow for 
## biologically feasible conditions
TRANSFORM_DEATHRATE <- TRUE
TRANSFORM_BIRTHRATE <- FALSE  # todo: make sure that if TRANSOFRM_DEATHRATE is True, this is set to false

param_save_name <- glue::glue("data/param_list_wPip_{wpip_level}_TRANSFORM_DEATHRATE_{TRANSFORM_DEATHRATE}_TRANSFORM_BIRTHRATE_{TRANSFORM_BIRTHRATE}.Rds")

# if(isTRUE(MIN_VALUES & MAX_VALUES)){
#   stop("Choose one, minimum or maximum parameter values")
# }

if(isTRUE(TRANSFORM_DEATHRATE & TRANSFORM_BIRTHRATE)){
  stop("Choose one, transform birth or death rate")
}

# Load data ---------------------------------------------------------------

parameters <- read_excel(here("code","parameters.xlsx"),
                         sheet=1, col_names =TRUE, n_max = 11)
parameters$Description<-NULL
parameters$Source<-NULL

## Code wants dataframe with row names
parameters <- as.data.frame(parameters)
rownames(parameters) <- parameters$Parameter

parameters$Parameter <- NULL


# Wrangling ---------------------------------------------------------------

genotypeNames <- unique(unlist(parameters["genotypeNames",])) #The different unique genotypes of wolbachia types
populationNames <- parameters["populationNames",1] #fixme: don't hardcode assumption wAlbAB contains all # c("block01") #subpopulations in the metapopulation

genotype_selector_0 <- genotypeNames
if(isTRUE(MIN_VALUES)){
  genotype_selector <- genotype_selector_0
  genotype_selector[2] <- paste(genotype_selector_0[2], "_min", sep = "")
  death_selector <- genotype_selector_0
  death_selector[2] <- paste(genotype_selector_0[2], "_max", sep = "") # The larger death rate (smaller life span) goes with lower parameters and long birth rate: This makes small population
}else{
  genotype_selector <- genotype_selector_0
  genotype_selector[2] <- paste(genotype_selector_0[2], "_max", sep = "")
  death_selector <- genotype_selector_0
  death_selector[2] <- paste(genotype_selector_0[2], "_min", sep = "") # The smaller death rate (long life) goes with larger parameters: This makes large population
}

# Keep population parameters constant
# Change only birth and death rates
numMaleClasses <- as.numeric(unlist(parameters["numMaleClasses", genotype_selector_0])) # rep(20, length(genotypeNames)) #The number of age classes (days) to use in the model for males of the different genotypes (don't have to be the same for each genotype)
numFemaleClasses <- as.numeric(unlist(parameters["numFemaleClasses", genotype_selector_0])) # rep(1, length(genotypeNames)) #The number of age classes (days) to use in the model for females of the different genotypes (don't have to be the same for each genotype)
numImmatureClasses <- as.numeric(unlist(parameters["numImmatureClasses", genotype_selector_0])) #rep(10, length(genotypeNames)) #The number of immature classes (days) to use in the model.  These give rise to a lag between birth and emergence as an adult
names(numMaleClasses) <- genotypeNames
names(numFemaleClasses) <- genotypeNames
names(numImmatureClasses) <- genotypeNames #NB: when re-combining the male and female immature classes, I removed what appeared to be a duplicate entry on the now-nonexistent next line and current next line

numClassesBySexAndGenotype_matrix <- cbind(numMaleClasses, numFemaleClasses, numImmatureClasses) #Don't edit this line
colnames(numClassesBySexAndGenotype_matrix) <- c("m", "f", "imm") #don't edit this line
rownames(numClassesBySexAndGenotype_matrix) <- genotypeNames #don't edit this line

carryingCapacityByPopulation_vector <- rep(as.numeric(unlist(parameters["carryingCapacityByPopulation_vector", genotype_selector[1]])), length(populationNames)) #todo: not hardcode assumptions single population in file #The number of adult mosquitoes in each of the populations at equilibrium
names(carryingCapacityByPopulation_vector) <- populationNames #Don't edit this line

transitionRate_imm <- 1 #don't edit this line: 1 as assuming average time spent in each is 1 day (hence rate is 1/1)
transitionRate_m <- 1 #don't edit this line
transitionRate_f <- 1 #don't edit this line

maleDeathRate <- as.numeric(unlist(parameters["maleDeathRate", death_selector])) # 0.2 #the death rate for males (can be different for different genotypes if needed).  Note, this isn't the % mortality per day, it is the rate.
femaleDeathRate <- as.numeric(unlist(parameters["femaleDeathRate", death_selector])) # 0.1 #the death rate for females (can be different for different genotypes if needed). Note, this isn't the % mortality per day, it is the rate.
propMated <- as.numeric(unlist(parameters["propMated", genotype_selector])) #0.8 #the proportion of females in the wildtype population that you expect to be mated at equilibrium
femaleBirthRate <- as.numeric(unlist(parameters["femaleBirthRate", death_selector])) # 0.4 #the birth rate of females
propFemale <- as.numeric(unlist(parameters["propFemale", genotype_selector])) # 0.5 # the proportion of immature mosquitoes that mature into adult females
names(maleDeathRate) <- genotypeNames # Name columns
names(femaleDeathRate) <- genotypeNames
names(propMated) <- genotypeNames
names(femaleBirthRate) <- genotypeNames
names(propFemale) <- genotypeNames

if(isTRUE(TRANSFORM_DEATHRATE & any(femaleDeathRate/(propMated*propFemale*femaleBirthRate) >= 1))){
  femaleDeathRate <- round(propMated*propFemale*femaleBirthRate, 3) - 0.001 # biological feasibility condition + fudge factor to ensure inequality
}
if(isTRUE(TRANSFORM_BIRTHRATE & any(femaleDeathRate/(propMated*propFemale*femaleBirthRate) >= 1))){
  femaleBirthRate <- round(femaleDeathRate/(propMated*propFemale), 3) + 0.001 # biological feasibility condition - fudge factor to ensure inequality
}

propMale <- 1 - propFemale  # so don't have 1 - propFemale everywhere

# Check biologically feasible constraint that
## \mu_F / (p_F p_mated \lambda) < 1
if(any(femaleDeathRate/(propMated*propFemale*femaleBirthRate) >= 1)){
  stop("Parameters do not meet biologically feasible conditions")
}

# ::NEED CHANGING::
# carryingCapacityByPopulation_vector[1] <- 100

# Create initial population and parameters -----------------------------------------------

#Initial population.  Starting out as a vector of zeros (i.e. no mosquitoes)
stateNames <- getStateNames(numMaleClasses, numFemaleClasses, numImmatureClasses, genotypeNames, populationNames, numClassesBySexAndGenotype_matrix ) #don't edit this line
state_vector <- rep(0, length(stateNames)) #don't edit this line
names(state_vector) <- stateNames #don't edit this line

theta <- maleDeathRate*propFemale/(femaleDeathRate*propMale) #don't edit this line
phi <- transitionRate_m/(transitionRate_m+maleDeathRate)
F_bar <- carryingCapacityByPopulation_vector*theta*(1 - propMated)/(1 + theta) #don't edit this line
M_bar <- carryingCapacityByPopulation_vector/(1 + theta) #don't edit this line
F_bar_mated <- propMated*carryingCapacityByPopulation_vector*theta/(1 + theta) #don't edit this line
M1 <- M_bar / ((1-phi^(numMaleClasses - 1))/(1-phi) + (transitionRate_m*phi^(numMaleClasses - 2))/maleDeathRate)

I_bar <- transitionRate_m* M1/(phi * transitionRate_imm * propMale)
I_total <- numImmatureClasses*I_bar
I_max <- I_total/(1 - (transitionRate_imm*I_bar)/(femaleBirthRate*F_bar_mated))

femaleMatingRateByPopulation <- ((transitionRate_imm * propFemale) * I_bar - femaleDeathRate * F_bar)/(F_bar * M_bar)
# femaleMatingRateByPopulation <- 0.4
names(femaleMatingRateByPopulation) <- names(carryingCapacityByPopulation_vector) #don't edit this line 

matingRateByName_vector <- rep(NA, length(stateNames)) #don't edit this line
for(p in 1:length(genotypeNames))
{
  ind <- getIndicesForStateQuery(stateNames, populationName = populationNames[1], genotype = genotypeNames[p], sex = "f", ageClass = NA, mate = "Unmated") #don't edit this line
  matingRateByName_vector[ind] <- femaleMatingRateByPopulation[p] * M_bar[p] #don't edit this line
}
# for(p in 1:length(populationNames))
# {
#   ind <- getIndicesForStateQuery(stateNames, populationName = populationNames[p], genotype = NA, sex = "f", ageClass = NA, mate = "Unmated") #don't edit this line
#   matingRateByName_vector[ind] <- femaleMatingRateByPopulation[p] * M_bar[p] #don't edit this line
# }
names(matingRateByName_vector) <- stateNames #don't edit this line

populate_m_values <- function(k, numMaleClasses, phi, maleDeathRate, M1, M_bar){
  # This function creates the steady state vector values for the male classes,
  # which depends on the wild type population.
  
  # Hardcoded that wild type is genotype 1 (wAlbAB)
  if(numMaleClasses[1]==1){ 
    return(ceiling(M_bar, 1))
  }
  if(k <= numMaleClasses[1]){
    return(ceiling(phi[1]^(k-1) * M1[1]))
  }else{
    return(ceiling(transitionRate_M*phi[1]^(k-2)/maleDeathRate[1] * M1[1]))
  }
}
populate_f_m_values <- function(k, femaleMatingRateByPopulation, femaleDeathRate, state_vector, F_bar){
  return(ceiling(femaleMatingRateByPopulation[1] * state_vector[sprintf("%s_wAlbAB_m_%d", populationNames, k)] * F_bar[1] / femaleDeathRate[1] ))
}
# Populate the state-vector here
print(length(state_vector))
state_vector["block01_wPip_m_1"] <- 0 #make this whatever you want the population of wPip males (first age class) to be in the population named "block01"
state_vector["block01_wPip_f_1_Unmated_None"  ] <- 0 #make this whatever you want the population of unmated wPip females (first age class) to be in the population named "block01"
# state_vector["block01_wAlbAB_m_1"  ] <-  as.numeric(round(M_bar)) #this is set to be the wildtype population of males in the first age class
state_vector["block01_wAlbAB_f_1_Unmated_None"  ] <- as.numeric(ceiling(F_bar[1])) #this is set to be the wildtype population of unmated females in the first age class
# state_vector["block01_wAlbAB_f_1_wAlbAB_1"  ] <- as.numeric(round(F_bar_mated[1])) #this is set to be the wildtype population of mated females (by wildtype males) in the first age class
for (k in seq(numMaleClasses[1])){
  state_vector[sprintf("%s_wAlbAB_m_%d", populationNames, k)] <- populate_m_values(k, numMaleClasses, phi, maleDeathRate, M1, M_bar)
  state_vector[sprintf("%s_wAlbAB_f_1_wAlbAB_%d", populationNames, k)] <- populate_f_m_values(k, femaleMatingRateByPopulation, femaleDeathRate, state_vector, F_bar)
}
for (imm_ind in seq(numImmatureClasses[1])) {
  state_vector[sprintf("%s_wAlbAB_imm_%d", populationNames, imm_ind)] <- ceiling(I_bar[1])  # todo: explore making this ceiling not round?
}
print(length(state_vector))

# state_vector[state_vector>0] <- state_vector[state_vector>0] + 10

ImaxByPopulation_vector <- rep(I_max[1], length(populationNames)) # fixme:this has assumed single genotype per population or something
names(ImaxByPopulation_vector) <- populationNames


#Creating the CI parameter vector.  Start by giving every type of individual 0 CI
## If a wAlbAB female is mated with a wPip, they cannot give birth (CI=1)
## If a wPip female is mated with a wAlbAB male age 1:15, the cannot give birth (CI=1)
## If a wPip female is mated with a wAlbAB male age 16:19, the can give birth at a reduced rate (CI=0.68)
## If a wPip female is mated with a wAlbAB male age 20, the can give birth (CI=0)
## This means that, although the CI is bi-directional, it is in favour of the wPip
## Taking into account that we assume wolbachia strain is passed on from the females (wPip mother means wPip baby),
## this means there is a slight reproductive advantage (probably wrong terminology) to the wPip.
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
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wAlbAB", sex = "m", ageClass = NA, mate = NA) #don't edit
deathRatesByName_vector[ind] <- maleDeathRate[1] #don't edit
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wPip", sex = "m", ageClass = NA, mate = NA) #don't edit
deathRatesByName_vector[ind] <- maleDeathRate[2] #don't edit
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wAlbAB", sex = "f", ageClass = NA, mate = NA) #don't edit
deathRatesByName_vector[ind] <- femaleDeathRate[1] #don't edit
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wPip", sex = "f", ageClass = NA, mate = NA) #don't edit
deathRatesByName_vector[ind] <- femaleDeathRate[2] #don't edit

birthRateByName_vector <- rep(0, length(stateNames)) #don't edit
names(birthRateByName_vector) <- stateNames #don't edit
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wAlbAB", sex = "f", ageClass = NA, mate = genotypeNames) #don't edit
birthRateByName_vector[ind] <- femaleBirthRate[1] #don't edit
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = "wPip", sex = "f", ageClass = NA, mate = genotypeNames) #don't edit
birthRateByName_vector[ind] <- femaleBirthRate[2] #don't edit

# This variable is defined as 
# a named vector that gives the per capita rate of transitions through stages for all male, female and immature state names
# Why is it a vector of ones?
transitionRateByName_vector <- rep(0, length(stateNames)) #don't edit
names(transitionRateByName_vector) <- stateNames #don't edit
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = NA, sex = "m", ageClass = NA, mate = NA) #don't edit
transitionRateByName_vector[ind] <- 1 #don't edit
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = NA, sex = "f", ageClass = NA, mate = NA) #don't edit
transitionRateByName_vector[ind] <- 1 #don't edit
ind <- getIndicesForStateQuery(stateNames, populationName = populationNames, genotype = NA, sex = "imm", ageClass = NA, mate = NA) #don't edit
transitionRateByName_vector[ind] <- 1 #don't edits


#This matrix defines the rate at which each individual in a block moves to other blocks.
# e.g. emmigrationRate_matrix <- matrix(c(0, 1, 1, 0), 2, 2) for when you need to include other blocks 
emmigrationRate_matrix <- matrix(0.0, length(populationNames), length(populationNames)) #This is where you can specify the per capita rates of emmigration from one subpopulation to another.  Applied to males and females in the same way.
colnames(emmigrationRate_matrix) <- populationNames #don't edit
rownames(emmigrationRate_matrix) <- populationNames #don't edit
immigrationRate_vector <- rep(0.0, length(stateNames)) #This is where you can specify the rate at which new individuals enter each subpopulation from outside (i.e. not from one of the other subpopulations)
overallMaleImmigrationRate = 2/7 #two per week
overallFemaleImmigrationRate = 2/7 #two per week

immigration_and_emigration_rates <- create_immigration_emmigration_vectors(state_vector = state_vector, overallMaleImmigrationRate = overallMaleImmigrationRate, 
                                                                           overallFemaleImmigrationRate = overallFemaleImmigrationRate, maleDeathRate = maleDeathRate, 
                                                                           numMaleClasses = numMaleClasses, phi = phi, theta = theta, M1 = M1, carryingCapacityByPopulation_vector=carryingCapacityByPopulation_vector)
immigrationRate_vector <- immigration_and_emigration_rates$immigrationRate_vector
emigrationRate_vector <- immigration_and_emigration_rates$emigrationRate_vector
# for(i in 1:numMaleClasses["wAlbAB"])
# {
#   ind_m_wAlb_i <- getIndicesForStateQuery(stateNames, populationName = NA, genotype = "wAlbAB", sex = "m", ageClass = i, mate = NA, mateStage = NA)
#   if(i < numMaleClasses["wAlbAB"])
#   {
#     immigrationRate_vector[ind_m_wAlb_i] <- overallMaleImmigrationRate*(pexp((1/transitionRate_m)*i, rate = maleDeathRate[1]) - pexp((1/transitionRate_m)*(i - 1), rate = maleDeathRate[1]))
#   }
#   else
#   {
#     immigrationRate_vector[ind_m_wAlb_i] <- overallMaleImmigrationRate*(1 - pexp((1/transitionRate_m)*i, rate = maleDeathRate[1]))
#   }
#   
# }
# 
# for(i in 1:numFemaleClasses["wAlbAB"])
# {
#   ind_f_wAlb_i <- getIndicesForStateQuery(stateNames, populationName = NA, genotype = "wAlbAB", sex = "f", ageClass = i, mate = c("Unmated", "wAlbAB"), mateStage = NA)
#   if(numFemaleClasses["wAlbAB"] == 1)
#   {
#     immigrationRate_vector[ind_f_wAlb_i] <- overallFemaleImmigrationRate/length(ind_f_wAlb_i)
#   }
#   else
#   {
#     if(i < numFemaleClasses["wAlbAB"])
#     {
#       immigrationRate_vector[ind_f_wAlb_i] <- overallFemaleImmigrationRate*(pexp((1/transitionRate_f)*i, rate = femaleDeathRate[1]) - pexp((1/transitionRate_f)*(i - 1), rate = femaleDeathRate[1]))/length(ind_f_wAlb_i)
#     }
#     else
#     {
#       immigrationRate_vector[ind_f_wAlb_i] <- overallFemaleImmigrationRate*(1 - pexp((1/transitionRate_f)*i, rate = femaleDeathRate[1]))/length(ind_f_wAlb_i)
#     }
#   }
# }


# Save parameters ---------------------------------------------------------

parameter_set <- glue::glue(
  "wPip_{wpip_level}_numMaleClasses_{numMaleClasses[1]}_numFemaleClasses_{numFemaleClasses[1]}_numImmatureClasses_{numImmatureClasses[1]}_carryingCapacityByPopulation_vector_{carryingCapacityByPopulation_vector[1]}_maleDeathRate_{maleDeathRate[2]}_femaleDeathRate_{femaleDeathRate[2]}_propMated_{propMated[2]}_femaleBirthRate_{femaleBirthRate[2]}"
)

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
  emigrationRate_vector = emigrationRate_vector,
  carryingCapacityByPopulation_vector = carryingCapacityByPopulation_vector,
  ImaxByPopulation_vector = ImaxByPopulation_vector,
  populationNames = populationNames,
  genotypeNames = genotypeNames,
  maleDeathRate = maleDeathRate, 
  numMaleClasses = numMaleClasses,
  phi = phi,
  theta = theta,
  M1 = M1,
  parameter_set = parameter_set
)


write_rds(param_list, param_save_name)