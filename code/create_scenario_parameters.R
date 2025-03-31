## Matt Ryan
## 13/12/2023
## Create the parameter csv for a fully factorial scenario exploration
## Requires: pacman, tidyverse


# libraries ---------------------------------------------------------------
pacman::p_load(tidyverse)


# Set parameters ----------------------------------------------------------

FCP <- c(0.01)
OF <- c(5)
female_immigration_rate <- male_immigration_rate <- c(0, 2, 10)/7
control_proportion <-  c(0.4) # Chosen from cage experiments
intervention_type <- c("complete_stop", "maintain", "elimination")

# Complete stop - When wPip reach certain level, stop releasing wPip
## Aim to suppress wAlb (Completely stop intervention)
## Stop wPip release when reach threshold A
# Maintain - Try to keep wPip at certain level, aim to suppress wAlb
## Aim to suppress wAlb (Completely stop intervention)
## Stop wPip release when reach threshold A
## Start wPip release if drop below threshold B
# Never stop until suppression - Aim to suppress wAlb
## Aim to suppress wAlb (Completely stop intervention)
## What proportion of sims at stopping point do we consider wPip established at suppression
## What proportion are reversible after 90 day period due to immigration/emigration

# Want to look at: 
## Probability of wPip establishment or reversibility (given immigration/emigration)
### After 3 month (90 day) period
## Probability of wAlb suppression (10% of steady state value)
### After 90 day period
## Cost of intervention (number released)


# create csv --------------------------------------------------------------

# Considering three imigration types:
## No immigration
## Low immigration (2 per week for F and M)
## High immigration (10 per week for F and M)
df <- expand_grid(
  FCP, OF, female_immigration_rate, male_immigration_rate, control_proportion, intervention_type
) %>% 
  filter(female_immigration_rate==male_immigration_rate) 

print(str_c("There are ", nrow(df), " scenarios under consideration"))

write_csv(df, "data/scenario_parameters.csv")

