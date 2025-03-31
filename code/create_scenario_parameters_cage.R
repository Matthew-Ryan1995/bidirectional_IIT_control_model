## Matt Ryan
## 13/12/2023
## Create the parameter csv for a fully factorial scenario exploration
## Requires: pacman, tidyverse


# libraries ---------------------------------------------------------------
pacman::p_load(tidyverse)


# Set parameters ----------------------------------------------------------

initialPopulation <- c(420)
propPip <- c(seq(0.05, 0.5, by = 0.05), 0.22, 0.8)
femaleProportion <- c(0.5)#, 0.67) # Even split and steady state split
direction <- c("bi", "uni", "bi_no_decay")

# create csv --------------------------------------------------------------

df <- expand_grid(
  initialPopulation, propPip, femaleProportion, direction
)

print(str_c("There are ", nrow(df), " scenarios under consideration"))

write_csv(df, "data/scenario_parameters_cage.csv")

