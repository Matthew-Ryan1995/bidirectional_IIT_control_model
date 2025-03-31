## Matt Ryan
## 12-11-2023
## Collect data for plotting
## Requires: tidyverse, glue, here


# libraries ---------------------------------------------------------------

pacman::p_load(tidyverse)


# Functions ---------------------------------------------------------------

# Need indices query function
source(here::here("code/model_functions.R"))

get_data <- function(f){
  column_names <- c("days", "block01_wAlbAB_m", "block01_wAlbAB_f", "block01_wAlbAB_immM",
                    "block01_wAlbAB_immF", "block01_wPip_m", "block01_wPip_f", "block01_wPip_immM",
                    "block01_wPip_immF", "total_released_to_date")
  df <- read_csv(f,
                 col_types = cols(), col_names = column_names, skip = 1,
                 id = "file_path")
  df <- df %>%
    mutate(OF = str_extract(file_path, "OF_[0-9]*"),
           OF = str_remove(OF, "OF_"),
           FCP = str_extract(file_path, "FCP_[0|1]*[e|.][-]*[0-9]*"),
           FCP = str_remove(FCP, "FCP_"),
           FCP = ifelse(FCP=="0.01", "1e-02", FCP),
           run = str_extract(file_path, "run [0-9]* "),
           run = str_trim(str_remove(run, "run"))) %>% 
    select(-file_path)
  return(df)
}

get_reduced_data <- function(dat, supressed = FALSE){
  dat <- dat %>% 
    group_by(OF, FCP, run) 
  if(supressed){
    reduced_dat <- dat %>% 
      filter(wAlbAB_prop < 0.1)
  }else{
    reduced_dat <- dat %>% 
      filter(wAlbAB_adults == 0)
  }
  
  reduced_dat <- reduced_dat %>% 
    slice(1) %>% 
    ungroup()
  
  return(reduced_dat)
}


# parameters --------------------------------------------------------------

OF <- c(1, 5, 15)
FCP <- c(1e-2, 1e-4)

label_df <- expand_grid(OF, FCP)

folder_paths <- here::here(glue::glue("outputs/scenario_OF_{label_df$OF}_FCP_{label_df$FCP}"))
file_paths <- map(folder_paths, list.files, full.names=TRUE) %>% 
  unlist()

# Find the initial population
param_list <- readRDS(here::here("data/param_list.Rds"))
m_f_ind <- getIndicesForStateQuery(stateNames = names(param_list$state_vector),
                                   sex = c("m", "f"), genotype = "wAlbAB")
init_pop <- sum(param_list$state_vector[m_f_ind])


full_data <- map_dfr(file_paths, get_data) %>% 
  mutate(wAlbAB_adults = block01_wAlbAB_m + block01_wAlbAB_f,
         wPip_adults = block01_wPip_m + block01_wPip_f,
         wAlbAB_prop = wAlbAB_adults/init_pop,
         wPip_prop = wPip_adults/init_pop) %>% 
  select(-contains("block")) 
write_rds(full_data, here::here("data/all_scenarios_all_runs.Rds"))

supressed_data <- get_reduced_data(full_data, supressed = TRUE)
eliminated_data <- get_reduced_data(full_data, supressed = FALSE)

write_rds(supressed_data,  here::here("data/all_scenarios_supressed.Rds"))
write_rds(eliminated_data,  here::here("data/all_scenarios_eliminated.Rds"))

rm("full_data", "supressed_data", "eliminated_data")
