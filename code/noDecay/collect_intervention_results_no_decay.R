## Matt Ryan
## 09/02/2024
### Collect the data from the intervention simulations
### Requires: tidyverse, pacman, here
### TO DO: Make code more general for other parameter choices


# packages ----------------------------------------------------------------

pacman::p_load(tidyverse)
mc.cores <- 5

# Functions ---------------------------------------------------------------

get_data <- function(f){
  column_names <- c("days", "block01_wAlbAB_m", "block01_wAlbAB_f", "block01_wAlbAB_imm",
                    "block01_wPip_m", "block01_wPip_f", "block01_wPip_imm", "total_released_to_date")
  df <- read_csv(f,
                 col_types = cols(), col_names = column_names, skip = 1,
                 id = "file_path")
  df <- df %>%
    mutate(immigration = str_extract(file_path, "femaleImmigration_[0-9]*[.]*[0-9]*"),
           immigration = str_remove(immigration, "femaleImmigration_"),
           controlProp = str_extract(file_path, "controlProportion_0.[0-9]*"),
           controlProp = str_remove(controlProp, "controlProportion_"), 
           intervention = str_extract(file_path, "interventionType_[a-z]*"),
           intervention = str_remove(intervention, "interventionType_"),
           run = str_extract(file_path, "run[0-9]*"),
           run = str_trim(str_remove(run, "run"))) %>% 
    select(-file_path)
  
  return(df)
}


# Load and clean data -----------------------------------------------------

# Get directories
dirs <- list.files(here::here("outputs/nodecay_intervention_simulations"))
dirs <- dirs[-1] # get rid of archive

param_combs <- c("Expected")

df <- list()

for(idx in 1:length(dirs)){
  # List all directories
  dr <- list.dirs(here::here(glue::glue("outputs/nodecay_intervention_simulations/{dirs[idx]}/")), 
                  full.names = T)
  dr <- dr[-1]
  
  dat <- parallel::mclapply(dr,
             function(d){
               fl <- list.files(d, full.names = TRUE)
               dat <- map_df(fl, get_data)
               dat <- dat %>% 
                 mutate(wAlb_adult = block01_wAlbAB_m + block01_wAlbAB_f,
                        wPip_adult = block01_wPip_m + block01_wPip_f,
                        wAlb = wAlb_adult + block01_wAlbAB_imm,
                        wPip = wPip_adult + block01_wPip_imm,
                        params = param_combs[idx])
               
               # Rename interventions
               dat <- dat %>%
                 mutate(intervention=ifelse(intervention=="elimination", "naive", intervention),
                        intervention=ifelse(intervention=="complete", "complete stop", intervention),
                        intervention=factor(intervention, levels=c("naive", "complete stop", "maintain")))
               
               return(dat)
             },
             mc.cores = mc.cores)
  
  dat <- do.call(bind_rows, dat)
  
  
  
  write_rds(dat, here::here(glue::glue("data/intervention_results_nodecay_{param_combs[idx]}.Rds")))
}

