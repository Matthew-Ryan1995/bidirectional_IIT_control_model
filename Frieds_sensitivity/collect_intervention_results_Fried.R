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
           run = str_trim(str_remove(run, "run")),
           fried = str_extract(file_path, "Frieds_[0-9]*.[0-9]*"),
           fried = str_remove(fried, "Frieds_")) %>% 
    select(-file_path)
  
  return(df)
}


# Load and clean data -----------------------------------------------------

# Get directories
dirs <- list.files(here::here("outputs/Fried/intervention_simulations"))
dirs <- dirs[-1] # get rid of archive

Fried_values <- c(
  "low" = "0.36", # Caputo 2023 ( DOI 10.1002/ps.7495)
  "mid" = "0.73" # Caputo 2023 ( DOI 10.1002/ps.7495)
  # "high" = "1.2", # Caputo 2019 (DOI 10.1002/ps.5643), Puggioli 2016 (https://doi.org/10.1016/j.actatropica.2016.10.014)
  # "extreme" = "1.71" # Moretti 2018 (doi: 10.1371/journal.pntd.0006626)
)

df <- list()

for(idx in 1:length(dirs)){
  # List all directories
  dr <- list.dirs(here::here(glue::glue("outputs/Fried/intervention_simulations/{dirs[idx]}/")), 
                  full.names = T)
  dr <- dr[-1]
  
  save_label <- case_when(
    str_detect(dirs[idx], Fried_values[1]) ~ names(Fried_values)[1],
    str_detect(dirs[idx], Fried_values[2]) ~ names(Fried_values)[2],
    str_detect(dirs[idx], Fried_values[3]) ~ names(Fried_values)[3],
    str_detect(dirs[idx], Fried_values[4]) ~ names(Fried_values)[4],
    TRUE ~ NA
  )
  
  dat <- parallel::mclapply(dr,
             function(d){
               fl <- list.files(d, full.names = TRUE)
               dat <- map_df(fl, get_data)
               dat <- dat %>% 
                 mutate(wAlb_adult = block01_wAlbAB_m + block01_wAlbAB_f,
                        wPip_adult = block01_wPip_m + block01_wPip_f,
                        wAlb = wAlb_adult + block01_wAlbAB_imm,
                        wPip = wPip_adult + block01_wPip_imm,
                        params = save_label)
               
               # Rename interventions
               dat <- dat %>%
                 mutate(intervention=ifelse(intervention=="elimination", "naive", intervention),
                        intervention=ifelse(intervention=="complete", "complete stop", intervention),
                        intervention=factor(intervention, levels=c("naive", "complete stop", "maintain")))
               
               return(dat)
             },
             mc.cores = mc.cores)
  
  dat <- do.call(bind_rows, dat)
  
  
  write_rds(dat, here::here(glue::glue("data/intervention_results_Fried_{save_label}.Rds")))
}

