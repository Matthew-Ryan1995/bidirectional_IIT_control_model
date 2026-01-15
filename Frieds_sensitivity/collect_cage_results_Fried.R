## Matt Ryan
## 09/02/2024
### Collect the data from the cage simulations
### Requires: tidyverse, pacman, here
### TO DO: Make code more general for other parameter choices


# packages ----------------------------------------------------------------

pacman::p_load(tidyverse)
mc.cores <- 5

# functions ---------------------------------------------------------------

get_data <- function(f){
  column_names <- c("days", "block01_wAlbAB_m", "block01_wAlbAB_f", "block01_wAlbAB_imm",
                    "block01_wPip_m", "block01_wPip_f", "block01_wPip_imm", "total_released_to_date")
  df <- read_csv(f,
                 col_types = cols(), col_names = column_names, skip = 1,
                 id = "file_path")
  
  idd <- which((rowSums(df[, -c(1:2)]>0))==0)[1] - 1
  
  if(is.na(idd)){
    idd <- nrow(df)
  }
  if(idd < 1){
    idd <- 1
  }
  df <- df %>% 
    slice(1:idd)
  
  # original ---------------------------------------------------
  df <- df %>%
    mutate(propPip = str_extract(file_path, "propPip_0.[0-9]*"),
           propPip = str_remove(propPip, "propPip_"),
           femaleProportion = str_extract(file_path, "femaleProportion_0.[0-9]*"),
           femaleProportion = str_remove(femaleProportion, "femaleProportion_"),
           # original ----------------------------------------------
           # direction = str_extract(file_path, "direction_[a-z]*"),
           # direction = str_remove(direction, "direction_"),
           # -------------------------------------------------------
           direction = str_match(file_path, "direction_\\s*(.*?)\\s*_start_")[ ,2],
           run = str_extract(file_path, "run[0-9]*"),
           run = str_trim(str_remove(run, "run")),
           fried = str_extract(file_path, "Frieds_[0-9]*.[0-9]*"),
           fried = str_remove(fried, "Frieds_")) %>%
    select(-file_path)

  df <- df %>% 
    mutate(pAlb = apply(df[, 2:7], 1, function(x) sum(x[1:3])/sum(x)),
           pAlb = ifelse(is.nan(pAlb), NA, pAlb),
           pPip = 1-pAlb,
           pAlb_female = (block01_wAlbAB_f)/(block01_wAlbAB_f + block01_wPip_f),
           pAlb_female = ifelse(is.nan(pAlb_female), NA, pAlb_female),
           pPip_female = 1-pAlb_female) 
  
  return(df)
}


# Get fdata ---------------------------------------------------------------

# Get directories
dirs <- list.files(here::here("outputs/Fried/cage_simulations"))
# dirs <- dirs[-1] # get rid of archive

Fried_values <- c(
  "low" = "0.36", # Caputo 2023 ( DOI 10.1002/ps.7495)
  "mid" = "0.73" # Caputo 2023 ( DOI 10.1002/ps.7495)
  # "high" = "1.2", # Caputo 2019 (DOI 10.1002/ps.5643), Puggioli 2016 (https://doi.org/10.1016/j.actatropica.2016.10.014)
  # "extreme" = "1.71" # Moretti 2018 (doi: 10.1371/journal.pntd.0006626)
)

df <- list()

for(idx in 1:length(dirs)){
  dr <- list.dirs(here::here(glue::glue("outputs/Fried/cage_simulations/{dirs[idx]}/")), full.names = T)
  # dr <- dr[-1]
  
  save_label <- case_when(
    str_detect(dirs[idx], Fried_values[1]) ~ names(Fried_values)[1],
    str_detect(dirs[idx], Fried_values[2]) ~ names(Fried_values)[2],
    str_detect(dirs[idx], Fried_values[3]) ~ names(Fried_values)[3],
    str_detect(dirs[idx], Fried_values[4]) ~ names(Fried_values)[4],
    TRUE ~ NA
  )
  
  # Get and load files
  print("start_collect_data")
  dat <- parallel::mclapply(dr,
                            function(d){
                              fl <- list.files(d, full.names = TRUE)
                              dat <- map_dfr(fl, get_data)
                              dat <- dat %>% 
                                mutate(wAlb_adult = block01_wAlbAB_m + block01_wAlbAB_f,
                                       wPip_adult = block01_wPip_m + block01_wPip_f,
                                       wAlb = wAlb_adult + block01_wAlbAB_imm,
                                       wPip = wPip_adult + block01_wPip_imm,
                                       params = save_label)
                              return(dat)
                            },
                            mc.cores = mc.cores)
  
  dat <- do.call(bind_rows, dat)
  
  write_rds(dat, here::here(glue::glue("data/cage_experiment_results_Fried_{save_label}.Rds")))
  
}
 
