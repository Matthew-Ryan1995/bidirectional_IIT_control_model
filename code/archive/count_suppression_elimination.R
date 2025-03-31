## Matt Ryan
## 14-11-2023
## Count number that reach elimination and suppression
## Requires: tidyverse


# libraries ---------------------------------------------------------------

pacman::p_load(tidyverse)


# data --------------------------------------------------------------------

df <- read_rds("data/all_scenarios_all_runs.Rds")


# elimination -------------------------------------------------------------

df %>% 
  group_by(run, FCP, OF) %>% 
  filter(wAlbAB_adults == 0) %>% 
  slice(1) %>% 
  ungroup() %>% 
  count(FCP, OF)


# suppression -------------------------------------------------------------

df %>% 
  group_by(run, FCP, OF) %>% 
  filter(wAlbAB_prop < 0.1) %>% 
  slice(1) %>% 
  ungroup() %>% 
  count(FCP, OF)
