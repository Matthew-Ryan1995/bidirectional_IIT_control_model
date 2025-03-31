# Matt Ryan
# 09/02/2024
## Create figures of results from intervention simulations
## Requires: pacman, here, tidyverse


# Figure 4
# Percentage of simulations with "success" (no establishment of wAlbAB and wPip ?) 
# according to intervention type and 0, 2 or 10 per week (Immigration rate) 



# packages ----------------------------------------------------------------

pacman::p_load(tidyverse)

# data and parameters -----------------------------------------------------

dat_full <- read_rds(here::here("data/intervention_results_full.Rds"))

dat_full <- dat_full %>% 
  mutate(immigration = case_when(
    immigration=="0" ~ "0 per week",
    immigration=="0.286" ~ "2 per week",
    TRUE ~ "10 per week"
  ),
  immigration = factor(immigration, levels = c("0 per week", "2 per week", "10 per week")))

height <- 8
text_size <- 20

lab.Alb <- expression(paste(Alb," (",italic("Wolbachia "),"cleared out)"))
lab.wAlbAB <- expression(paste(italic("w"),AlbAB))
lab.wPip <- expression(paste(italic(w),Pip))


# Find reduced dataset where sims stopped
# dat_endtime_full <- dat_full %>% 
#   filter(wPip > 0 | wAlb > 0) %>% 
#   filter(days <= 365*3) %>% 
#   group_by(run, intervention, immigration, params) %>% 
#   slice_max(days) %>% 
#   ungroup() 



# plots -------------------------------------------------------------------

for(pp in unique(dat_full$params)){
  
  # pp <- "Expected"
  
  dat <- dat_full %>% filter(params==pp)

  # dat_endtime <- dat_endtime_full %>% 
  #   filter(params==pp)

  lifespan <- case_when(
    pp == "High" ~ "Large population parameters",
    pp == "Low" ~ "Small population parameters",
    TRUE ~ "Expected population parameters"
  )

  #todo: Don't hardcode carrying capacity
  cc <- case_when(
    pp == "High" ~ 800,
    pp == "Low" ~ 400,
    TRUE ~ 420
  )
  
  # Compute stop time for each intervention, immigration and run combination
  
  # add row indices to the original dataframe
  aux <- dat %>%
    mutate(row_index = row_number())
 
  # compute the row index of the second-to-last time 'total_released_to_date' changed within each group
  dat_stoptime <- aux %>%
    group_by(intervention, immigration, run) %>%
    mutate(change_flag = total_released_to_date != lag(total_released_to_date, default = first(total_released_to_date))) %>%
    filter(change_flag == TRUE) %>%
    summarise(antepenultimate_change_row = nth(row_index, -2)) %>%
    ungroup()
  
  # add a 'days' column based on 'antepenultimate_change_row'
  dat_stoptime <- dat_stoptime %>%
    left_join(aux %>% select(row_index, days), by = c("antepenultimate_change_row" = "row_index")) %>% 
    mutate(days = days + 7) # add extra week as they first check how things and a week later the final decision is made
  

  tmp <- left_join(dat, rename(dat_stoptime, endtime=days), by=c("run", "intervention", "immigration"))
  tmp <- tmp %>% mutate(after_6_months = endtime + 180)
  
 
  # Function to compute success at different times (stop time and after 6 months) for both wAlbAB and wPip
  calculate_establishment <- function(data, strain, threshold, time_var) {
    data %>%
      group_by(run, intervention, immigration) %>%
      filter(days == !!sym(time_var)) %>%
      pivot_longer(c(wPip_adult, wAlb_adult)) %>%
      mutate(name = ifelse(name=="wPip_adult", "wPip", "wAlbAB")) %>% 
      filter(name == strain) %>%
      group_by(intervention, immigration, name) %>%
      count(established = value > threshold * cc) %>%
      mutate(n = (n / sum(n) * 100)) %>%
      ungroup() %>%
      mutate(n = ifelse(established, n, 0)) %>%
      select(-established) %>%
      group_by(across(-n)) %>%
      summarise(n = sum(n)) %>%
      mutate(n_success = 100 - n, time = time_var)
  }
  
  ## Establishment proportions -----------------------------------------------
  dat_wAlbAB_stop <- calculate_establishment(tmp, "wAlbAB", 0.1, "endtime")
  dat_wPip_stop <- calculate_establishment(tmp, "wPip", 0.4, "endtime")
  dat_wAlbAB_6months <- calculate_establishment(tmp, "wAlbAB", 0.1, "after_6_months")
  dat_wPip_6months <- calculate_establishment(tmp, "wPip", 0.4, "after_6_months")
  
  # combine data sets
  dat_both <- bind_rows(dat_wAlbAB_stop, dat_wPip_stop, dat_wAlbAB_6months, dat_wPip_6months)
  
  # ensure 'time' is a factor 
  dat_both$time <- factor(dat_both$time, levels = c("endtime", "after_6_months"))
  
  
  # Updated the facet labels for better clarity
  p_all <- ggplot(dat_both, aes(x = immigration, y = n_success, fill = intervention)) +
    geom_col(position = "dodge") +
    geom_hline(yintercept = 50, lty = 2) +
    facet_grid(name ~ time,
               labeller = labeller(time = c("endtime" = "Within stopping time", 
                                            "after_6_months" = "After six months"),
                                   name = c("wAlbAB" = lab.wAlbAB, "wPip" = lab.wPip)
                                   )
    ) +
    labs(x = "Immigration rate", 
         y = "Percentage of simulations with management success",
         fill = "Intervention type", 
         title = "") +
    harrypotter::scale_fill_hp_d("ravenclaw", direction = -1, labels = c("Naïve", "Complete stop", "Maintain")) +
    theme_bw() + 
    theme(legend.position = "bottom",
          plot.title = element_text(size = 18),
          text = element_text(size = text_size))
  
  ggsave(here::here(glue::glue("img/interventions_prop_establishment_{pp}_MM.png")), 
         plot = p_all,
         height = 1.5*height, width = 1.75*height, dpi = 300)
}

