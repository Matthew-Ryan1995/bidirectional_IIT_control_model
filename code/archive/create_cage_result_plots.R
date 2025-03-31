# Matt Ryan
# 09/02/2024
## Create figures of results from cage simulations
## Requires: pacman, here, tidyverse


# packages ----------------------------------------------------------------

pacman::p_load(tidyverse)


# data and parameters -----------------------------------------------------

# dat <- read_rds(here::here("data/cage_experiment_results_expected_parameters.Rds"))
dat_full <- read_rds(here::here("data/cage_experiment_results_full.Rds"))

dat_full <- dat_full %>% 
  filter(!(propPip %in% c(0.22, 0.8)))

height <- 8
text_size <- 20

# Find reduced dataset where sims stopped
dat_endtime_full <- dat_full %>% 
  filter(wPip >0 | wAlb > 0) %>% 
  group_by(run, propPip, direction, params) %>% 
  slice_max(days) %>% 
  ungroup() 


# plots -------------------------------------------------------------------

for(pp in unique(dat_full$params)){
  
  dat <- dat_full %>% 
    filter(params==pp)
  dat_endtime <- dat_endtime_full %>% 
    filter(params==pp)
  
  lifespan <- case_when(
    pp=="High" ~ "Large population parameters",
    pp=="Low" ~ "Small population parameters",
    TRUE ~ "Expected population parameters"
  )
  
  ## Time until wAlb suppressed ----------------------------------------------
  
  p_stop <- dat_endtime %>% 
    ggplot(aes(x = as.factor(propPip), y=days, fill=direction)) +
    geom_boxplot() +
    labs(x = "Initial proportion of wPip", 
         y = "Days until simulation stopped",
         fill="Direction of CI",
         # title = str_c("Parameter set: ", pp, "\n", lifespan)
         ) +
    harrypotter::scale_fill_hp_d("ravenclaw") +
    theme_bw() + 
    theme(legend.position = "bottom",
          text = element_text(size = text_size)) +
    coord_flip()
  
  ggsave(here::here(glue::glue("img/cage_time_until_stop_params_{pp}.png")), plot = p_stop,
         height = height, width = 1.618*height, dpi = 300)
  
  ## Proportion established --------------------------------------------------
  
  p_props <- dat_endtime %>% 
    pivot_longer(c(wPip_adult, wAlb_adult)) %>% 
    mutate(name = ifelse(name=="wPip_adult", "wPip", "wAlbAB")) %>% 
    group_by(propPip, name, direction) %>% 
    count(established = value > 0.1*420) %>% 
    mutate(n = (n/sum(n)*100)) %>%
    ungroup() %>% 
    mutate(n=ifelse(established, n, 0)) %>% 
    select(-established) %>% 
    # pivot_wider(names_from=name, values_from = n, values_fill = 0) %>% 
    ggplot(aes(x=propPip, y = n, fill=direction)) +
    geom_col(position="dodge") +
    geom_hline(yintercept = 50, lty=2) +
    facet_wrap(~name, ncol = 1) +
    labs(x = "Initial proportion of wPip", 
         y = "Percentage of simulations with establishment",
         fill="Direction of CI", 
         # caption="Dotted line = 50% of simulation",
         # title = str_c("Parameter set: ", pp, "\n", lifespan)
         ) +
    harrypotter::scale_fill_hp_d("ravenclaw") +
    theme_bw() + 
    theme(legend.position = "bottom",
          text = element_text(size = text_size))
  
  ggsave(here::here(glue::glue("img/cage_prop_establishment_params_{pp}.png")), plot = p_props,
         height = height, width = 1.618*height, dpi = 300)
  
  
  # Stormcloud --------------------------------------------------------------
  
  plot_data <- dat %>% 
    filter(wPip > 0 | wAlb > 0) %>% 
    mutate(wAlb_adult = wAlb_adult/420,
           wPip_adult = wPip_adult/420)
  plot_data_medians <- plot_data %>% 
    group_by(days, propPip, direction) %>% 
    ggdist::median_qi(.exclude=c("run", "femaleProportion", "params")) %>% 
    ungroup()
  
  p_storm <- plot_data %>% 
    filter(propPip %in% c("0.1", "0.2", "0.3", "0.4")) %>% 
    ggplot(aes(x=days)) +
    geom_line(aes(y=wAlb_adult, group = run), colour="black", alpha=0.1) +
    geom_line(data=plot_data_medians %>% 
                filter(propPip %in% c("0.1", "0.2", "0.3", "0.4")), aes(y=wAlb_adult), colour = "blue") +
    geom_ribbon(data=plot_data_medians %>% 
                  filter(propPip %in% c("0.1", "0.2", "0.3", "0.4")), aes(ymin=wAlb_adult.lower, ymax=wAlb_adult.upper), 
                fill = "blue", alpha=0.3) +
    geom_line(aes(y=wPip_adult, group = run), colour="black", alpha=0.1)+
    geom_line(data=plot_data_medians %>% 
                filter(propPip %in% c("0.1", "0.2", "0.3", "0.4")), aes(y=wPip_adult), colour = "red") +
    geom_ribbon(data=plot_data_medians %>% 
                  filter(propPip %in% c("0.1", "0.2", "0.3", "0.4")), aes(ymin=wPip_adult.lower, ymax=wPip_adult.upper), 
                fill = "red", alpha=0.3) +
    geom_hline(yintercept = 0.1, lty=2) +
    facet_grid(propPip~direction, labeller = "label_both") +
    theme_bw() +
    ylim(0, 1) +
    labs(y="Proportion of steady state", 
         # caption="Blue=wAlbAB\nred=wPip",
         # title = str_c("Parameter set: ", pp, "\n", lifespan)
         ) +
    theme(text = element_text(size = text_size))
  
  ggsave(here::here(glue::glue("img/cage_stormcloud_params_{pp}.png")), plot = p_storm,
         height = height, width = 1.618*height, dpi = 300)
}
