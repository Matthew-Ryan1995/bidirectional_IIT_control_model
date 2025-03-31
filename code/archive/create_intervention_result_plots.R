# Matt Ryan
# 09/02/2024
## Create figures of results from intervention simulations
## Requires: pacman, here, tidyverse


# packages ----------------------------------------------------------------

pacman::p_load(tidyverse)

# data and parameters -----------------------------------------------------

# dat <- read_rds(here::here("data/intervention_results_expected_parameters.Rds"))
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

# Find reduced dataset where sims stopped
dat_endtime_full <- dat_full %>% 
  filter(wPip >0 | wAlb > 0) %>% 
  group_by(run, intervention, immigration, params) %>% 
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
  
  #todo: Don't hardcode carrying capacity
  cc <- case_when(
    pp=="High" ~800,
    pp=="Low" ~400,
    TRUE ~ 420
  )
  
  
  ## Time until stop  -------------------------------------------------------
  
  p_stop <- dat_endtime %>% 
    ggplot(aes(x = as.factor(intervention), y=days, fill=as.factor(immigration))) +
    geom_boxplot() +
    labs(y = "Days until simulation stopped",
         fill="Immigration rate",
         # title = str_c("Parameter set: ", pp, "\n", lifespan)
         ) +
    harrypotter::scale_fill_hp_d("ravenclaw", direction = -1) +
    theme_bw() + 
    theme(legend.position = "bottom",
          text = element_text(size = text_size)) +
    coord_flip()
  
  ggsave(here::here(glue::glue("img/interventions_time_until_stop_{pp}.png")), plot = p_stop,
         height = height, width = 1.618*height, dpi = 300)
  
  ## Establishment proportions -----------------------------------------------
  
  p_props <- dat_endtime %>% 
    pivot_longer(c(wPip_adult, wAlb_adult)) %>% 
    mutate(name = ifelse(name=="wPip_adult", "wPip", "wAlbAB")) %>% 
    group_by(intervention, immigration, name) %>% 
    count(established = value > 0.1*cc) %>%
    mutate(n = (n/sum(n)*100)) %>%
    ungroup() %>% 
    mutate(n=ifelse(established, n, 0)) %>%
    select(-established) %>% 
    ggplot(aes(x=immigration, y = n, fill=intervention)) +
    geom_col(position = "dodge") +
    geom_hline(yintercept = 50, lty=2) +
    facet_wrap(~name, ncol = 1) +
    labs(x = "Immigration rate", 
         y = "Percentage of simulations with establishment",
         fill="Intervention type", 
         # caption="Dotted line = 50% of simulation",
         # title = str_c("Parameter set: ", pp, "\n", lifespan)
         ) +
    harrypotter::scale_fill_hp_d("ravenclaw", direction = -1) +
    theme_bw() + 
    theme(legend.position = "bottom",
          text = element_text(size = text_size))
  
  ggsave(here::here(glue::glue("img/interventions_prop_establishment_{pp}.png")), plot = p_props,
         height = height, width = 1.618*height, dpi = 300)
  
  # Stormclouds -------------------------------------------------------------
  
  ## Create the median curves
  plot_data <- dat %>% 
    filter(wPip > 0 | wAlb > 0) %>% 
    mutate(wAlb_adult = wAlb_adult/cc,
           wPip_adult = wPip_adult/cc) # Change everything to proportion land
  plot_data_medians <- plot_data %>% 
    group_by(days, intervention, immigration) %>% 
    ggdist::median_qi(.exclude=c("run", "controlProp", "params")) %>% 
    ungroup()
  
  p_storm <- plot_data %>% 
    ggplot(aes(x=days)) +
    geom_line(aes(y=wAlb_adult, group = run), colour="black", alpha=0.1) +
    geom_line(data=plot_data_medians, aes(y=wAlb_adult), colour = "blue") +
    geom_ribbon(data=plot_data_medians, aes(ymin=wAlb_adult.lower, ymax=wAlb_adult.upper), fill = "blue", alpha=0.3) +
    geom_line(aes(y=wPip_adult, group = run), colour="black", alpha=0.1)+
    geom_line(data=plot_data_medians, aes(y=wPip_adult), colour = "red") +
    geom_ribbon(data=plot_data_medians, aes(ymin=wPip_adult.lower, ymax=wPip_adult.upper), fill = "red", alpha=0.3) +
    geom_hline(yintercept = 0.1, lty=2) +
    geom_hline(yintercept = 0.4, lty=2, colour = "gray70") +
    facet_grid(intervention ~ immigration, labeller = "label_both") +
    theme_bw() +
    labs(y="Proportion of steady state population", 
         x = "Days",
         # caption="Blue=wAlbAB\nred=wPip",
         # title = str_c("Parameter set: ", pp, "\n", lifespan)
         ) +
    theme(text = element_text(size = text_size))
  
  ggsave(here::here(glue::glue("img/interventions_stormcloud_{pp}.png")), plot = p_storm,
         height = 11, width = 1.618*11, dpi = 300)
  
  
  
  # cost - boxplot ----------------------------------------------------------
  
  p_cost_box <- dat %>% 
    group_by(run, intervention, immigration) %>% 
    slice_max(total_released_to_date) %>% 
    slice(1) %>% 
    ggplot(aes(x = total_released_to_date, y = immigration, fill = intervention)) +
    geom_boxplot() + 
    labs(x = "Total adult mosquitos released", y = "Immigration rate",
         # caption = "",
         fill="Intervention type",
         # title = str_c("Parameter set: ", pp, "\n", lifespan)
         ) +
    harrypotter::scale_fill_hp_d("ravenclaw", direction=-1) +
    theme_bw() + 
    theme(legend.position = "bottom",
          text = element_text(size = text_size))
  
  ggsave(here::here(glue::glue("img/interventions_cost_boxplot_{pp}.png")), plot = p_cost_box,
         height = height, width = 1.618*height, dpi = 300)
  
  
  # cost - stormcloud -------------------------------------------------------
  
  p_cost_storm <- plot_data %>% 
    ggplot(aes(x=days, colour=intervention, fill = intervention)) +
    geom_line(aes(y=total_released_to_date, group = run), colour="black", alpha=0.1) +
    geom_line(data=plot_data_medians, aes(y=total_released_to_date)) +
    geom_ribbon(data=plot_data_medians, aes(ymin=total_released_to_date.lower, 
                                            ymax=total_released_to_date.upper), 
                colour=NA, alpha=0.3, show.legend = FALSE) +
    facet_wrap( ~ immigration, labeller = "label_both") +
    theme_bw() +
    harrypotter::scale_colour_hp_d("ravenclaw", direction = -1) +
    harrypotter::scale_fill_hp_d("ravenclaw", direction = -1) +
    labs(y="Total wPip released", 
         colour = "Intervention type",
         # title = str_c("Parameter set: ", pp, "\n", lifespan)
         ) +
    theme(legend.position = "bottom",
          text = element_text(size = text_size))
  
  ggsave(here::here(glue::glue("img/interventions_cost_stormcloud_{pp}.png")), plot = p_cost_storm,
         height = height, width = 1.618*height, dpi = 300)
  
}