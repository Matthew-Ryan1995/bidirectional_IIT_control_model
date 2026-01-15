## 13/08/2024
## Figures for Matt's IIT paper
## Based on his code from cage simulations
## create_cage_result_plots.R
## Requires: pacman, here, tidyverse
## author: Manuela Mendiolar


# Figure 2
# Percentage of simulations with establishment of WAlbAB and wPip 
# for bi and uni directions


# packages ----------------------------------------------------------------

pacman::p_load(tidyverse, patchwork)


# data and parameters -----------------------------------------------------

# initial population size (wildtype)
popsize.ini <- 420 

# proportion of the initial population size where the pop becomes established
prop.establish <- 0.1 

# basic setup for plots
height <- 8
text_size <- 20

col.Alb <- "gray60"
col.wAlbAB <- "royalblue3"
col.wPip <- "red3"

# lab.Alb <- expression(paste(Alb," (",italic("Wolbachia "),"cleared out)"))
lab.wAlbAB <- expression(paste(italic("w"),AlbAB))
lab.wPip <-  expression(paste(AR,italic(w),P))

# plots -------------------------------------------------------------------
param_combs <- c("Low", "Expected", "High")
# for(pp in unique(param_combs)){
  
  dat_full <- list(
  df_expected = read_rds(here::here(glue::glue("data/cage_experiment_results_Expected.Rds"))) %>%
    filter(direction=="bi") %>%
    mutate(fried = "1.0"),
  df_low = read_rds(here::here("data/cage_experiment_results_Fried_low.Rds")),
  df_mid = read_rds(here::here("data/cage_experiment_results_Fried_mid.Rds"))
  )
  dat_full <- do.call(bind_rows, dat_full)
  
  dat_full <- dat_full %>%
    filter(!(propPip %in% c(0.22, 0.8)))

  # Find reduced dataset where sims stopped
  dat_endtime <- dat_full %>%
    filter(wPip >0 | wAlb > 0) %>%
    group_by(run, propPip, direction, params) %>%
    slice_max(days) %>%
    ungroup()

  ## Proportion established --------------------------------------------------
  
  # CI: uni-directional
  
  # prop_uni <- dat_endtime %>% 
  #   pivot_longer(c(wPip_adult, wAlb_adult)) %>% 
  #   mutate(name = ifelse(name=="wPip_adult", "wPip", "wAlbAB_cleared")) %>%
  #   group_by(propPip, name, direction) %>% 
  #   count(established = value > prop.establish*popsize.ini) %>% 
  #   mutate(n = (n/sum(n)*100)) %>%
  #   ungroup() %>% 
  #   mutate(n=ifelse(established, n, 0)) %>% 
  #   select(-established) %>% 
  #   filter(direction %in% c("uni")) 
  
  prop_bi <- dat_endtime %>% 
    pivot_longer(c(wPip_adult, wAlb_adult)) %>% 
    mutate(name = ifelse(name=="wPip_adult", "wPip", "wAlbAB")) %>% 
    group_by(propPip, name, direction, fried) %>% 
    count(established = value > prop.establish*popsize.ini) %>% 
    mutate(n = (n/sum(n)*100)) %>%
    ungroup() %>% 
    mutate(n=ifelse(established, n, 0)) %>% 
    select(-established) %>% 
    filter(direction %in% c("bi")) %>% 
    mutate(name = factor(name, levels = c("wAlbAB", "wPip")),
           fried = factor(fried,
                          levels = c("1.0", "0.73", "0.36"),
                          labels = c("Fried's Index: 1.0", "Fried's Index: 0.73", "Fried's Index: 0.36")),
           fried = fct_rev(fried))
  
  # prop_bi_no_decay <- dat_endtime %>% 
  #   pivot_longer(c(wPip_adult, wAlb_adult)) %>% 
  #   mutate(name = ifelse(name=="wPip_adult", "wPip", "wAlbAB")) %>% 
  #   group_by(propPip, name, direction) %>% 
  #   count(established = value > prop.establish*popsize.ini) %>% 
  #   mutate(n = (n/sum(n)*100)) %>%
  #   ungroup() %>% 
  #   mutate(n=ifelse(established, n, 0)) %>% 
  #   select(-established) %>% 
  #   filter(direction %in% c("bi_no_decay"))
  
  # prop_all <- bind_rows(prop_uni, prop_bi, prop_bi_no_decay) %>%
  #   mutate(
  #     name = factor(name, levels = c("wAlbAB_cleared", "wAlbAB", "wPip")),
  #     direction = factor(
  #       direction, 
  #       levels = c("uni", "bi", "bi_no_decay"), 
  #       labels = c("Uni-directional CI", "Bi-directional CI", "Bi-directional CI (no age-decay)")
  #     )
  #   )
  
  # Plot
  p_props <- ggplot(prop_bi, aes(x = propPip, y = n, fill = name)) +
    geom_col(colour="black", position = "dodge") +
    geom_hline(yintercept = 50, linetype = "dashed") +
    facet_wrap(~ direction, nrow = 1) +
    scale_y_continuous(limits = c(0, 100)) +
    scale_x_discrete(labels = c(" ","0.10"," ","0.20"," ","0.30"," ","0.40"," ","0.50" )) +
    scale_fill_manual(
      values = c("wAlbAB" = col.wAlbAB, "wPip" = col.wPip),
      labels = c("wAlbAB" = lab.wAlbAB, "wPip" = lab.wPip),
      name = "Population"
    ) +
    labs( 
      x = bquote("Initial proportion of " ~ .(lab.wPip[[1]])),
      y = "Percentage of simulations\nwith establishment"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      axis.title.x = element_text(size = 18), 
      axis.title.y = element_text(size = 18), 
      axis.text.x = element_text(size = 14.5, hjust = 0.5),
      axis.text.y = element_text(size = 14.5),
      strip.text = element_text(size = 17.5),
      text = element_text(family = "Times New Roman")
    ) +
    facet_wrap(~fried)
  
  # Save
  ggsave(here::here(glue::glue("img/fried_cage_prop_establishment.png")),
         plot = p_props,
         height = 6, width = 14, dpi = 300)
  
# }