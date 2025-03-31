# Matt Ryan
# 09/02/2024
## Create figures of results from intervention simulations
## Requires: pacman, here, tidyverse



# Figure 3
# Stormclouds + Zoom


# packages ----------------------------------------------------------------

pacman::p_load(tidyverse, patchwork, posterior)

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

col.simruns <- "black"
col.medianwPip <- "red"
col.medianwAlbAB <- "blue"


lab.Alb <- expression(paste(Alb," (",italic("Wolbachia "),"cleared out)"))
lab.wAlbAB <- expression(paste(italic("w"),AlbAB))
lab.wPip <-  expression(paste(AR,italic(w),P))

lab.intervention <- c(`naive` = "Intervention:\nnaïve",
                      `complete stop` = "Intervention:\ncomplete stop",
                      `maintain` = "Intervention:\nmaintain")

lab.immigration <- c(`0 per week` = "Immigration: 0 per week",
                     `2 per week` = "Immigration: 2 per week",
                     `10 per week` = "Immigration: 10 per week")

# plots -------------------------------------------------------------------
param_combs <- c("Expected")
for(pp in unique(param_combs)){
  
  # pp <- "Expected" 
  dat_full <- read_rds(here::here(glue::glue("data/intervention_results_nodecay_{pp}.Rds")))
  
  dat <- dat_full %>% 
    mutate(immigration = case_when(
      immigration=="0" ~ "0 per week",
      immigration=="0.286" ~ "2 per week",
      TRUE ~ "10 per week"
    ),
    immigration = factor(immigration, levels = c("0 per week", "2 per week", "10 per week")))
  
  
  # Find reduced dataset where sims stopped
  dat_endtime <- dat_full %>% 
    filter(wPip >0 | wAlb > 0) %>% 
    group_by(run, intervention, immigration, params) %>% 
    slice_max(days) %>% 
    ungroup() 
  
  
  
  lifespan <- case_when(
    pp=="High" ~ "Large population parameters",
    pp=="Low" ~ "Small population parameters",
    TRUE ~ "Expected population parameters"
  )
  
  #todo: Don't hardcode carrying capacity
  cc <- case_when(
    pp=="High" ~ 800,
    pp=="Low" ~ 400,
    TRUE ~ 420
  )
  
  # Computing stopping time
  
  # add row indices to the original dataframe
  width <- 0.95
  
  plot_data <- dat %>%
    mutate(wAlb_adult = wAlb_adult/cc,
           wPip_adult = wPip_adult/cc) # Change everything to proportion land
  
  
  plot_data_medians <- plot_data %>%
    group_by(days, intervention, immigration) %>% 
    ggdist::curve_interval(.exclude=c("run", "controlProp", "params"), 
                           .width = width) %>% 
    ungroup()
  
  aux <- plot_data_medians %>%
    mutate(row_index = row_number())
  # aux <- dat %>%
  #   mutate(row_index = row_number())
  
  
  dat_stop_time <- aux %>%
    group_by(intervention, immigration) %>%
    mutate(change_flag = total_released_to_date != lag(total_released_to_date, default = first(total_released_to_date))) %>%
    filter(change_flag == TRUE) %>%
    summarise(antepenultimate_change_row = nth(row_index, -2)) %>%
    ungroup()
  
  # add a 'days' column based on 'antepenultimate_change_row'
  dat_stop_time <- dat_stop_time %>%
    left_join(aux %>% select(row_index, days), by = c("antepenultimate_change_row" = "row_index"))
  
  # view the result
  #print(dat_stop_time)
  dat_stop_time$days <- dat_stop_time$days + 7
  
  # Computing average stopping time
  dat_avg_stop_time <- dat_stop_time %>% 
    group_by(intervention, immigration) %>% 
    summarise(avg_stopping_time = round(mean(days))) %>% 
    mutate(after_6_months = avg_stopping_time + 180)
  
  # view the result
  #print(dat_avg_stop_time)
  
  
  
  # Stormclouds -------------------------------------------------------------
  
  ## 3 x 3
  
  ## Create the median curves (first 500 simulation days)
  
  width <- 0.95
  
  
  tmp_plot_data <- plot_data %>% 
    filter(days<=500)
  tmp_plot_data_medians <- plot_data_medians %>% 
    filter(days<=500)
  
  
  # Random sample of 10 actual simulation runs
  set.seed(1234)
  runs_subset <- sample(1:length(unique(dat_full$run)), 10, replace = F)
  
  
  
  p_storm_all <- 
    ggplot(data = subset(tmp_plot_data, run %in% runs_subset), aes(x=days)) +
    geom_line(aes(y = wAlb_adult, group = run), col = col.simruns, alpha=0.1) +
    geom_line(aes(y = wPip_adult, group = run), col = col.simruns, alpha=0.1) +
    geom_ribbon(data = tmp_plot_data_medians, 
                aes(ymin=wAlb_adult.lower, ymax=wAlb_adult.upper), 
                fill = col.medianwAlbAB, alpha=0.2) +
    geom_ribbon(data = tmp_plot_data_medians, 
                aes(ymin = wPip_adult.lower, ymax=wPip_adult.upper), 
                fill = col.medianwPip, alpha=0.2) +
    geom_line(data = tmp_plot_data_medians, 
              aes(y=wAlb_adult, col = "wAlbAB"), 
              key_glyph = "point") +  # Ensure labels here match in scale_color_manual
    geom_line(data = tmp_plot_data_medians, 
              aes(y=wPip_adult, col = "wPip"), 
              key_glyph = "point") +    # Ensure labels here match in scale_color_manual
    geom_hline(yintercept = 0.1, lty=2) +
    geom_hline(yintercept = 0.4, lty=2, colour = "gray60") +
    facet_grid(intervention ~ immigration,
               labeller = labeller(intervention = lab.intervention, immigration = lab.immigration)) +
    scale_x_continuous(breaks = seq(0,500,100), limits = c(0,500)) +
    scale_color_manual(
      values = c(col.medianwAlbAB, col.medianwPip), 
      labels = c(lab.wAlbAB, lab.wPip),
      guide = guide_legend(override.aes = list(linetype = 0, shape = 15, size = 5))  # Changes the legend to squares
    ) +
    theme_bw() +
    labs(y = "Proportion of steady state population", 
         x = "Days",
         col = "Population",
         fill = "") +
    theme(text = element_text(size = text_size),
          legend.position = "bottom",
          legend.title = element_text(size = 16),
          strip.text = element_text(size = 12))
  
  auxx <- plot_data_medians %>% 
    select(days, intervention, immigration, wAlb_adult, wPip_adult) 
  
  auxx %>% 
    group_by(intervention, immigration) %>% 
    filter(wAlb_adult < 0.1) %>% 
    summarise(day_first_below_0.1 = min(days), .groups = 'drop')
  
  ## Create the median curves for only 10 per week migration type intervention (after first 200 simulation days)
  
  tmp_plot_data <- plot_data %>% 
    filter(days>= 200)
  tmp_plot_data_medians <- plot_data_medians %>% 
    filter(days>= 200)
  
  tmp_plot_data <- tmp_plot_data %>% filter(intervention =="maintain")
  tmp_plot_data_medians <- tmp_plot_data_medians %>% filter(intervention =="maintain")
  
  
  p_storm_maintain <- ggplot(data = subset(tmp_plot_data, run %in% runs_subset), aes(x=days)) +
    geom_line(aes(y=wAlb_adult, group = run), col=col.simruns, alpha=0.1) +
    geom_line(aes(y=wPip_adult, group = run), col=col.simruns, alpha=0.1)+
    geom_ribbon(data = tmp_plot_data_medians, 
                aes(ymin=wAlb_adult.lower, ymax=wAlb_adult.upper), 
                fill = col.medianwAlbAB, alpha=0.3) +
    geom_ribbon(data = tmp_plot_data_medians, 
                aes(ymin=wPip_adult.lower, ymax=wPip_adult.upper), 
                fill = col.medianwPip, alpha=0.3) +
    geom_line(data = tmp_plot_data_medians, 
              aes(y=wAlb_adult, col = "wAlbAB"), 
              key_glyph = "point") +  # Ensure labels here match in scale_color_manual
    geom_line(data = tmp_plot_data_medians, 
              aes(y=wPip_adult, col = "wPip"), 
              key_glyph = "point") +    # Ensure labels here match in scale_color_manual
    geom_hline(yintercept = 0.1, lty=2) +
    geom_hline(yintercept = 0.4, lty=2, colour = "gray60") +
    facet_grid(intervention ~ immigration,
               labeller = labeller(intervention = lab.intervention, immigration = lab.immigration)) +
    scale_x_continuous(breaks = seq(200,910,200), limits = c(200,910)) +
    scale_color_manual(
      values = c(col.medianwAlbAB, col.medianwPip), 
      labels = c(lab.wAlbAB, lab.wPip),
      guide = guide_legend(override.aes = list(linetype = 0, shape = 15, size = 5))  # Changes the legend to squares
    ) +
    theme_bw() +
    labs(y = "Proportion of steady\nstate population", 
         x = "Days",
         col = "Population",
         fill = ""
         
    ) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 16),
          plot.title = element_text(size = 18),
          text = element_text(size = text_size))
  
  
  ## together
  p1 <- p_storm_all      + theme(strip.text = element_text(size = 17), axis.title.y = element_text(size = text_size), axis.title.x = element_text(size = text_size), legend.position = "")
  p2 <- p_storm_maintain + theme(strip.text = element_text(size = 17), axis.title.y = element_text(size = text_size), axis.title.x = element_text(size = text_size), legend.title = element_text(size = text_size), legend.text = element_text(size = text_size))
  
  p_storm <- p1 + p2 + plot_layout(ncol = 1, heights = c(3.25, 1.25), axis_titles = "collect") + plot_annotation(tag_levels = 'A') 
  
  ggsave(here::here(glue::glue("img/intervention_stormcloud_{pp}_nodecay.png")), 
         plot = p_storm,
         height = 15, width = 15, dpi = 300)
}
