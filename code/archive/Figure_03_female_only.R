# Matt Ryan
# 09/02/2024
## Create figures of results from intervention simulations
## Requires: pacman, here, tidyverse



# Figure 3
# Stormclouds + Zoom


# packages ----------------------------------------------------------------

pacman::p_load(tidyverse, patchwork)

# data and parameters -----------------------------------------------------

dat_full <- read_rds(here::here("data/intervention_results_full.Rds"))

dat_full <- dat_full %>% 
  mutate(immigration = case_when(
    immigration=="0" ~ "0 per week",
    immigration=="0.286" ~ "2 per week",
    TRUE ~ "10 per week"
  ),
  immigration = factor(immigration, levels = c("0 per week", "2 per week", "10 per week")))


# Find reduced dataset where sims stopped
dat_endtime_full <- dat_full %>% 
  filter(wPip >0 | wAlb > 0) %>% 
  group_by(run, intervention, immigration, params) %>% 
  slice_max(days) %>% 
  ungroup() 


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
lab.wPip <- expression(paste(italic(w),Pip))

lab.intervention <- c(`naive` = "Intervention: naïve",
                      `complete stop` = "Intervention: complete stop",
                      `maintain` = "Intervention: maintain")

lab.immigration <- c(`0 per week` = "Immigration: 0 per week",
                     `2 per week` = "Immigration: 2 per week",
                     `10 per week` = "Immigration: 10 per week")

# plots -------------------------------------------------------------------

for(pp in unique(dat_full$params)){
  
  #pp <- "Expected" 
  
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
  # cc <- case_when(
  #   pp=="High" ~ 800,
  #   pp=="Low" ~ 400,
  #   TRUE ~ 420
  #   )
  cc <- dat %>% 
    filter(days==1) %>% 
    slice(1) %>% 
    pull(block01_wAlbAB_f)
  
  
  # Computing stopping time
  
  # add row indices to the original dataframe
  aux <- dat %>%
    mutate(row_index = row_number())
  
  # compute the row index of the second-to-last time 'total_released_to_date' changed within each group
  dat_stop_time <- aux %>%
    group_by(intervention, immigration, run) %>%
    mutate(change_flag = total_released_to_date != lag(total_released_to_date, default = first(total_released_to_date))) %>%
    filter(change_flag == TRUE) %>%
    summarise(antepenultimate_change_row = nth(row_index, -2)) %>%
    ungroup()
  
  # add a 'days' column based on 'antepenultimate_change_row'
  dat_stop_time <- dat_stop_time %>%
    left_join(aux %>% select(row_index, days), by = c("antepenultimate_change_row" = "row_index"))
  
  # view the result
  #print(dat_stop_time)
  
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
  
  plot_data <- dat %>%
    filter(days <= 500) %>%
    filter(wPip > 0 | wAlb > 0) %>%
    mutate(wAlb_adult = block01_wAlbAB_f/cc,
           wPip_adult = block01_wPip_f/cc) # Change everything to proportion land
  
  plot_data_medians <- plot_data %>%
    filter(days <= 500) %>% 
    group_by(days, intervention, immigration) %>% 
    ggdist::median_qi(.exclude=c("run", "controlProp", "params")) %>% 
    ungroup()
  
  # Random sample of 10 actual simulation runs
  set.seed(1234)
  runs_subset <- sample(1:length(unique(dat_full$run)), 10, replace = F)
  col.avg_stop_time <- "forestgreen"
  col.after_6_months <- "orange2"
  p_storm_all <- ggplot(data = subset(plot_data, run %in% runs_subset), aes(x=days)) +
    geom_line(aes(y = wAlb_adult, group = run), col = col.simruns, alpha=0.1) +
    geom_line(aes(y = wPip_adult, group = run), col = col.simruns, alpha=0.1) +
    geom_ribbon(data = plot_data_medians, aes(ymin=wAlb_adult.lower, ymax=wAlb_adult.upper), fill = col.medianwAlbAB, alpha=0.2) +
    geom_ribbon(data = plot_data_medians, aes(ymin = wPip_adult.lower, ymax=wPip_adult.upper), fill = col.medianwPip, alpha=0.2) +
    geom_line(data = plot_data_medians, aes(y=wAlb_adult, col = col.medianwAlbAB)) +
    geom_line(data = plot_data_medians, aes(y=wPip_adult, col = col.medianwPip)) +
    geom_hline(yintercept = 0.1, lty=2) +
    geom_hline(yintercept = 0.4, lty=2, colour = "gray70") +
    geom_vline(data = dat_avg_stop_time, aes(xintercept = avg_stopping_time), lty=1, colour = "forestgreen") +
    geom_vline(data = dat_avg_stop_time, aes(xintercept = after_6_months), lty=1, colour = "orange2") +
    facet_grid(intervention ~ immigration,
               labeller = labeller(intervention = lab.intervention, immigration = lab.immigration)
               ) +
    scale_x_continuous(breaks = seq(0,500,100), limits = c(0,500)) +
    scale_color_manual(values = c(col.medianwAlbAB, col.medianwPip, col.avg_stop_time, col.after_6_months), 
                       labels = c(lab.wAlbAB, lab.wPip)) +
    theme_bw() +
    labs(y = "Proportion of steady state population", 
         x = "Days",
         col = "Population",
         fill = ""
         ) +
    theme(text = element_text(size = text_size),
          legend.position = "bottom",
          legend.title = element_text(size = 16),
          strip.text = element_text(size = 12))
  #p_storm_all 
  
  
  ## Create the median curves for only 10 per week migration type intervention (after first 200 simulation days)
  plot_data <- dat %>% 
    filter(days >= 200) %>% 
    filter(wPip > 0 | wAlb > 0) %>% 
    mutate(wAlb_adult = block01_wAlbAB_f/cc,
           wPip_adult = block01_wPip_f/cc) # Change everything to proportion land
  
  plot_data_medians <- plot_data %>% 
    filter(days >= 200) %>% 
    group_by(days, intervention, immigration) %>% 
    ggdist::median_qi(.exclude=c("run", "controlProp", "params")) %>% 
    ungroup()

  plot_data <- plot_data %>% filter(intervention =="maintain")
  plot_data_medians <- plot_data_medians %>% filter(intervention =="maintain")
  
  
  p_storm_maintain <- ggplot(data = subset(plot_data, run %in% runs_subset), aes(x=days)) +
    geom_line(aes(y=wAlb_adult, group = run), col=col.simruns, alpha=0.1) +
    geom_line(aes(y=wPip_adult, group = run), col=col.simruns, alpha=0.1)+
    geom_ribbon(data = plot_data_medians, aes(ymin=wAlb_adult.lower, ymax=wAlb_adult.upper), fill = col.medianwAlbAB, alpha=0.3) +
    geom_ribbon(data = plot_data_medians, aes(ymin=wPip_adult.lower, ymax=wPip_adult.upper), fill = col.medianwPip, alpha=0.3) +
    geom_line(data = plot_data_medians, aes(y=wAlb_adult, col = col.medianwAlbAB)) +
    geom_line(data = plot_data_medians, aes(y=wPip_adult, col = col.medianwPip)) +
    geom_hline(yintercept = 0.1, lty=2) +
    geom_hline(yintercept = 0.4, lty=2, colour = "gray70") +
    geom_vline(data = subset(dat_avg_stop_time, intervention == "maintain" & immigration != "0 per week"), aes(xintercept = avg_stopping_time), lty=1, colour = "forestgreen") +
    geom_vline(data = subset(dat_avg_stop_time, intervention == "maintain" & immigration != "0 per week"), aes(xintercept = after_6_months), lty=1, colour = "orange2") +
    facet_wrap(~immigration, labeller = labeller(immigration = lab.immigration)) +
    scale_x_continuous(breaks = seq(200,910,200), limits = c(200,910)) +
    scale_color_manual(values = c(col.medianwAlbAB, col.medianwPip), 
                       labels = c(lab.wAlbAB, lab.wPip)) +
    theme_bw() +
    labs(y = "Proportion of steady state population", 
         x = "Days",
         col = "Population",
         fill = "",
         title = "Intervention: maintain"
         ) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 16),
          plot.title = element_text(size = 18),
          text = element_text(size = text_size))

    
    ## together
    p1 <- p_storm_all + theme(legend.position = "")
    p2 <- p_storm_maintain + theme(axis.title.y = element_text(size = 14.5))
    
    p_storm <- p1 + p2 + plot_layout(ncol = 1, heights = c(3, 1.25),axis_titles = "collect") + plot_annotation(tag_levels = 'A') 
    
    ggsave(here::here(glue::glue("img/femaleOnly_interventions_stormcloud_{pp}_MM.png")), 
           plot = p_storm,
           height = 15, width = 15, dpi = 300)

}




############# 
# # To play with the computation of the average stopping time
# dat_full <- read_rds(here::here("data/intervention_results_full.Rds"))
# 
# dat_full <- dat_full %>% 
#   mutate(immigration = case_when(
#     immigration=="0" ~ "0 per week",
#     immigration=="0.286" ~ "2 per week",
#     TRUE ~ "10 per week"
#   ),
#   immigration = factor(immigration, levels = c("0 per week", "2 per week", "10 per week"))) %>% 
#   filter(params=="Expected" )
# 
# names(dat_full)
# 





# 
# # Example dataframe
# df <- data.frame(
#   col = c(1, 1, 2, 2, 3, 3, 3, 4)
# )
# 
# # Find the indices where the column values change
# 
# 
# # Add the first row as a starting point (if needed)
# change_indices <- c(1, change_indices)
# 
# # Show the result
# change_indices

