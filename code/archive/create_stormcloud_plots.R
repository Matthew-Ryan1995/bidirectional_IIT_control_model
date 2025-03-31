## Matt Ryan
## 14-11-2023
## Create stormcloud plots
## Requires: tidyverse, ggdist, here


# libraries ---------------------------------------------------------------

pacman::p_load(tidyverse)


# plotting parameters -----------------------------------------------------

height <- 7
dpi <- 200
fontsize <- 30


# Load data ---------------------------------------------------------------

df <- read_rds(here::here("data/all_scenarios_all_runs.Rds"))

# remove to eliminated cases
df <- df %>% 
  filter(wAlbAB_adults > 0 & wPip_adults > 0) %>% 
  select(-total_released_to_date)

median_curves <- df %>% 
  group_by(OF, FCP, days) %>% 
  ggdist::median_qi(.exclude=c("run")) %>% 
  ungroup()

fcp_list <- unique(df$FCP)
of_list <-  unique(df$OF)
for(fcp in fcp_list){
  for(of in of_list){
    df_tmp <- df %>% 
      filter(OF==of, FCP==fcp)
    median_tmp <- median_curves %>% 
      filter(OF==of, FCP==fcp)
    
    supression_date <- median_tmp %>% 
      filter(wAlbAB_prop <0.1) %>% 
      slice(1) %>% 
      pull(days)
    
    p <- median_tmp %>% 
      ggplot(aes(x=days)) + # Everything has common x-axis
      geom_line(data = df_tmp, # Plot the individual runs for each mosquito type
                aes(y = wAlbAB_prop, group = run), 
                colour = "black", alpha = 0.1) +
      geom_line(data = df_tmp , 
                aes(y = wPip_prop, group = run), 
                colour = "black", alpha = 0.1) +
      geom_ribbon(aes(ymin = wAlbAB_prop.lower, ymax= wAlbAB_prop.upper), # Plot confidence bands and median for each mosquito
                  alpha = 0.3, fill = "blue") +
      geom_line(aes(y = wAlbAB_prop), colour = "blue") +
      geom_ribbon(aes(ymin = wPip_prop.lower, ymax= wPip_prop.upper), 
                  alpha = 0.3, fill = "red") +
      geom_line(aes(y = wPip_prop), colour = "red") +
      geom_vline(xintercept = supression_date, colour="grey70", lty=2) + # add median supression date
      theme_bw() +
      theme(text = element_text(size=fontsize)) +
      labs(x = "Time (days)",
           y = "Proportion of initial population",
           # Caption to be removed
           #caption = "red = median wPip, blue = median wAlbAB, 95% CI band\nDotted line = median supression date"
           )
    
    ggsave(here::here(glue::glue("img/stormcloud_OF_{of}_FCP_{fcp}.png")), plot = p, 
           height = height, width = 1.618*height, dpi = dpi)
  }
}
