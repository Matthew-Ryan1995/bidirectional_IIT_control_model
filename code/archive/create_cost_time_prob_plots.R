## Matt Ryan
## 14-11-2023
## Create cost and time to probability plots
## Requires: tidyverse, here, harrypotter


# libraries ---------------------------------------------------------------

pacman::p_load(tidyverse)


# plotting parameters -----------------------------------------------------

height <- 7
dpi <- 200
fontsize <- 30


# Data --------------------------------------------------------------------

suppression_data <- read_rds(here::here("data/all_scenarios_supressed.Rds"))

suppression_data <- suppression_data %>% 
  mutate(OF = str_c("1:", OF),
         OF = factor(OF, levels = c("1:1", "1:5", "1:15")))


# Time plot ---------------------------------------------------------------

p <- suppression_data %>% 
  ggplot(aes(x = days, colour = OF, lty = FCP)) +
  stat_ecdf(geom = "line") +
  labs(x = "Length of intervention (days)",
       y = "Cum. prop'n suppressed",
       lty = "Female contamination probability",
       colour = "Overflooding ratio") +
  theme_bw() +
  guides(colour=guide_legend(ncol=2)) +
  harrypotter::scale_color_hp_d("ravenclaw") +
  theme(legend.position = "bottom",
        legend.direction = "vertical", text=element_text(size=fontsize))

ggsave(here::here(glue::glue("img/time_probability.png")), plot = p, 
       height = height, width = 1.618*height, dpi = dpi)

# Cost plot ---------------------------------------------------------------

p <- suppression_data %>% 
  ggplot(aes(x = total_released_to_date, colour = OF, lty = FCP)) +
  stat_ecdf(geom = "line") +
  labs(x = "Number of adults released",
       y = "Cum. prop'n suppressed",
       lty = "Female contamination probability",
       colour = "Overflooding ratio") +
  theme_bw() +
  guides(colour=guide_legend(ncol=2)) +
  harrypotter::scale_color_hp_d("ravenclaw") +
  theme(legend.position = "bottom",
        legend.direction = "vertical", text=element_text(size=fontsize))

ggsave(here::here(glue::glue("img/cost_probability.png")), plot = p, 
       height = height, width = 1.618*height, dpi = dpi)
