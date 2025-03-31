## Matt Ryan
## 14-11-2023
## Create cost and time to suppression plots
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

# Time to suppression plot ------------------------------------------------

p <- suppression_data %>% 
  ggplot(aes(x = days, y = OF, fill = FCP)) +
  geom_boxplot() +
  labs(fill = "Female contamination probability",
       x = "Days until supression",
       y = "Overflooding ratio") +
  harrypotter::scale_fill_hp_d("ravenclaw") +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size=fontsize))

ggsave(here::here(glue::glue("img/time_to_suppression.png")), plot = p, 
       height = height, width = 1.618*height, dpi = dpi)

# Cost of suppression -----------------------------------------------------

p <- suppression_data %>% 
  ggplot(aes(x = total_released_to_date, y = OF, fill = FCP)) +
  geom_boxplot() +
  labs(fill = "Female contamination probality",
       x = "Total adult mosquitoes released",
       y = "Overflooding ratio") +
  harrypotter::scale_fill_hp_d("ravenclaw") +
  theme_bw() +
  theme(legend.position = "bottom", text=element_text(size=fontsize))


ggsave(here::here(glue::glue("img/cost_of_suppression.png")), plot = p, 
       height = height, width = 1.618*height, dpi = dpi)

