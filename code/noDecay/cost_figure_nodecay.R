# Matt Ryan
# 09/09/2024
## Create cost results from interventions
## Requires: pacman, here, tidyverse, patchwork


# Libraries ---------------------------------------------------------------
pacman::p_load(tidyverse, patchwork)


# Params ------------------------------------------------------------------
height <- 7
width <- 14
dpi <- 600

text_size <- 20
x_axis_text <- 8

runtime <- 2*365




pp <- "Expected"

dat_full <- read_rds(here::here(glue::glue("data/intervention_results_nodecay_Expected.Rds")))

dat_full <- dat_full %>%
  mutate(immigration = case_when(
    immigration=="0" ~ "0 per week",
    immigration=="0.286" ~ "2 per week",
    TRUE ~ "10 per week"
  ),
  immigration = factor(immigration,
                       levels = c("0 per week", "2 per week", "10 per week"))) %>%
  rename(Immigration=immigration)

dat <- dat_full %>% 
  # filter(params==pp) %>% 
  filter(days<=runtime)

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

cost_data <-  dat %>% 
  mutate(day_groups = floor(days/100)) %>% 
  group_by(run, intervention, Immigration, day_groups) %>% 
  slice_max(total_released_to_date) %>% 
  slice(1) %>% 
  select(days, total_released_to_date, `run`, `intervention`, `Immigration`, `day_groups`) %>% 
  ungroup(`day_groups`) %>% 
  mutate(lag=lag(total_released_to_date),
         lag=ifelse(is.na(lag), 0, lag),
         diff = (total_released_to_date-lag)/1000) %>% 
  filter(diff>=0, day_groups<7)

cost_avg <- dat %>% 
  mutate(day_groups = floor(days/100)) %>% 
  group_by(run, intervention, Immigration, day_groups) %>% 
  slice_max(total_released_to_date) %>% 
  slice(1) %>% 
  select(days, total_released_to_date, `run`, `intervention`, `Immigration`, `day_groups`) %>% 
  filter(total_released_to_date > 0) %>% 
  mutate(total_released_to_date = total_released_to_date/1000) %>% 
  group_by(intervention, Immigration, day_groups) %>% 
  summarise(mean = mean(total_released_to_date),
            sd = sd(total_released_to_date),
            .groups = "drop") %>% 
  filter(day_groups<7)

p1 <- cost_data %>% 
  ggplot(aes(x=as.factor(day_groups),
             fill = intervention)) +
  # geom_jitter(position=position_dodge(width=0.75),
  #             aes(colour=intervention,
  #                 y = (diff)),
  #             show.legend = FALSE,
  #             size = 1)+
  geom_boxplot(alpha=1,#0.7,
               aes(y = (diff)), 
               show.legend = FALSE) +
  facet_wrap(~Immigration, nrow=1, labeller = "label_both") +
  labs(x=NULL, y="Ongoing release\n('000s of mosquitoes)",
       fill = "Intervention strategy",
       colour = "Intervention strategy") +
  harrypotter::scale_fill_hp_d("ravenclaw", direction = -1, labels = c("Naïve","Complete stop","Maintain")) +
  harrypotter::scale_colour_hp_d("ravenclaw", direction = -1, labels = c("Naïve","Complete stop","Maintain")) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text=element_text(size=text_size))

if(pp=="High"){
  p1 <- cost_data %>% 
    mutate(is_low = diff<15) %>% 
    ggplot(aes(x=as.factor(day_groups),
               fill = intervention)) +
    geom_jitter(position=position_dodge(width=0.75),
                aes(colour=intervention,
                    y = (diff)),
                show.legend = FALSE,
                size = 1)+
    geom_boxplot(alpha=0.7,
                 aes(y = (diff)), 
                 show.legend = FALSE) +
    scale_y_continuous(minor_breaks = NULL) +
    facet_grid(is_low~Immigration, scales="free_y", labeller = "label_both") +
    labs(x=NULL, y="Ongoing release\n('000s of mosquitoes)",
         fill = "Intervention strategy",
         colour = "Intervention strategy") +
    harrypotter::scale_fill_hp_d("ravenclaw", direction = -1, labels = c("Naïve","Complete stop","Maintain")) +
    harrypotter::scale_colour_hp_d("ravenclaw", direction = -1, labels = c("Naïve","Complete stop","Maintain")) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          text=element_text(size=text_size),
          strip.background.y = element_blank(),
          strip.text.y = element_blank())
}

p2 <- cost_avg %>% 
  ggplot() +
  geom_ribbon(data=cost_avg, 
              aes(x=day_groups+1, 
                  ymin=mean - 2*sd,
                  ymax=mean + 2*sd,
                  fill=intervention),
              # position=position_dodge(width=0.8), 
              alpha=0.1, 
              show.legend = TRUE) +
  geom_line(data=cost_avg, 
            aes(x=day_groups+1, 
                y=mean,
                colour=intervention),
            # position=position_dodge(width=0.8), 
            show.legend = FALSE)+
  facet_wrap(~Immigration, nrow=1, labeller = "label_both") +
  labs(x="Time (days)",
       y = "Cumulative release\n('000s of mosquitoes)",
       fill = "Intervention strategy",
       colour = "Intervention strategy") +
  scale_x_continuous(breaks = 1:7,
                     labels = c("0-99",
                                "100-199",
                                "200-299",
                                "300-399",
                                "400-499",
                                "500-599",
                                "600-699"),
                     minor_breaks = NULL) +
  harrypotter::scale_fill_hp_d("ravenclaw", direction = -1, labels = c("Naïve","Complete stop","Maintain")) +
  harrypotter::scale_colour_hp_d("ravenclaw", direction = -1, labels = c("Naïve","Complete stop","Maintain")) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +
  theme(legend.position = "bottom",
        text=element_text(size=text_size),
        axis.text.x=element_text(size=x_axis_text))

p <- (p1/plot_spacer()/p2) + 
  plot_layout(guides="collect", heights = c(4,0.3,4)) & 
  theme(legend.position = "bottom")

ggsave(here::here(glue::glue("img/intervention_cost_by_dayblock_nodecay_{pp}.png")),
       plot=p, 
       height = height, 
       width=width, 
       dpi=dpi)
