
# Table 1 - Main manuscript and Tables S03 and S04 in supp.
# Percentage of sims with management success for wAlbAB and wPip populations
# under different scenarios for immigration rates and intervention strategies. 
# Management success is defined as keeping the wAlbAB population below 10% 
# of the initial population size, and wPip below 40% (unstable eq. threshold). 
# The data is shown both within the stopping time and six months after. 



# packages ----------------------------------------------------------------

pacman::p_load(tidyverse, dplyr, knitr, kableExtra)


# data and parameters -----------------------------------------------------

dat_full <- read_rds(here::here("data/intervention_results_nodecay_Expected.Rds"))

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


# Tables -------------------------------------------------------------------

# Pre-allocate list to store LaTeX tables
latex_code_list <- list()

for(pp in unique(dat_full$params)){
  
  #pp <- "Expected"
  
  dat <- dat_full %>% filter(params==pp)

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
  aux <- dat %>%
    mutate(row_index = row_number()) # add row indices to the original dataframe
  
  # Compute the row index of the second-to-last time 'total_released_to_date' changed within each group
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
  
 
  # ============================================================================
  # Function to compute success at different times (stop time and after 6 months) 
  # for both wAlbAB and wPip
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
  # ============================================================================
  
  
  ## Establishment proportions -----------------------------------------------
  dat_wAlbAB_stop <- calculate_establishment(tmp, "wAlbAB", 0.1, "endtime")
  dat_wPip_stop <- calculate_establishment(tmp, "wPip", 0.4, "endtime")
  dat_wAlbAB_6months <- calculate_establishment(tmp, "wAlbAB", 0.1, "after_6_months")
  dat_wPip_6months <- calculate_establishment(tmp, "wPip", 0.4, "after_6_months")
  
  # combine data sets
  dat_both <- bind_rows(dat_wAlbAB_stop, dat_wPip_stop, dat_wAlbAB_6months, dat_wPip_6months)
  
  # ensure 'time' is a factor 
  dat_both$time <- factor(dat_both$time, levels = c("endtime", "after_6_months"))

  
  # Arrange and reshape the data before generating the table
  dat_table <- dat_both %>%
    mutate(
      immigration = factor(immigration, levels = c("0 per week", "2 per week", "10 per week")),
      time = factor(time, levels = c("endtime", "after_6_months"))
      ) %>%
    arrange(immigration, time) %>%
    select(intervention, immigration, name, n_success, time) %>%
    pivot_wider(names_from = c(intervention, name), values_from = n_success) %>%
    select(immigration, time, 
           naive_wAlbAB, naive_wPip, 
           `complete stop_wAlbAB`, `complete stop_wPip`, 
           maintain_wAlbAB, maintain_wPip) %>% 
    mutate(across(where(is.numeric), ~ round(.x, 1)))

  # Start LaTeX table
  latex_code <- paste0("\\begin{table}[!h]
      \\caption{Percentage of simulations with management success for \\walbab\\ and \\wpip\\ populations under different scenarios (", pp, ") for immigration rates and intervention strategies. Management success is defined as keeping the \\walbab\\ population below $10\\%$ of the initial population size, and \\wpip\\ below $40\\%$ (unstable equilibrium threshold). The data is shown both within the stopping time and six months after. \\textbf{The model parameters are the ", pp, " values defined in Table \\ref{tab:parameters}.}}
      \\label{tab:intervention_establishment_", pp, "}
      \\centering
      \\begin{tabular}[t]{lrrrrrr}
      \\toprule
      Scenario &  \\multicolumn{2}{c}{Na\\\"ive} & \\multicolumn{2}{c}{Complete stop}  & \\multicolumn{2}{c}{Maintain}\\\\
      \\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}
      & \\walbab & \\wpip & \\walbab & \\wpip & \\walbab & \\wpip \\\\
      \\midrule")
  
  # Loop through each scenario to add data rows
  scenarios <- unique(dat_table$immigration)
  
  for (scenario in scenarios) {
    latex_code <- paste0(latex_code, "
      \\addlinespace[0.3em]
      \\multicolumn{7}{l}{\\textbf{", scenario, "}}\\\\")
    
    subset_data <- dat_table %>% filter(immigration == scenario)
    
    first_row <- subset_data[1, ]
    second_row <- subset_data[2, ]
    
    latex_code <- paste0(
      latex_code,
      "
      \\rowcolor{gray!6} {\\hspace{1em}Within stopping time} & ",
      first_row$naive_wAlbAB, " & ", first_row$naive_wPip, " & ",
      first_row$`complete stop_wAlbAB`, " & ", first_row$`complete stop_wPip`, " & ",
      first_row$maintain_wAlbAB, " & ", first_row$maintain_wPip, "\\\\",
      
      "
      \\hspace{1em}After six months & ",
      second_row$naive_wAlbAB, " & ", second_row$naive_wPip, " & ",
      second_row$`complete stop_wAlbAB`, " & ", second_row$`complete stop_wPip`, " & ",
      second_row$maintain_wAlbAB, " & ", second_row$maintain_wPip, "\\\\"
    )
  }
  
  # Finalize the LaTeX code
  latex_code <- paste0(latex_code, "
        \\bottomrule
      \\end{tabular}
    \\end{table}")
  
  # Store LaTeX code in list
  latex_code_list[[pp]] <- latex_code
  
}

# Print LaTeX code for each scenario - It needs "small" adjustments...
for (pp in names(latex_code_list)) {
  cat("\n\n%%% LaTeX Table for", pp, "Scenario %%%\n")
  cat(latex_code_list[[pp]], "\n")
}



