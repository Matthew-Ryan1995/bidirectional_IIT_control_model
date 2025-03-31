# reid-wolbachia

This is the code to reproduce the results in **Modelling Aedes albopictus management, incorporating immigration and bidirectional Wolbachia interactions**.

# To replicate results

1. `code/00_runMe.R` - runs main intervention simulations
  a. `code/create_model_parameters.R` - Creates the model parameters to use in `00_runME.R` based on `code/parameters.xlsx`
  b. `code/create_scenario_parameters.R` - Creates the scenario parameters to run different scenarios in `code/00_runME.R`
  c. `code/collect_intervention_results.R` - Collates results produced by `00_runME.R`
  d. `code/run_simulations.R` - Wrapper to run the intervention simulations
  e. `code/model_functions.R` - All model functions to run experiments
2. `code/01_runME_cage_experiments.R` - Run in silico cage experiments
  a. `code/create_scenario_parameters_cage.R` - Creates the scenario parameters to run different scenarios in  `code/01_runME_cage_experiements.R`
  b. `code/collect_cage_results.R` - Collates results produced by `01_runME_cage_experiements.R`
  c. `code/run_simulations_cage.R` - Wrapper to run the cage simulations
  d. `code/model_functions.R` - All model functions to run experiments
3. Figures for the paper are produced using the scripts in `code/paper_figures/` 
4. Simulation results are saved in the `outputs` folder.
5. Collected data is saved in the `data` folder.
6. The scripts in `code/HPC_versions` are adapted versions of `code/00_runMe.R` and `code/01_runME_cage_experiments.R` to run on a SLURM HPC system.  The corresponding batch scripts are in the `batch_scripts` folder.
7. The folder `code/noDecay` contains the code to reproduce the results ith age-related CI decay turned off.


