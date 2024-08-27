This is a folder for simulations studies and data analysis for the manuscript "Semiparametric Change Point Model for Survival Outcomes in the Presence of
a U-Shaped Risk".

For simulation studies, relevant files are:
1. "simulation_one_rep.R", which is the code for running one replication of one scenario of simulation under a correctly specified model;
3. "simulation_correct_runbatch.R", which is provided as an illustration of the whole simulation study and can be used to run the whole simulation study under correctly specified models, but a parallel computing equivalence on a computing cluster is probably preferred;
4. "simulation_correct_summary.R", which summarizes simulation results obtained by running "simulation_one_rep.R" and "simulation_run_array.R", and formats them into tables shown in the manuscript;

File "simulation_correct_trueROC.R" should be run before files "simulation_correct_summary.R" and "simulation_misspecified_summary.R" can be run.

For real data analysis, the relevant file is "DREAM_analysis.R".
