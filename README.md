This is a folder for simulations studies and data analysis for the manuscript "Semiparametric Change Point Model for Survival Outcomes in the Presence of
a U-Shaped Risk".

For simulation studies, relevant files are:
1. "simulation_one_rep.R", which is the code for running one replication of one scenario of simulation;
2. "run_array_simulation.R", which is provided as an illustration of the whole simulation study and can be used to run the whole simulation study, but a parallel computing equivalence on a computing cluster is probably preferred;
3. "C_S_test.R", which is the code for calculating estimated C-index and survival probability on testing dataset for one scenario of simulation;
4. "run_array_test_data.R", which is provided as an illustration of the whole batch running process to calculate the estimated C-index and survival probability on the testing dataset for every scenario;
5. "simulation_summary.R", which summarizes simulation results obtained by running "simulation_one_rep.R" and "simulation_run_array.R", and formats them into tables shown in the manuscript;

File "simulation_correct_trueROC.R" should be run before files "simulation_correct_summary.R" and "simulation_misspecified_summary.R" can be run.

For real data analysis, the relevant file is "DREAM_analysis.R".
