This is a folder for simulations studies and data analysis for the manuscript "Semiparametric Change Point Model for Survival Outcomes in the Presence of a U-Shaped Risk".

For simulation studies, relevant files are:
1. "simulation_one_rep.R", which is the code for running one replication of one scenario of simulation;
2. "run_array_simulation.R", which is provided as an illustration of the whole simulation study and can be used to run the whole simulation study, but a parallel computing equivalence on a computing cluster is probably preferred;
3. "C_S_test.R", which is the code for calculating estimated C-index and survival probability on testing dataset for one scenario of simulation;
4. "run_array_test_data.R", which is provided as an illustration of the whole batch running process to calculate the estimated C-index and survival probability on the testing dataset for every scenario;
5. "simulation_summary.R", which summarizes simulation results obtained by running "simulation_one_rep.R" and "simulation_run_array.R", and formats them into tables shown in the manuscript;


For real data analysis, the relevant files are:
1. "real data_lilly.R", which is corresponding to the real data analysis part in the manuscript;
2. "bootstrap_realdata_lilly.R", which is used to realize the bootstrap for estimates in "real data_lilly.R".
