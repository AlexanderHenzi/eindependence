# Anytime Valid Tests of Conditional Independence Under Model-X

This repository contains replication material and additional figures for the article

Grünwald, Peter, Alexander Henzi, and Tyron Lardy. "Anytime Valid Tests of Conditional Independence Under Model-X." arXiv preprint [arXiv:2209.12637](https://arxiv.org/abs/2209.12637) (2022).
 
Code for the data application is given in the folder *data_application*. The file *simulation_functions.R* contains all functions used for the simulations in Section 4 of the article. The folder *simulation_1* contains the files to run the simulations for rejection rates under null and alternative hypothesis, and *simulation_2* the files for rejection rates under misspecification of the conditional distribution of X given Z. The code was run on a HPC cluster, with reproducible seed for the simulations, and the files *simulation_1.R* and *simulation_2.R* generate the output for a single run of the simulations, with parameters depending on the variable *id* and other manually specified parameters as explained in comments in the files.

To reproduce the simulations under the alternative (files in *simulation_1/2*) (i) select the desired parameters in the file *simulation_1.R* (positive/negative correlations, dimension, penalization for RMLE) (ii) submit the batch job *sim_logistic_1_job.sh* on a HPC cluster with SLURM (adjust the email address, partition name etc. in the shell scripts and upload all required files, including *simulation_functions.R*); this will automatically collect the results in a file called *simulation_1/2.rda*, which (iv) is read in *simulation_1_plots_table.R* to reproduce the results. The same workflow applies to *simulation_2* and the file *simulation_2_plots.R* for generating the plots concerning robustness under violations of the Model-X assumptions.

The table below gives an overview of the parameters in the files for reproducing all figures and tables in the article and supplement. Figures for other parameter combinations are contained in the folder *additional_figures*.

| Objects            | Files                                | Parameters                                                        | 
| ------------------ | ------------------------------------ | --------------------------------------------------------          |
| Table 1            | simulation_1.R                       | *q = 4*, *correlation = pos*, *not_penalize_rmle = TRUE*          |
|                    | simulation_1_plots_table.R           | *eps = 0.05*                                                      |
| Figures 1, 2, S4   | simulation_1.R                       | *q = 4*, *correlation = pos*, *not_penalize_rmle = TRUE*          |
|                    | simulation_1_plots_table.R           | *eps = 0.05*                                                      |
| Figure 3           | data_application.R                   | none                                                              |
| Figure S1          | simulation_1.R                       | *q = 4*, *correlation = pos*, *not_penalize_rmle = TRUE*          |
|                    | simulation_1_plots_table.R           | *eps = 0*                                                         |
| Figure S2          | simulation_2_tests.R                 | none                                                              |
| Figure S3          | simulation_2.R                       | *q = 4*, *correlation = pos*                                      |
|                    | simulation_2_plots.R                 | none                                                              |
| Figures S5, S7     | simulation_1.R                       | *q = 8*, *correlation = pos*, *not_penalize_rmle = TRUE*          |
|                    | simulation_1_plots_table.R           | *eps = 0.05*                                                      |
| Figures S6, S8     | simulation_1.R                       | *q = 4*, *correlation = negative_cor*, *not_penalize_rmle = TRUE* |
|                    | simulation_1_plots_table.R           | *eps = 0.05*                                                      |
| Figure S9          | simulation_2.R                       | *q = 8*, *correlation = pos*                                      |
|                    | simulation_2_plots.R                 | none                                                              |
| Figure S10         | simulation_2.R                       | *q = 4*, *correlation = negative_cor*                             |
|                    | simulation_2_plots.R                 | none                                                              |