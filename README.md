# Anytime Valid Tests of Conditional Independence Under Model-X

This repository contains replication material and additional figures for the article

Peter Grünwald, Alexander Henzi, and Tyron Lardy. "Anytime Valid Tests of Conditional Independence Under Model-X." arXiv preprint [arXiv:2209.12637](https://arxiv.org/abs/2209.12637) (2022).
 
Code for the data application is given in the folder `data_application`. The file `simulation_functions.R` contains all functions used for the simulations in Section 4 of the article. The folder `simulation_1` contains the files to run the simulations for rejection rates under null and alternative hypothesis, and `simulation_2` the files for rejection rates under misspecification of the conditional distribution of X given Z. The code was run on a HPC cluster, with reproducible seed for the simulations, and the files `simulation_1.R` and `simulation_2.R` generate the output for a single run of the simulations, with parameters depending on the variable `id` and other parameters that can be passed as command line arguments as explained in comments in the files.

To reproduce all the simulations in the article and more (i) clone this repository on an HPC cluster (ii) adjust the email address, partition name, desired R module etc. in the files `sim_logistic_{1,2}_{array,collect,plot}.sh` in the directories `simulation_{1,2}/`
(iii) make sure that the following R packages are installed for the R module previously specified: ggpubr, ggthemes, ldbounds, MASS and tidyverse
(iv) from the main directory, submit the batch job `figure_reproduction.sh`, i.e. `sbatch figure_reproduction.sh`
(v) when the job finishes, the tables and figures can be found in the main directory. The default parameters for the figures in the article are `correlation = "pos"`, `not_penalize_rmle = FALSE`, and `q = 4`. For other values of these parameters (`correlation = "negative_cor"`, `not_penalize_rmle = TRUE`, or `q = 8`), the corresponding parameter is appended to the file name. For example, Table 1 corresponds to `table_1_eps05.txt`, Figure 1 to `simulation_alternative_eps05.pdf`, Figure 2 to `simulation_group_sequential_eps05.pdf`, etc.

To further experiment with simulations under the null and alternative, parameters can be manually specified in the bash file `simulation_1/sim_logistic_1_job.sh` by adding e.g. `--export=q=8,corr="pos",pen="TRUE"`.
The workflow for `simulation_2` to generate plots about robustness under violations of the model-X assumption is similar (with only parameters `q` and `corr`).

The table below gives an overview of the parameters in the files for reproducing all figures and tables in the article and supplement. Figures for other parameter combinations are contained in the directory `additional_figures`.

| Objects            | Files                                | Parameters                                                           | 
| ------------------ | ------------------------------------ | --------------------------------------------------------             |
| Table 1            | simulation_1.R                       | *q = 4*, *correlation = "pos"*, *not_penalize_rmle = FALSE*          |
|                    | simulation_1_plots_table.R           | *eps = 0.05*                                                         |
| Figures 1, 2, S4   | simulation_1.R                       | *q = 4*, *correlation = "pos"*, *not_penalize_rmle = FALSE*          |
|                    | simulation_1_plots_table.R           | *eps = 0.05*                                                         |
| Figure 3           | data_application.R                   | none                                                                 |
| Figure S1          | simulation_1.R                       | *q = 4*, *correlation = "pos"*, *not_penalize_rmle = FALSE*          |
|                    | simulation_1_plots_table.R           | *eps = 0*                                                            |
| Figure S2          | simulation_2_tests.R                 | none                                                                 |
| Figure S3          | simulation_2.R                       | *q = 4*, *correlation = "pos"*                                       |
|                    | simulation_2_plots.R                 | none                                                                 |
| Figures S5, S7     | simulation_1.R                       | *q = 8*, *correlation = "pos"*, *not_penalize_rmle = FALSE*          |
|                    | simulation_1_plots_table.R           | *eps = 0.05*                                                         |
| Figures S6, S8     | simulation_1.R                       | *q = 4*, *correlation = "negative_cor"*, *not_penalize_rmle = FALSE* |
|                    | simulation_1_plots_table.R           | *eps = 0.05*                                                         |
| Figure S9          | simulation_2.R                       | *q = 8*, *correlation = "pos"*                                       |
|                    | simulation_2_plots.R                 | none                                                                 |
| Figure S10         | simulation_2.R                       | *q = 4*, *correlation = "negative_cor"*                              |
|                    | simulation_2_plots.R                 | none                                                                 |