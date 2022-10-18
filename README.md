# Anytime valid tests of conditional independence under model-X

This repository contains replication material and additional figures for the article

Grünwald, P., Henzi, A., and Lardy, T. (2022). Anytime Valid Tests of Conditional Independence Under Model-X. arXiv e-prints [arXiv:2209.12637](https://arxiv.org/abs/2209.12637)
 
 The file *simulation_functions.R* contains all functions used for the simulations in Section 4 of the article. The folder *simulation_1* contains the files to run the simulations for rejection rates under null and alternative hypothesis, and *simulation_2* the files for rejection rates under misspecification of the conditional distribution of X given Z. The code was run on a HPC cluster, and the files *simulation_1.R* and *simulation_2.R* generate the output for a single run of the simulations, with parameters depending on the variable *id*. Additional figures with different parameters for the simulations are contained in the folder *additional_figures*. The figures are generated with the file *simulation_plots.R*.
