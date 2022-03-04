# Pollen Competition Modeling

### Code

Netlogo Folder - In this folder you'll find the one .nlogo file that generates all of our synthetic data. The Pollen_Competiton.nlogo file contains BehaviorSpace experiments to:

1. Generate LHS data from single accession simulations
2. Generate LHS data from double accession simulations
3. Generate synthetic experimental data for ler dominance based on the best parameter combinations

Germination Fitting - In this folder you'll find the .m files to fit the germination rate curves and caclulate AICc scores.

Main Directory - Here you'll find 
1. Gen_LHS.m which is how we generated our LHS (note to get the exact LHS we had reference the available data section below)
2. A Jupyter Notebook to identify which parameters are best
3. A Jupyter Notebook used to generate the figures in our report

### Available data: 

Repo 1: https://data.mendeley.com/v1/datasets/9xd5cc2nbv/draft?a=1a9cecda-1f44-4a19-91fb-663bc86c3a61 

* best_distr_params.csv - a file containing the 100 parameter combinations from the LHS according to the DFO criteria.  
* best_mean_params.csv - a file containing the 100 parameter combinations from the LHS according to the MLFO criteria. 
* col_germ.csv - a file containing the Gompertz fit for each minute for col. 
* ler_germ.csv - a file containing the Gompertz fit for each minute for ler. 
* emp_results.xlsx -  a file containing the wet lab experiments for ler dominance. 
* LHS.csv - the parameter combinations used in our LHS. 
* ler_dominance_data_MLFO/DFO - synthetic data generated from ler dominance experiments in NetLogo using MLFO/DFO critera for the best parameters.

Repo 2: https://data.mendeley.com/datasets/n6vwc2prv4/draft?a=b3795357-ae0b-461f-bce4-d1cf6b006561

* calibrate_data_double_accession - output file from NetLogo for model calibration (many realizations of many parameter combinations) in double accession experiments. 
* calibrate_data_single_accession - output file from NetLogo for model calibration (many realizations of many parameter combinations) in single accession experiments. 
