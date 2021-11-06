![alt text](KX_Project_Logo_2.png "Interpreting regulatory metabolite-protein interactions in Synechocystis with kinetic metabolic modelling.")

# KX - Interpreting regulatory metabolite-protein interactions in Synechocystis with kinetic metabolic modelling.

Code for computational analysis of the kinetic metabolic model used in Sporre & Karlsen et al., 2021

---

### 1. Flux Calculations

##### Scripts

'Flux_Calculations/KX_Flux_Calculations.m' - Script to load the genome-scale metabolic model of Synechocystis 6803 [Sarkar et al., 2019](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006692), and perform FBA and subsequent flux sampling.

Requires installed toolboxes COBRA [Link](https://opencobra.github.io/cobratoolbox/latest/) and Raven 2.0 [Link](https://github.com/SysBioChalmers/RAVEN)

##### Input data: 
'Models/GEMs/SYN_GEM_iSyn731_CycleSyn.mat' - model file


##### Output:
Flux distribution from flux sampling
Note: Requires manual curation to adjust fluxes to small-scale model. Adjusted fluxes used in the parameter sampling are found in the model files 'Model/KX_Expanded_Model_Base.mat' and others.


---

### 2. Metabolite Concentration Variability Analysis
Find initial concentrations for Hit-and-Run Sampling of feasible metabolite concentration sets (fMCSs).

##### Scripts:
'Metabolite_Sampling/Metabolite_Variability_Analysis/KX_MDF_MetaboliteVariabilityAnalysis.m' - Main file to perform MDF analysis and subsequent metabolite concentration variability analysis to find many metabolite concentration sets to start HnR-Sampling from.

Support scripts:
- 'Metabolite_Sampling/Metabolite_Variability_Analysis/KX_FindIndex.m' - Helpful script to find indeces of variables in a list
- 'Metabolite_Sampling/Metabolite_Variability_Analysis/KX_MDF_SolveLP.m' - Script to solve the linear programming problem for MDF and max/min metabolite concentrations
- 'Metabolite_Sampling/Metabolite_Variability_Analysis/KX_MDF_VA.m' - Script to perform actual variability analysis

##### Input data:
'Metabolite_Sampling/Metabolite_Variability_Analysis/concentration_ranges_MDF.txt' - Metabolite concentration ranges, adapted from Asplund-Samuelsson et al. 2018 [link](https://www.sciencedirect.com/science/article/pii/S1096717617301076)

'Model/KX_Expanded_Model_Base.mat' - Kinetic model structure, needed for reactions, S-matrix, Keq and directionality

'Metabolite_Sampling/Metabolite_Variability_Analysis/Met_Order_Model.txt' - Metabolite names in the order required for the HnR-Sampling


##### Output:
Initial metabolite concentrations sets: a set of concentrations for all metabolites is created by maximizing and minimizing each metabolite concentration without violating the MDF value (to a certain degree, "Leeway")


---

### 3. Hit-And-Run Metabolite Sampling

##### Scripts:
'Metabolite_Sampling/HnR_Sampling/KX_HnR_1_MetaboliteSampling.py' - Hit-And-Run Sampling for metabolite concentrations, needs to be performed for each initial metabolite concentration set from 2. individually (can be looped in bash). Each sampling runs for 1 million steps, records about every 1000th step, until about 5000 fMCS are created. Resulting in total in about 400K fMCS. 

'Metabolite_Sampling/HnR_Sampling/KX_HnR_2_Create_Final_fMCSs.py' - Take in all metabolite concentration sets (about 400K) and subsample them randomly to create a final set of 5000 metabolite concentrations for each metabolite.

'Metabolite_Sampling/HnR_Sampling/KX_HnR_3_PoolSampling.py' - Sample the pool sizes for unbalanced metabolite supply, based on the metabolite concentration they supply.


Support scripts:
- 'Metabolite_Sampling/HnR_Sampling/read_concentration_ranges.py' - Reads in concentration ranges
- 'Metabolite_Sampling/HnR_Sampling/read_equilibrator_dg.py' - Reads in equilibrator results
- 'Metabolite_Sampling/HnR_Sampling/SMatrix_from_reactions.py' - Creates S-matrix from reactions


##### Input data:
'Metabolite_Sampling/HnR_Sampling/concentration_ranges.txt' - Metabolite concentration ranges in which to sample

'Metabolite_Sampling/HnR_Sampling/equilibrator_results.tsv' - Equilibrium constants/dG0 retrieved from eQuilibrator

'Metabolite_Sampling/HnR_Sampling/reactions.txt' - Metabolic reactions of the kinetic model

Additionally, the initial metabolite concentration sets from 2. are required. as an example, there is 'Metabolite_Sampling/HnR_Sampling/ConcSet_Init_EXAMPLE.txt'


##### Output:
Output after performing all three scripts consecutively is a file containing 5000 fMCSs, serving as input for the parameter sampling. An example output is found in 'Data/Metabolite_Concentration_Sets.txt'.


---

### 4. Parameter Sampling/Ensemble Modelling
Perform parameter sampling around the metabolic state defined by the fluxes and metabolite concentrations.
There are four models, and the scripts perform similar for each model, but are individual for the names they call. The description of the scripts is therefore exemplified for the Base model, but can be applied for all model and sampling versions.

##### Scripts:
'Parameter_Sampling/Base/KX_ScriptToCallParameterSampling_III.m' - Script that loads the model, metabolite concentrations and defines number of samplings ("iterations"). Can be looped in bash, or in the matlab file itself, to go through all fMCS.

'Parameter_Sampling/Base/KX_Parameter_Sampling_III.m' - Actual parameter sampling script, performs sampling of Km values around metabolite concentrations, calculates eigenvalues of Jacobian and flux control coefficients

'Parameter_Sampling/Base/KX_CalDFODC.m' - Script that calculates the values of the Jacobian.

'Parameter_Sampling/Base/KX_ManualCuration_Promiscuity.m' - Script that ensures that promiscuous enzymes have the same Km values.

'Parameter_Sampling/KX_WriteDataToText.m' - Script to transform output (FCCs and stability indicators) into text-files for each fMCS


##### Input data:
'Model/KX_Expanded_Model_Base.mat' - Kinetic model structure, containing all equations to be parameterized, fluxes and S-matrix

Metabolite concentration sets derived from 2. (Example: 'Data/Metabolite_Concentration_Sets.txt').


##### Output:
For each fMCS:
- one text file indicating the percentage of stable parameter sets (of the 1000 samplings)
- one text file with all flux control coefficients of all stable parameter sets

---

### 5. Data Analysis I - Concatenate Stability

##### Scripts:
'Data_Analysis/KX_1_TransformStabilityPercentage.m' - Script that loops through all stability-indicator files and transforms them into one single file

##### Input data:
Stability-indicator text files from the parameter sampling performed in 4.

##### Output:
A .tab-file (for each model) containing percentage of stable parameter sets for each fMCS. Results are found under 'Data/KX_Results_1_met_set_vs_percent_steady_Base.tab' (also for the other 3 model variants).


---

### 6. Data Analysis II - Concatenate FCCs and Calculate Median/MAD
First, the reactions are extracted from the model, required to create the header for the concatenated FCC files.

Second, the FCC files created for each fMCS are read in, and all FCCs for each reaction are combined in one file each, resulting in a FCC-file for each reaction, containing all FCCs of this reaction.

Third, the median and associated MAD over all FCCs for each reaction is calculated.

##### Scripts:
'Data_Analysis/KX_2_Extract_Reaction_Headers.m' - Script that extracts the reactions from the kinetic model input structure.

'Data_Analysis/KX_3_Concatenate_FCCs_Add_Header.py' - Script that concatenates all FCCs of a reactions into individual FCC files for each reaction.

'Data_Analysis/KX_4_FCCs_MEDMAD.R' - Script that reads in the individual FCC files for each reaction and calulates mthe median and MAD over all values.


##### Input:
For KX_2_Extract_Reaction_Headers.m: model structure, 'Model/KX_Expanded_Model_Base.mat'

For KX_3_Concatenate_FCCs_Add_Header.py: FCC files from parameter sampling in 4., needs to be performed for each reaction. Can be looped.

For KX_4_FCCs_MEDMAD.R: FCC files for each reaction. Looped internally.

##### Output:
For KX_2_Extract_Reaction_Headers.m: 'KX_Reaction_Header.txt' (Example found in Data-folder)

For KX_3_Concatenate_FCCs_Add_Header.py: For each reaction, one text file containing all FCCs for this reaction.

For KX_4_FCCs_MEDMAD.R: One .csv-file containing all median and MAD values for all reactions. Results are found in the Data-folder ('KX_Results_3_FCCs_Med_MAD_Base.csv').

---

### 7. Data Analysis III - Plot Metabolome Stability and perform/plot Kolmogorov-Smirnov test for all models

##### Scripts:
'Data_AnalysisKX_5_Plot_Concs_vs_Stability_And_KS-Test.R' - Script that reads in .tab-files containing stability percentages for all model variants and plots the distribution of metabolite concentrations for the top and bottom 10% in stability percentages. Furthermore, a Kolmogorov-Smirnov (KS) Test was performed for the distributions of the top 10% in stability percentages between the model variants. "Significance" was decided as having a p-value below 0.05.

##### Input:
'Data/Metabolite_Concentration_Sets.txt' - Metabolite concentration sets (fMCSs)

'Data/KX_Results_1_met_set_vs_percent_steady_Base.tab' - percentage of stable parameter sets for each fMCS for the base model

'Data/KX_Results_1_met_set_vs_percent_steady_FSBPase.tab' - percentage of stable parameter sets for each fMCS for the FSBPase model

'Data/KX_Results_1_met_set_vs_percent_steady_TKT.tab' - percentage of stable parameter sets for each fMCS for the TKT model

'Data/KX_Results_1_met_set_vs_percent_steady_BOTH.tab' - percentage of stable parameter sets for each fMCS for the BOTH model


##### Output:
Ridgeline plot (Example found in 'Data/KX_Results_4_Conc_Stab_All.pdf')

KS Test results for each comparison as a large, colored table (Example found in Data/KX_Results_5_KS-Test_Stable.pdf)

---

### 8. Data Analysis III - Plot FCCs for all Models

##### Scripts:
'Data_Analysis/KX_6_Plot_FCCs_Checkered.R' - Script that reads in .csv-files containing median and MAD values for FCCs for all four model variants and plotting them in a checkered, subsquared heatmap.

##### Input:
.csv-files containing median and MAD values for the FCCs:
'Data/KX_Results_3_FCCs_Med_MAD_Base.csv',
'Data/KX_Results_3_FCCs_Med_MAD_BOTH.csv',
'Data/KX_Results_3_FCCs_Med_MAD_FSBPase.csv',
'Data/KX_Results_3_FCCs_Med_MAD_TKT.csv'

##### Output:
Heatmap with 2x2 subsquares representing each model variant (Example found in 'Data/KX_Results_6_All_FCCs_Supplement.pdf')


## Note


---

## Dependencies

CORBA Toolbox [Link](https://opencobra.github.io/cobratoolbox/latest/) 

RAVEN 2.0 [Link](https://github.com/SysBioChalmers/RAVEN)

Mathworks Matlab (Performed on MATLAB_R2018a)

Gurobi solver for FBA, free academic license (https://www.gurobi.com)

Python 3 (3.8.5) with the following libraries:
- sys, os, numpy, time, math, scipy

R version 3.6.1, with the libraries:
- data.table, reshape2, tidyverse, foreach, doMC, ggridges, scales, optparse, ggrepel, ggplot2

---

## Author

Markus Janasch, SciLifeLab/KTH (markus.janasch@scilifelab.se)