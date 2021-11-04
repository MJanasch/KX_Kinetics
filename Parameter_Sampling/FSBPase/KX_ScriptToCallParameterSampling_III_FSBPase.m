% This script calls the parameter sampling and coefficient calculating 
% functions for the kinetic models
% Markus Janasch, Ph.D. Student, KTH;
% Created: 2017-03-22, last modified: 2021-10-11


% Get Folder
currentFolder = pwd;

%% Add path to input data
addpath('PATH_TO_MODEL');
addpath('PATH_TO_METABOLITE_DATA');

%% Input data
Iterations = 1000;
InputDataStructure = '../Models/KX_Expanded_Model_FSBPase.mat';
MetConcSamplingData = '../Data/Metabolite_Concentration_Sets.txt';
Date = ''; % Date of Sampling
NrOfMetDataSet = 1; % Number of feasible metabolite concentration set, 1-5000 (10000), can be put into loop

%% Load Metabolite Concentration Sets
MetConcData_RAW = importdata(MetConcSamplingData);

%% Define Output
Output_Folder_Name = ['OUTPUT_FOLDER_NAME_',Date,'_FSBPase'];

[DataOut] = KX_Parameter_Sampling_III_FSBPase(Iterations,InputDataStructure,MetConcData_RAW.data(NrOfMetDataSet,:),MetConcData_RAW.textdata);

cd(Output_Folder_Name);
%% KX: Save output not into matlab files, but rather text files
% This simplifies the data analysis as matlab is not required for it.
KX_WriteDataToText(DataOut,Date,NrOfMetDataSet)
cd(currentFolder)

