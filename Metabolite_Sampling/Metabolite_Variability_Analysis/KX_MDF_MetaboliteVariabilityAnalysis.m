%% KX Metabolite Variability Analysis
% Markus Janasch
% For the KX metabolic model of core metabolism of Synechocystis 6803:
% 1. Perform MDF analysis
% 2. Use MDF analysis result to perform metabolite variability analysis as
% performed in Janasch et al., 2021
% created:21/04-23, last modified: 21/11-2

%% Load Data

R = 8.3145e-3;  % [kJ/mol*K]
T = 303.15;     % [K]
RT = R*T;

% Load model
InputDataStructure = '../Models/KX_Expanded_Model_Base.mat';
load(InputDataStructure)
clear InputDataStructure

%% Number of actual metabolic reactions
Nr_act_Rxns = 0;
for i=1:length(N.reaction)
    if ~strncmp(N.reaction(i).name,'Sink',4) && ~strncmp(N.reaction(i).name,'Supply',6)
        Nr_act_Rxns = Nr_act_Rxns + 1;
    end
end

%% Number of actual metabolites (Excludes biomass sink metabolites and supply pools)
Nr_act_Mets = 0;
for i=1:length(N.species)
    if ~startsWith(N.species(i).name,'BioM') && ~endsWith(N.species(i).name,'Pool')
        Nr_act_Mets = Nr_act_Mets + 1;
    end
end

%% Extract K_eq from model
for Rxn = 1:Nr_act_Rxns
    for p=1:length(N.reaction(Rxn).kineticLaw.localParameter)
        if strncmp(N.reaction(Rxn).kineticLaw.localParameter(p).id,'K_eq',4)
            K_eq(Rxn) = N.reaction(Rxn).kineticLaw.localParameter(p).value;            
        end
    end
end

% Manual curation of equilibrium constants for Rubisc/o and the lumped
% reactions for photosynthesis
K_eq(1) = 38327.61674;
K_eq(36) = 1.6681068660407923e+91;
K_eq(21) = 1.3123E+04;
K_eq(22) = 21.8;
K_eq(39) = exp(-(-8.1)/RT);
K_eq(40) = exp(-(0.7)/RT);

for Rxn = 1:Nr_act_Rxns
    dG0(Rxn)=(-RT*log(K_eq(Rxn))); 
end


%% Load metabolite ranges
ConcRanges_Raw_MDF = importdata('concentration_ranges_MDF.txt'); % Read in as M
ConcRanges_MDF = ConcRanges_Raw_MDF.data/1000; % Transform from mM to M in the process
ConcRanges_Names_MDF = ConcRanges_Raw_MDF.textdata(2:end,1);


% Extract cofactor ratios
NrRatios = 0;
for i=1:length(ConcRanges_Names_MDF)
    if contains(ConcRanges_Names_MDF(i),'/')
        NrRatios = NrRatios + 1;
        Ratios.RatioRanges(NrRatios,:) = ConcRanges_MDF(i,:)*1000;
        Ratios.RatioNames(NrRatios,1) = ConcRanges_Names_MDF(i);
    end
end

ConcRanges_Names_MDF(46) = []; % Remove PPool
ConcRanges_MDF(46,:) = []; % Remove PPool

ConcRanges_Names_MDF(32) = []; % Remove NADPH/NADP ratio
ConcRanges_MDF(32,:) = []; % Remove NADPH/NADP ratio

ConcRanges_Names_MDF(31) = []; % Remove ATP/ADP ratio
ConcRanges_MDF(31,:) = []; % Remove ATP/ADP ratio

ConcRanges_Names_MDF(30) = []; % Remove NADH/NAD ratio
ConcRanges_MDF(30,:) = []; % Remove NADH/NAD ratio

% Bring concentration ranges into model order
for i = 1:Nr_act_Mets
    Met_Names_Model{i} = N.species(i).id;
end
Met_Names_Model{end+1} = 'H2O';     % H2O is not in the model structure, so added here
Nr_act_Mets = length(Met_Names_Model);

ConcRanges_Model = ones(length(ConcRanges_Names),2);
ConcRanges_Model_MDF = ones(length(ConcRanges_Names_MDF),2);
for i = 1:length(ConcRanges_Model_MDF)
    Index_Met = KX_FindIndex(Met_Names_Model,ConcRanges_Names{i});
    ConcRanges_Model_MDF(Index_Met,:) = ConcRanges_MDF(i,:);
end

% Get proper S-matrix
S_Used = SFullFull(1:Nr_act_Mets,1:Nr_act_Rxns);
[m,n] = size(S_Used);


%% Max-Min-Driving Force calculation
% Set Solver
changeCobraSolver('matlab','LP'); % Employ COBRA's function to easily change LP solvers (I can choose between matlab, gurobi and mosek)

%% Employ MDF Algorithm

solutionLP = KX_MDF_SolveLP(dG0,S_Used,RT,ConcRanges_Model_MDF,Ratios,Met_Names_Model);%,beq);

A = full(solutionLP.A);
b = solutionLP.b;
c = solutionLP.c;

dG = transpose(dG0)+RT*transpose(S_Used)*solutionLP.z(1:m);
%dGtrans = transpose(dG);
MDF = solutionLP.z(end,1)*RT;
conc = round(exp(solutionLP.z(1:m))*1000, 5); % in [mM] via *1000
sum(conc/1000);
dG_rounded = transpose(dG0)+RT*transpose(S_Used)*log(conc/1000);
%% Variability Analysis
Leeway = 0.9; % Factor to which degree the metabolite ranges are to represent the MDF, 0 < Leeway =< 1

[VA_MetRanges,VA_MetRanges_All] = KX_MDF_VA(S_Used,MDF,RT,A,b,Leeway,ConcRanges_Model_MDF); % Metabolite concentration ranges from Variability analysis
% Can remove fixed concentrations: H20, CO2_cax, CO2_cyt, O2, which are the last 4
%VA_MetRanges_All contains 86 concentration sets, where each metabolite
%concentration has been maximized and minimized while keeping the MDF (or a
%lower MDF value, depending on leeway). These can be used as starting 
% points for metabolite sampling via HnR

%% Need to order the metabolites according to python script for metabolite sampling
Met_Order_HnR_Script = importdata('Met_Order_Model.txt');
for i = 1:length(Met_Order_HnR_Script)
    Index_New(i,1) = KX_FindIndex(Met_Order_HnR_Script,Met_Names_Model(i));
end


%%
Met_Names_Model = transpose(Met_Names_Model);

for i=1:m
    Met_Names_Model_Sorted(Index_New(i),1) = Met_Names_Model(i);
    S_NewOrder(Index_New(i),:) = S_Used(i,:);
end
    
for i = 1:length(Met_Order_HnR_Script)
    Temp_Concs_Min(Index_New,1) = VA_MetRanges_All(i).Conc_Out_All_Min;
    Temp_Concs_Max(Index_New,1) = VA_MetRanges_All(i).Conc_Out_All_Max;
    sum(Temp_Concs_Min);
    sum(Temp_Concs_Max);
    
    
    VA_MetRanges_ALL_NewOrder(Index_New(i)).Conc_Out_All_Min = Temp_Concs_Min;
    VA_MetRanges_ALL_NewOrder(Index_New(i)).Conc_Out_All_Max = Temp_Concs_Max;
end

%% Output metabolite concentration range files to be input for random sampling
mkdir Output_for_HnR_ConcSet_Init
cd('Output_for_HnR_ConcSet_Init/');
    for i = 1:m
        filename = ['ConcSet_Init_',Met_Names_Model_Sorted{i},'_Min','.txt'];
        fileID = fopen(filename,'w');
        fprintf(fileID,'[');
        fprintf(fileID,' %g,',VA_MetRanges_ALL_NewOrder(i).Conc_Out_All_Min(1:end-1));
        fprintf(fileID,' %g]',VA_MetRanges_ALL_NewOrder(i).Conc_Out_All_Min(end));
        fclose(fileID);
        
        filename = ['ConcSet_Init_',Met_Names_Model_Sorted{i},'_Max','.txt'];
        fileID = fopen(filename,'w');
        fprintf(fileID,'[');
        fprintf(fileID,' %g,',VA_MetRanges_ALL_NewOrder(i).Conc_Out_All_Max(1:end-1));
        fprintf(fileID,' %g]',VA_MetRanges_ALL_NewOrder(i).Conc_Out_All_Max(end));
        fclose(fileID);
    end