% Function to sample the parameter space and calculate coefficients 
% Markus Janasch, Ph.D. Student, KTH;
% Created: 2021-08-19, last modified: 2021-08-19


function KX_WriteDataToText(DataOut,Date,NrOfMetDataSet)

% INPUT:
%
% * DataOut - Matlab structure containing all output data, including:
%
% * CJ_rec  - 3-dimensional matrix of the scaled flux control coefficients
%             for all parameter sets resulting in the steady-state being
%             stable. CJ_rec(i,j,z) is the control exerted by enzyme j upon
%             flux i given the parameters sampled at iterazion z.
%
% * CS_rec  - 3-dimensional matrix of the scaled concentration control
%             coefficients for all parameter sets resulting in the steady-
%             state being stable.
%
% * E_rec   - 3-dimensional matrix of the scaled elasticity coefficients
%             for all parameter sets resulting in the steady-state being
%             stable.
%
% * MaxRealEigens -
%             Vector containing the maximal real part of the eigenvalues
%             for each parameter sampling iteration.
%
% * Parameters - 
%             2-dimensional matrix containing all the system's parameters
%             for each sampling iteration. Each row correspond to a
%             sampling iteration and each column to a parameter.
%
% OUTPUTS:
%
% * Text-file containing the percentage of stable steady states for the
% corresponding fMCS
%
% * Text-file containing the Flux Control Coefficients (FCCs) for each
% parameter set for the corresponding fMCS




%% Write Stability Percentage
Stability_Percent = (sum(DataOut.StabilityIndicator)/length(DataOut.StabilityIndicator))*100;

filename_1 = ['Stability_Percent_',Date,'_fMCS_',int2str(NrOfMetDataSet),'.txt'];
fileID_1 = fopen(filename_1,'w');
fprintf(fileID_1,'Stability_Percent:\t%g',Stability_Percent);
%fprintf(fileID,'\n');
%fprintf(fileID,'%g',Stability_Percent);
fclose(fileID_1);


%% Write FCCs
%filename_2 = ['FCCs_',Date,'_fMCS_',int2str(NrOfMetDataSet),'.tab'];
%fileID_2 = fopen(filename_2,'w');

FCC_Data = DataOut.CJ_rec;
[m,n,NrStableSets] = size(FCC_Data);

%% First column is ID of fMCS
fMCS_col = ones(n,1)*NrOfMetDataSet;

%% Second column is ID of stable parameter set
% going to be defined in the loop below

%% Third column is ID of Reaction
Rxn_col(:,1) = 1:n;
Output = [];
for i=1:NrStableSets % Iterate over stable parameter sets
    SPS_col = ones(n,1)*i;
    
    %% Round FCCs to make output files smaller
    for p = 1:m
        for q = 1:n
            FCC_Data_rounded(p,q) = round(FCC_Data(p,q,i),4);
        end
    end

    Output = [Output;fMCS_col SPS_col Rxn_col FCC_Data_rounded];
    
    
end
    outfile = ['FCCs_',Date,'_fMCS_',int2str(NrOfMetDataSet),'.tab'];
    dlmwrite(outfile, Output,'delimiter', '\t');
end