%% Transform Stability percentage output into right format

%% Read in Data


%a = dir([indir,'/*.txt']);
%Nr_fMCSs = size(a,1)-2;
tic
Nr_fMCSs = 5000;
for i = 1:Nr_fMCSs
    dataraw = importdata([indir,'/Stability_Percent_',date,'_fMCS_',int2str(i),'.txt']);
    Stability_Percent_Results(i,1) = i;
    Stability_Percent_Results(i,2) = dataraw.data;
    clear dataraw
    toc
end

%% Write data
filename = 'KX_Results_1_met_set_vs_percent_steady.tab';
fileID = fopen(filename,'w');
for i = 1:Nr_fMCSs
    fprintf(fileID,'%g\t%g\n',Stability_Percent_Results(i,1),Stability_Percent_Results(i,2));
end
fclose(fileID);