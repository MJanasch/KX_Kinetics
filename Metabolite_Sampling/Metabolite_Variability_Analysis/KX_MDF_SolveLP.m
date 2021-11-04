function [solution]=KX_MDF_SolveLP(dG0,S,RT,ConcRanges_Model,Ratios,Met_Names_Model)%,beq)

[m,n] = size(S);
% m metabolite
% n reactions

%% Define Objective
c = [
    zeros(m,1);
    1
    ];

%% Define S_Ratio matrices
% S_ratios has to be Nr_Rations x Nr_metabolites, here 4xn
Nr_Ratios = length(Ratios.RatioNames);

S_ratios = zeros(m,Nr_Ratios);
for j = 1:Nr_Ratios
    Ratio_Split = split(Ratios.RatioNames(j),"/");
    Over = Ratio_Split(1);
    Under = Ratio_Split(2);
    S_ratios(KX_FindIndex(Met_Names_Model,Over),j) = 1;
    S_ratios(KX_FindIndex(Met_Names_Model,Under),j) = -1;
end

%% Define max and min vector for ratios

% Define A
A = [
     transpose(S),          ones(n,1);
     eye(m),                zeros(m,1);
    -eye(m),                zeros(m,1);
     transpose(S_ratios),   zeros(Nr_Ratios,1);
    -transpose(S_ratios),   zeros(Nr_Ratios,1)
    ];

% Define b
b = [
    -transpose(dG0)/RT;
     log(ConcRanges_Model(:,2));
    -log(ConcRanges_Model(:,1));
     log(Ratios.RatioRanges(:,2));
    -log(Ratios.RatioRanges(:,1))
    ];

solution.c=c;
solution.A=A;
solution.b=b;


% Primal Problem: Finding MDF
[X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = linprog(transpose(-c),A,b);

solution.z=X;

end
