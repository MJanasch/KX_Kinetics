%% Check EVA for concentrations

function [Met_Ranges,Met_Ranges_All] = KX_MDF_VA(S_Used,MDF,RT,A,b,Leeway,Original_Ranges)

    MDF_EVA_Orig = MDF;
    MDF_EVA_Leeway = MDF*Leeway; % Leeway = 0.9
    [m,n] = size(S_Used); % m metabolites, n reactions
    
    % loop through metabolites, maximize/minimize them
    for Met_Index = 1:m
        % Define b
        b_new=b;
        b_new(n+Met_Index) = log(Original_Ranges(Met_Index,2));
                
        % Maximize [Met]
        mima = 1; 
        
        % MDF at optimum
        %[Results_Min_Orig] = FindExtreme_EVA(RT,MDF_EVA_Orig,mima,m,Met_Index,A,b);
        
        % MDF with margin
        [Results_Min_Leeway,Conc_Out_All_Min] = FindExtreme_EVA(RT,MDF_EVA_Orig,mima,m,Met_Index,A,b_new);
        
        %Conc_Out_Min = min(Results_Min_Orig,Results_Min_Leeway);
        Conc_Out_Min = Results_Min_Leeway;
        
        
        % Minimize [Met]
        mima = -1;
        
        % MDF at optimum
        %[Results_Max_Orig] = FindExtreme_EVA(RT,MDF_EVA_Orig,mima,m,Met_Index,A,b);
        
        % MDF with margin
        [Results_Max_Leeway,Conc_Out_All_Max] = FindExtreme_EVA(RT,MDF_EVA_Leeway,mima,m,Met_Index,A,b);
        
        %Conc_Out_Max = min(Results_Max_Orig,Results_Max_Leeway);
        Conc_Out_Max = Results_Max_Leeway;
        
        
        Met_Ranges(Met_Index,1) = Conc_Out_Min;
        Met_Ranges(Met_Index,2) = Conc_Out_Max;
        
        Met_Ranges_All(Met_Index).Conc_Out_All_Min = Conc_Out_All_Min;
        Met_Ranges_All(Met_Index).Conc_Out_All_Max = Conc_Out_All_Max;
    end

end

function [Conc_Out,Conc_Out_All] = FindExtreme_EVA(RT,MDF_EVA,mima,Nr_Met_EVA,Met_Index,A,b)   
%% Define conc vector to min/max
    Conc_Vector = zeros(Nr_Met_EVA,1);
    
    if mima == 1
        % Set value in c to 1
        Conc_Vector(Met_Index,1) = 1;

        c = [Conc_Vector;
            0];
        
    elseif mima == -1
        % Set value in c to -1
        Conc_Vector(Met_Index,1) = -1;
        
        c = [Conc_Vector;
            0];
    else
        disp('minmax needs to either 1 or -1 to find the maximum or minimum of the concentration range, respectively.')
    end
    
    %% Set Aeq and beq so that B = MDF
    Aeq = [zeros(1,Nr_Met_EVA) 1];
   
    % Define beq
    beq = MDF_EVA/RT;
    
    % Solve
    options = optimoptions('linprog','Display','none');
    X = linprog(transpose(c),A,b,Aeq,beq,[],[],options);
    Conc_Out = exp(X(Met_Index)); % [M]
    Conc_Out_All = exp(X(1:Nr_Met_EVA)); % [M]
end