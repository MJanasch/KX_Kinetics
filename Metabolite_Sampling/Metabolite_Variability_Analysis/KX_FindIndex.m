%% FindIndex Function
% Modified from Avlant Nilsson
% Markus Janasch, Ph.D. Student
% Microbial Metabolic Engineering Group
% KTH School of Biotechnology, Science For Life Laboratory, Stockholm


function [n] = KX_FindIndex(haystack, needle)
    n=find(ismember(haystack,needle));
end
