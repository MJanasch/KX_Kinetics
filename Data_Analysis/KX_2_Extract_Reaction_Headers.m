% Extract the reactions header from a CBB SKM model file

load(Model_Input)
header = '';
for i = 1:length(N.reaction)
    name = char(N.reaction(i).id);
    header{i} = name;
    clear name
end
filename = 'KX_Reaction_Header.txt';
fileID = fopen(filename,'w');
for i=1:length(N.reaction)
    fprintf(fileID,'%s\t',header{i});
end
fclose(fileID);
