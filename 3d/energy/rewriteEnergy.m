%just a script that handles the painful task of reformating the energy code
%generated with the symbolic lib

fileToChange = "NeoHookeanPsiGradHess.c";

missingDouble = "  t";
addDouble = "double t";

Gridsize = 9;
HessianName = "C";
GradientName = "dpsidF";

useIndices = false;

counter = 0;
fileID = fopen(fileToChange);
newFile = fopen("output.c","w");
while ~feof(fileID)
    tline = fgetl(fileID);
    tline = strrep(tline,missingDouble, addDouble);

    A0pos = strfind(tline,'A0');
    if ~isempty(A0pos)
        tline = regexprep(tline,"\[(.*?)\]", ""); %removes the bracket stuff
        if counter < Gridsize*Gridsize
            x = mod(counter,Gridsize);
            y = floor(counter/Gridsize);
            if x == 0
                fprintf(newFile,"\n");
            end
            if useIndices == true
                tline = strrep(tline,"A0", HessianName + "["+string(x)+"]"+"["+string(y)+"]");
            else
                tline = strrep(tline,"A0", HessianName + "[i++]");
            end
        elseif counter < Gridsize*Gridsize + Gridsize
            x = counter - Gridsize*Gridsize;
            if useIndices == true
                tline = strrep(tline,"A0", GradientName + "["+string(x)+"]");
            else
                tline = strrep(tline,"A0", GradientName + "[j++]");
            end
        else
            tline = strrep(tline,"A0", "psi");
        end
        counter = counter + 1;
    end

    fprintf(newFile,tline+"\n");
end
fclose(fileID);
fclose(newFile);