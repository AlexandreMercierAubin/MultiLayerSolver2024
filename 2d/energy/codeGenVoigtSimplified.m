clear;
x = sym( 'x', [6,1], 'real' );
DmInv = sym( 'DmInv', [2,2], 'real' );
minv = sym( 'minv', [1,6], 'real' );

p0 = x(1:2);
p1 = x(3:4);
p2 = x(5:6);

Ds = [p0-p2,p1-p2];
F = Ds*DmInv;
E = 0.5 *( F' * F - eye(2) );
constraintC = [E(1,1);E(2,2);E(1,2)];

% b1 = (DmInv(1,:).*F)';
% b2 = (DmInv(2,:).*F)';
% b3 = simplify(-b1-b2); %-(DmInv(2,:)+DmInv(1,:)).*F

% permDmInv = DmInv(:,[2,1]);
% r31 = (0.5*sum(permDmInv(1,:).*F,2))';
% r32 = (0.5*sum(permDmInv(2,:).*F,2))';
% sumDm = sum(permDmInv,1);
% r33 = -(0.5*sum(sumDm.*F,2))';
% blockOfGradC = [b1, b2, b3;
%                 r31, r32, r33];

%need to generate code for all 3 constraints of constraintC
%here I just show the first one as an example
gradC = [gradient(constraintC(1),x)';gradient(constraintC(2),x)';gradient(constraintC(3),x)'];
gradCM = gradC.*minv;

out = [reshape(gradC',[],1); constraintC]; %reshape(gradCM',[],1); 
% fileOutName = 'STVKVoigtConstraintGrad.c';
% ccode( out, 'File', fileOutName );
matlabFunction(out, "File","StVKconstraint.m");

counter = 0;
fileID = fopen("StVKconstraint.m");
newFile = fopen("StVKconstraintOut.m","w");
gradCcount = 18;
gradCMcount = gradCcount*2;
while ~feof(fileID)
    tline = fgetl(fileID);

    for i =1:6
        tline = strrep(tline,"x"+string(i), "x("+string(i)+",:,:)");
    end
    for i =1:6
        tline = strrep(tline,"minv"+string(i), "minv("+string(i)+",:,:)");
    end
    for i =1:4
        tline = strrep(tline,"DmInv"+string(mod(i-1,2)+1)+"_"+string(floor((i-1)/2)+1), "DmInv("+string(mod(i-1,2)+1)+","+string(floor((i-1)/2)+1)+",:)");
    end
   
    fprintf(newFile,tline+"\n");
end
fclose(fileID);
fclose(newFile);

% counter = 0;
% fileID = fopen(fileOutName);
% newFile = fopen("STVKVoigtConstraintGradRewritten.c","w");
% gradCcount = 18;
% gradCMcount = gradCcount*2;
% while ~feof(fileID)
%     tline = fgetl(fileID);
% 
%     for i =1:6
%         tline = strrep(tline,"x"+string(i), "x["+string(i-1)+"]");
%     end
%     for i =1:4
%         tline = strrep(tline,"DmInv"+string(mod(i-1,2)+1)+"_"+string(floor((i-1)/2)+1), "DmInv["+string(i-1)+"]");
%     end
% 
%     %reading directly from the input source is likely a bit faster than
%     %copying so here we go!
%     tline = strrep(tline,"  t", "double t");
% 
%     A0pos = strfind(tline,'A0');
%     if ~isempty(A0pos)
%         if counter < gradCcount
%             tline = strrep(tline,"A0["+string(counter)+"][0]", "gradC"+"["+string(counter)+"]");
%         elseif counter < gradCMcount
%             tline = strrep(tline,"A0["+string(counter)+"][0]", "gradCM"+"["+string(counter)+"]");
%         elseif counter < gradCMcount+3
%             tline = strrep(tline,"A0["+string(counter)+"][0]", "constraintC"+"["+string(counter-gradCcount)+"]");
%         end
%         counter = counter + 1;
%     end
% 
%     fprintf(newFile,tline+"\n");
% end
% fclose(fileID);
% fclose(newFile);