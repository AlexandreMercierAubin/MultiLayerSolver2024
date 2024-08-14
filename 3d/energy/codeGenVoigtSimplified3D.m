clear;
x = sym( 'x', [12,1], 'real' );
DmInv = sym( 'DmInv', [3,3], 'real' );
minv = sym( 'minv', [1,6], 'real' );

p0 = x(1:3);
p1 = x(4:6);
p2 = x(7:9);
p3 = x(10:12);

Ds = [p1-p0,p2-p0,p3-p0];
F = Ds*DmInv;
E = 0.5 *( F' * F - eye(3) );
constraintC = [E(1,1);E(2,2);E(3,3);E(1,2);E(1,3);E(2,3)];

%need to generate code for all 3 constraints of constraintC
%here I just show the first one as an example
gradC = [gradient(constraintC(1),x)';gradient(constraintC(2),x)';gradient(constraintC(3),x)';gradient(constraintC(4),x)';gradient(constraintC(5),x)';gradient(constraintC(6),x)'];
minvfull = reshape([minv;minv],[],1);
gradCM = gradC*minvfull;

out = [reshape(gradC',[],1); constraintC]; %reshape(gradCM',[],1); 
% fileOutName = 'STVKVoigtConstraintGrad.c';
% ccode( out, 'File', fileOutName );
matlabFunction(out, "File","StVKconstraint3D.m");

counter = 0;
fileID = fopen("StVKconstraint3D.m");
newFile = fopen("StVKconstraintOut3D.m","w");
gradCcount = 18;
gradCMcount = gradCcount*2;
while ~feof(fileID)
    tline = fgetl(fileID);

    for i =12:-1:1
        tline = strrep(tline,"x"+string(i), "x("+string(i)+",:,:)");
    end
    for i =1:6
        tline = strrep(tline,"minv"+string(i), "minv("+string(i)+",:,:)");
    end
    for i =1:9
        tline = strrep(tline,"DmInv"+string(mod(i-1,3)+1)+"_"+string(floor((i-1)/3)+1), "DmInv("+string(mod(i-1,3)+1)+","+string(floor((i-1)/3)+1)+",:)");
    end
   
    fprintf(newFile,tline+"\n");
end
fclose(fileID);
fclose(newFile);