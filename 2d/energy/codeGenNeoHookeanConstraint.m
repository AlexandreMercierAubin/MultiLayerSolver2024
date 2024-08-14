clear;
x = sym( 'x', [2*3,1], 'real' );
DmInv = sym( 'DmInv', [2,2], 'real' );
mass = sym( 'mass', [2*3,1], 'real' );
materialMu = sym( 'materialMu' , 'real' );
materialLambda = sym( 'materialLambda' , 'real' );
area = sym( 'area' , 'real' );
h = sym( 'h' , 'real' );
lambda = sym( 'lambda' , 'real' );
gamma = sym( 'gamma' , 'real' );

p0 = x(1:2);
p1 = x(3:4);
p2 = x(5:6);

%computing deformation gradient
Ds  = [p0-p2,p1-p2];
F = Ds*DmInv;

%stvk energy
% gamma = 1 + (materialMu/materialLambda);
C_H = det(F)-gamma; %slightly after eq 14
alphaH = 1/(area*materialLambda); %equation 13 note that the paper is confusing and uses lambda for both the lam'e parameters and the Lagrange multiplyers

FtF = F'*F;
C_D = sqrt(trace(FtF)); %Eq 8
alphaD = 1/(area*materialMu);

%I don't need those for constraint based programming.
% C_H2 = C_H*C_H;
% psiH = 0.5 * materialLambda * area * C_H2;
% psiD = 0.5 * materialMu * area * C_D;

%grads from Supplemental to A Constraint-based Formulation of Stable Neo-Hookean Materials
gradCH = gradient(C_H, x);
%not sure what is the 2D equivalent for this one...
% f2xf3 = cross2D(F(:,2),F(:,3));
% f3xf1 = cross2D(F(:,3),F(:,1));
% f1xf2 = cross2D(F(:,1),F(:,2));
% gradCH = [f2xf3,f3xf1,f1xf2]*DmInv';

% gradCD = gradient(C_D, x);
rs = sqrt(sum(F(:).^2));
gradCDtop = (1/rs) * F * DmInv';
gradCDbot = -sum(gradCDtop,2);
gradCD = [gradCDtop(:);gradCDbot];

h2 = h*h;
alphaTildeH = alphaH/h2;
alphaTildeD = alphaD/h2;

minv = 1./mass;
Minv = diag(minv);

numeratorH = -C_H - alphaTildeH * lambda;
denumeratorH =  gradCH'*Minv*gradCH + alphaTildeH;
deltaLambdaH = numeratorH/denumeratorH;
deltaxH = minv.*gradCH*deltaLambdaH;

numeratorD = -C_D - alphaTildeD * lambda;
denumeratorD = gradCD'*Minv*gradCD + alphaTildeD;
deltaLambdaD = numeratorD/denumeratorD;
deltaxD = minv.*gradCD*deltaLambdaD;

% asking matlab to make a file for all parts simultaneously allows for a
% better use (elimination) of common sub expressions.  Here, the force is 
% put following Hessian so it is easy to pull off the end of the code gen.
outD = [ reshape(deltaxD,[],1); minv; reshape(gradCD,[],1); deltaLambdaD ; alphaTildeD; C_D]; 
fileOutNameD = 'neoHookeanDconstraint.c';
ccode( outD, 'File', fileOutNameD );

outH = [ reshape(deltaxH,[],1); minv; reshape(gradCH,[],1); deltaLambdaH ; alphaTildeH; C_H]; 
fileOutNameH = 'neoHookeanHconstraint.c';
ccode( outH, 'File', fileOutNameH );

rewriteFile(fileOutNameD);
rewriteFile(fileOutNameH);

function rewriteFile(fileOutName)
    counter = 0;
    fileID = fopen(fileOutName);
    revisedName = strcat('Revised',fileOutName);
    newFile = fopen(revisedName,"w");
    while ~feof(fileID)
        tline = fgetl(fileID);
        for i =1:6
            tline = strrep(tline,"x"+string(i), "x["+string(i-1)+"]");
        end
        
        for i =1:6
            tline = strrep(tline,"xOld"+string(i), "xOld["+string(i-1)+"]");
        end
    
        for i =1:6
            tline = strrep(tline,"mass"+string(i), "mass["+string(i-1)+"]");
        end
    
        for i =1:4
            tline = strrep(tline,"DmInv"+string(mod(i-1,2)+1)+"_"+string(floor((i-1)/2)+1), "DmInv["+string(i-1)+"]");
        end
    
        %reading directly from the input source is likely a bit faster than
        %copying so here we go!
        tline = strrep(tline,"materialMu", "materialMu[0]");
        tline = strrep(tline,"materialLambda", "materialLambda[0]");
        tline = strrep(tline,"volume", "volume[0]");
        tline = strrep(tline,"alpha", "alpha[0]");%no clue why alpha ends up alph[0]a[0]
        tline = strrep(tline,"beta", "beta[0]");
        tline = strrep(tline,"h", "h[0]");
        tline = strrep(tline,"lambda", "lambda[0]");
        tline = strrep(tline,"  t", "double t");
    
        A0pos = strfind(tline,'A0');
        if ~isempty(A0pos)
            if counter < 6
                tline = strrep(tline,"A0["+string(counter)+"][0]", "deltax"+"["+string(counter)+"]");
            elseif counter < 12
                tline = strrep(tline,"A0["+string(counter)+"][0]", "Minv"+"["+string(counter-6)+"]");
            elseif counter < 18
                tline = strrep(tline,"A0["+string(counter)+"][0]", "gradC"+"["+string(counter-12)+"]");
            elseif counter < 19
                tline = strrep(tline,"A0["+string(counter)+"][0]", "deltaLambda[0]");
            elseif counter < 20
                tline = strrep(tline,"A0["+string(counter)+"][0]", "alphaTilde[0]");
            elseif counter < 21
                tline = strrep(tline,"A0["+string(counter)+"][0]", "constraintC[0]");
            end
            counter = counter + 1;
        end
       
        fprintf(newFile,tline+"\n");
    end
    fclose(fileID);
    fclose(newFile);
end