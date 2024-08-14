clear;
x = sym( 'x', [2*3,1], 'real' );
mass = sym( 'mass', [2*3,1], 'real' );
DmInv = sym( 'DmInv', [2,2], 'real' );
materialMu = sym( 'materialMu' , 'real' );
materialLambda = sym( 'materialLambda' , 'real' );
h = sym( 'h' , 'real' );
lambda = sym( 'lambda' , [3,1], 'real' );

alpha = inv([materialLambda + 2*materialMu,materialLambda,0;
         materialLambda, materialLambda + 2*materialMu, 0;
         0,0, 2*materialMu]);

h2 = h*h;
alphaTilde = alpha./h2;

minv = 1./mass;
Minv = diag(minv);

p0 = x(1:2);
p1 = x(3:4);
p2 = x(5:6);

Ds = [p0-p2,p1-p2];
F = Ds*DmInv;
E = 0.5 *( F' * F - eye(2) );
constraintC = [E(1,1);E(2,2);E(1,2)];


%need to generate code for all 3 constraints of constraintC
%here I just show the first one as an example
gradC = [gradient(constraintC(1),x)';gradient(constraintC(2),x)';gradient(constraintC(3),x)'];

%anything from here should be computed numerically... maybe just generate
%code for grads instead of this

%The top is coupled
topDenominator = gradC(1:2,:)*Minv*gradC(1:2,:)' + alphaTilde(1:2,1:2);
topNumerator = -constraintC(1:2) - alphaTilde(1:2,1:2)*lambda(1:2);
topDeltaLambda = inv(topDenominator)*topNumerator;

%This is problematic... alphaTilde should be a scalar
denominator = gradC(3,:)*Minv*gradC(3,:)' + alphaTilde(3,3);
numerator = -constraintC(3) - alphaTilde(3,3)*lambda(3);

deltaLambda = [topDeltaLambda; numerator/denominator];
deltax = Minv*gradC'*deltaLambda;

out = [ reshape(deltax,[],1); minv; reshape(gradC',[],1);deltaLambda ;  reshape(alphaTilde,[],1); constraintC]; 
out = simplify(out);
fileOutName = 'STVKVoigtConstraintSimplify.c';
ccode( out, 'File', fileOutName );

counter = 0;
fileID = fopen(fileOutName);
newFile = fopen("STVKVoigtConstraintSimplifyRewritten.c","w");
gradCcount = 12+3*6;
dlCount= gradCcount+3;
alphaTildeCount = dlCount + 9;
while ~feof(fileID)
    tline = fgetl(fileID);
    for i =1:6
        tline = strrep(tline,"x"+string(i), "x["+string(i-1)+"]");
    end
    
    for i =1:3
        tline = strrep(tline,"lambda"+string(i), "lambda["+string(i-1)+"]");
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
    tline = strrep(tline,"alpha", "alpha[0]");%no clue why alpha ends up alph[0]a[0]
    tline = strrep(tline,"beta", "beta[0]");
    tline = strrep(tline,"h", "h[0]");
    tline = strrep(tline,"  t", "double t");

    A0pos = strfind(tline,'A0');
    if ~isempty(A0pos)
        
        if counter < 6
            tline = strrep(tline,"A0["+string(counter)+"][0]", "deltax"+"["+string(counter)+"]");
        elseif counter < 12
            tline = strrep(tline,"A0["+string(counter)+"][0]", "Minv"+"["+string(counter-6)+"]");
        elseif counter < gradCcount
            tline = strrep(tline,"A0["+string(counter)+"][0]", "gradC"+"["+string(counter-12)+"]");
        elseif counter < dlCount
            tline = strrep(tline,"A0["+string(counter)+"][0]", "deltaLambda"+"["+string(counter-gradCcount)+"]");
        elseif counter < alphaTildeCount
            while isempty(strfind(tline,"A0["+string(counter)+"][0]"))
                counter = counter + 1;
            end
            tline = strrep(tline,"A0["+string(counter)+"][0]", "alphaTilde"+"["+string(counter-dlCount)+"]");
        elseif counter < alphaTildeCount+3
            tline = strrep(tline,"A0["+string(counter)+"][0]", "constraintC"+"["+string(counter-alphaTildeCount)+"]");
        end
        counter = counter + 1;
    end
   
    fprintf(newFile,tline+"\n");
end
fclose(fileID);
fclose(newFile);