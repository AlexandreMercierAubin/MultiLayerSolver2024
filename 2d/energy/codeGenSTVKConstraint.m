clear;
x = sym( 'x', [2*3,1], 'real' );
mass = sym( 'mass', [2*3,1], 'real' );
DmInv = sym( 'DmInv', [2,2], 'real' );
materialMu = sym( 'materialMu' , 'real' );
materialLambda = sym( 'materialLambda' , 'real' );
volume = sym( 'volume' , 'real' );
alpha = sym( 'alpha' , 'real' );
h = sym( 'h' , 'real' );
lambda = sym( 'lambda' , 'real' );

p0 = x(1:2);
p1 = x(3:4);
p2 = x(5:6);

Ds = [p0-p2,p1-p2];
F = Ds*DmInv;
E = 0.5 *( F' * F - eye(2) );
EtE = E'*E;
froE2 = EtE(1,1) + EtE(2,2);
trE = E(1,1)+E(2,2);
psi = materialMu * froE2 + materialLambda * 0.5 * trE*trE;
constraintC = psi*volume;


gradC = gradient(constraintC,x); %-------
% piolaStress = F*(2*materialMu*E + materialLambda*trE*eye(2));
% dpsidx = piolaStress * DmInv';
% W = -volume * dpsidx;
% force3 = -W(:,1) - W(:,2);
% forceStackH = [W,force3];
% gradC = -reshape(forceStackH,[], 1);


h2 = h*h;
alphaTilde = alpha/h2;

numerator = -constraintC - alphaTilde*lambda;
minv = 1./mass;
Minv = diag(minv);
denominator = gradC'*Minv*gradC + alphaTilde;
deltaLambda = numerator/denominator;
deltax = Minv*gradC*deltaLambda;

% asking matlab to make a file for all parts simultaneously allows for a
% better use (elimination) of common sub expressions.  Here, the force is 
% put following Hessian so it is easy to pull off the end of the code gen.
out = [ reshape(deltax,[],1); minv; reshape(gradC,[],1); deltaLambda ; alphaTilde; constraintC]; 
fileOutName = 'STVKconstraint.c';
ccode( out, 'File', fileOutName );

counter = 0;
fileID = fopen(fileOutName);
newFile = fopen("STVKconstraintRevised.c","w");
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