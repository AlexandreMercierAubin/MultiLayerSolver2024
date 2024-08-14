function mesh = mesh2DAsShellLoader(fileNameRoot, attributes, materials, scale, perturbation)
    %generateCloth(attributes, materials, precision, scale)
    %generates some 2D rectangular mesh and project it to 3D 
    if nargin < 4
        scale = [1,1,1];
    end
    if nargin < 5
        perturbation = 0;
    end
    
    [ p, I ] = readNODE( sprintf('%s.node', fileNameRoot ) );
    
    % I tells us the indices... must be a sequence and start from 1 or we 
    % are in trouble!
    
    [ T, A ] = readELE( sprintf('%s.ele', fileNameRoot ) );
         
    rands = randi(100,size(p,1),1);
    pert = -perturbation + rands*perturbation*2;

    V = [p(:,1),zeros(size(p,1),1),p(:,2)];
    V(:,1) = scale(1)*V(:,1);
    V(:,2) = scale(2)*V(:,2)+pert;
    V(:,3) = scale(3)*V(:,3);
    J = [1:size(T,1)]';
    mesh = Mesh3D(V,T, A, materials, T, J,[],[],[],[]);
end