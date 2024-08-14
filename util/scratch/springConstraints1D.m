function [residualArray,xnew, lambda]= springConstraints1D(residualArray,elasticSprings,springVertices,xnew, lambda, mass, restLength, k, timeStep)
    residualArray = [residualArray,residualArray(:,end)];
    for j = 1:numel(elasticSprings)
        spring = elasticSprings(j);
        springVert = springVertices(spring,:);
        springx = xnew(springVert);
    
        mass1 = mass(springVert(1));
        mass2 = mass(springVert(2));
        w1 = 1/mass1;
        w2 = 1/mass2;
        if nargout > 1
            [deltaLambda, deltax1, deltax2, constraintC, alphaTilde]=springConstraint.project(timeStep,springx(1),springx(2),w1,w2,mass1,mass2,lambda(spring),0,0,restLength(spring),k(spring));
            
            xnew(springVert) = xnew(springVert)+[deltax1,deltax2];
            lambda(spring) = lambda(spring) + deltaLambda;
        
            residualArray(spring,end)=computeResidualh(constraintC,lambda(j),alphaTilde);
        else
            [constraintC, alphaTilde]=springConstraint.evaluate(timeStep,springx(1),springx(2), 0,restLength(spring),k(spring));
            residualArray(spring,end)=computeResidualh(constraintC,lambda(j),alphaTilde);
        end
    end
end