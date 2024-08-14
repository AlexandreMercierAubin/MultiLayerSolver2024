classdef NeoHookean3DEnergy < EnergyModel
    %CorotationalEnergy class that allows the computation of corotational
    %energy based on https://graphics.pixar.com/library/StableElasticity/paper.pdf
    properties
        
    end
    
    methods
        function obj = NeoHookean3DEnergy()
            obj@EnergyModel();
            obj.name = "neoHookean";
        end 
        
        function computeEnergy(obj, mesh3D, deformationGradientF)
            [ii, jj, CblockVals, dpsidF, psi] = mexNeoHookean3D( deformationGradientF, mesh3D.elV, mesh3D.elMu, mesh3D.elLambda );
            sizeC = 9*size(mesh3D.t,1);
            C = sparse( ii, jj, CblockVals,sizeC,sizeC);
            obj.elasticForces = mesh3D.B' * dpsidF; 
            obj.derivative1Gradient = dpsidF;
            obj.derivative2HessianC = C;
            obj.energy = psi;
        end
        
    end
    methods(Static)
        function psi = peekEnergy(mesh3D,deformationGradientF)
            F1_1 = deformationGradientF(1:9:end);
            F2_1 = deformationGradientF(2:9:end);
            F3_1 = deformationGradientF(3:9:end);
            F1_2 = deformationGradientF(4:9:end);
            F2_2 = deformationGradientF(5:9:end);
            F3_2 = deformationGradientF(6:9:end);
            F1_3 = deformationGradientF(7:9:end);
            F2_3 = deformationGradientF(8:9:end);
            F3_3 = deformationGradientF(9:9:end);
            t2 = F1_1.^2;
            t3 = F1_2.^2;
            t4 = F1_3.^2;
            t5 = F2_1.^2;
            t6 = F2_2.^2;
            t7 = F2_3.^2;
            t8 = F3_1.^2;
            t9 = F3_2.^2;
            t10 = F3_3.^2;
            perElementPsi = mesh3D.elV'.*((mesh3D.elMu'.*(t2+t3+t4+t5+t6+t7+t8+t9+t10-3.0))./2.0+(mesh3D.elLambda'.*((mesh3D.elMu'.*(3.0./4.0))./mesh3D.elLambda'-F1_1.*F2_2.*F3_3+F1_1.*F2_3.*F3_2+F1_2.*F2_1.*F3_3-F1_2.*F2_3.*F3_1-F1_3.*F2_1.*F3_2+F1_3.*F2_2.*F3_1+1.0).^2)./2.0-(mesh3D.elMu'.*log(t2+t3+t4+t5+t6+t7+t8+t9+t10+1.0))./2.0);

            psi = sum(perElementPsi);
        end
    end
                      
end

