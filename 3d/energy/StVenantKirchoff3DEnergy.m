classdef StVenantKirchoff3DEnergy < EnergyModel
    %StVenantKirchoff3DEnergy class that allows the computation of corotational
    %energy based on Saint-Venant Kirchoff energy formula
    properties
        
    end
    
    methods
        function obj = StVenantKirchoff3DEnergy()
            obj@EnergyModel();
            obj.name = "StVK";
        end 
        
        function computeEnergy(obj, mesh3D, deformationGradientF)
            [ii, jj, CblockVals, dpsidF, psi] = mexSTVK3D( deformationGradientF, mesh3D.elV, mesh3D.elMu, mesh3D.elLambda );
            sizeC = 9*size(mesh3D.t,1);
            C = sparse( ii, jj, CblockVals,sizeC,sizeC);
            obj.elasticForces = mesh3D.B' * dpsidF; 
            obj.elasticForces(mesh3D.pinnedDOFs) = 0;
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
            t11 = t2./2.0;
            t12 = t3./2.0;
            t13 = t4./2.0;
            t14 = t5./2.0;
            t15 = t6./2.0;
            t16 = t7./2.0;
            t17 = t8./2.0;
            t18 = t9./2.0;
            t19 = t10./2.0;
            perElementPsi = mesh3D.elV'.*(mesh3D.elMu'.*((t11+t14+t17-1.0./2.0).^2+(t12+t15+t18-1.0./2.0).^2+(t13+t16+t19-1.0./2.0).^2+((F1_1.*F1_2)./2.0+(F2_1.*F2_2)./2.0+(F3_1.*F3_2)./2.0).^2.*2.0+((F1_1.*F1_3)./2.0+(F2_1.*F2_3)./2.0+(F3_1.*F3_3)./2.0).^2.*2.0+((F1_2.*F1_3)./2.0+(F2_2.*F2_3)./2.0+(F3_2.*F3_3)./2.0).^2.*2.0)+(mesh3D.elLambda'.*(t11+t12+t13+t14+t15+t16+t17+t18+t19-3.0./2.0).^2)./2.0);
            psi = sum(perElementPsi);
        end
    end
                      
end

