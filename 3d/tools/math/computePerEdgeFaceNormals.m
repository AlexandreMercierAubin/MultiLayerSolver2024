function [normal,normalTilde] = computePerEdgeFaceNormals(e0,e2,eTilde2)
        normal = fastCross(e0',e2');
        n1norm = fastVectorNorm(normal);
        normal = (normal./n1norm)';

        normalTilde = fastCross(eTilde2',e0');
        n2norm = fastVectorNorm(normalTilde);
        normalTilde = (normalTilde./n2norm)';
end

