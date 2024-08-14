function [k] = materialToBendingStiffness(E,nu,thickness)
%[k] = materialToBendingStiffness(E,nu,thickness)
%  E: young's modulus
%  nu: poisson ratio
%  thickness: thickness of the shell
 k = (E*thickness^3)/(24*(1-nu^2));
end

