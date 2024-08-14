function [result]=fastDot(V1,V2)
	%fastDot calculates the dot product of two vectors without matlab's heavy validations.
	result=V1'*V2; 
end