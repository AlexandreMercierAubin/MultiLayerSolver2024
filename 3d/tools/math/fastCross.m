function [result]=fastCross(V1,V2)
	%fastCross calculates the cross product of two vectors without matlab's heavy validations.
	i=(V1(2)*V2(3) - V2(2)*V1(3));
	j=(V1(3)*V2(1) - V2(3)*V1(1));
	k=(V1(1)*V2(2) - V2(1)*V1(2));
	result=[i;j;k]; 
end