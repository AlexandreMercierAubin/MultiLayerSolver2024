function [result]=fastVectorNorm(V1)
	%fastVectorNorm norm without matlab's heavy checks assumes column
    %vectors
    result = sqrt(V1'*V1);
end