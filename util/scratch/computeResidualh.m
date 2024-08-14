function residual = computeResidualh(constraintC,alphaTilde,lambda)
    residual = constraintC + alphaTilde*lambda;
end