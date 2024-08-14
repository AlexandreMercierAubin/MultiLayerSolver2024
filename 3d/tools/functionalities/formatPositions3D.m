function x = formatPositions3D(p)
    %Get x in a n x 3 format compatible with most gptoolbox
    %functions
    x = reshape(p',3,[])';
end