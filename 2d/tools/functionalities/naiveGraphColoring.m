function elementColor = naiveGraphColoring(T)
    elementColor = zeros(size(T,1),1);
    %elements*color compute time... just assigns the first color
    %that works

    vertexUsedColors = dictionary("",false);

    for i = 1:numel(elementColor)
        verts = T(i,:);
        colorInd = 1;
        keyStrings = string(verts) +","+string(colorInd);
        while any(isKey(vertexUsedColors,keyStrings))
            colorInd = colorInd + 1;
            keyStrings = string(verts) +","+string(colorInd);
        end
        vertexUsedColors = vertexUsedColors.insert(keyStrings,true(size(keyStrings)));
        elementColor(i) = colorInd;
    end
end

