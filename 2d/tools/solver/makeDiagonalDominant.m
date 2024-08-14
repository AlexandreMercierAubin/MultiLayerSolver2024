function [A, b] = makeDiagonalDominant(A, b)
    n = size(b, 1);
    
    unfixedRows = 1:n;
    unfixedCols = 1:n;
    
    while size(unfixedRows, 2) * size(unfixedCols, 2) > 0
        
        maxElem = 0;
        row = 0;
        col = 0;
        
        for i = unfixedRows
            for j = unfixedRows
                if row == 0 || abs(A(i, j)) > maxElem
                   maxElem = abs(A(i, j));
                   row = i;
                   col = j;
                end
            end
        end
        
        swap = A(row, :);
        A(row, :) = A(col, :);
        A(col, :) = swap;
        unfixedRows = unfixedRows(unfixedRows ~= row);
        unfixedCols = unfixedCols(unfixedCols ~= col);
        
        swap = b(row);
        b(row) = b(col);
        b(col) = swap;
    end
end

