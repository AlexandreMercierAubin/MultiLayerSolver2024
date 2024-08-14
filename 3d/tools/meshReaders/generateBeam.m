function [V,T] = generateBeam(n, m, l)
    %Generates a mesh according to dimensions n, m, l
    
    T = zeros(n * m * l * 6, 4);
    V = zeros((n + 1) * (m + 1) * (l + 1), 2);

    pinned = zeros((l + 1) * (m + 1), 1);
    currentTet = 1;
    
    for k = 1:(l + 1)
        for i = 1:(n + 1)
            for j = 1:(m + 1)
                pos = (k - 1) * (n + 1) * (m + 1) + (i - 1) * (m + 1) + j;
                V(pos, 1) = j - 1;
                V(pos, 2) = i - 1;
                V(pos, 3) = k - 1;

                if i == 1
                    pinned((k - 1) * (m + 1) + j) = pos;
                end
                
                if k <= l && j <= m && i <= n
                   v1 = (k - 1) * (n + 1) * (m + 1) + (i - 1) * (m + 1) + j;
                   v2 = v1 + 1;
                   v3 = (k - 1) * (n + 1) * (m + 1) + i * (m + 1) + j;
                   v4 = v3 + 1;
                   v5  = k * (n + 1) * (m + 1) + (i - 1) * (m + 1) + j;
                   v6 = v5 + 1;
                   v7 = k * (n + 1) * (m + 1) + i * (m + 1) + j;
                   v8 = v7 + 1;
                   
                   T(currentTet, :) = [v1, v2, v3, v6];
                   T(currentTet + 1, :) = [v2, v3, v4, v6];
                   T(currentTet + 2, :) = [v3, v5, v6, v7];
                   T(currentTet + 3, :) = [v3, v6, v7, v8];
                   T(currentTet + 4, :) = [v1, v3, v5, v6];
                   T(currentTet + 5, :) = [v3, v4, v6, v8];
                   currentTet = currentTet + 6;
                end
            end
        end
    end
end