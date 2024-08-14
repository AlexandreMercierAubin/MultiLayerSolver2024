classdef BDF
    %BDF a bunch of BDF functions
    
    methods(Static = true)
        function v = BDF1(x2,x1,h)
            v = (x2-x1)./h;
        end

        function v = BDF2(x3,x2,x1,h)
            v = (3/(2*h)).*(x3 - (4/3).*x2 + (1/3).*x1);
        end

        function v = BDF3(x4,x3,x2,x1,h)
            v = (11/(6*h)).*(x4 - (18/11).*x3 + (9/11).*x2 -(2/11).*x1);
        end
    end
end

