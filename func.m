classdef func
    methods(Static)
        % functions
        % check G
        function G = calcG(p)
        G = [-p(2) p(1) p(4) -p(3); 
             -p(3) -p(3) p(1) -p(2);
             -p(4) p(2) p(2) -p(1)];     
        end
        
        function at = tilde(a) 
            at = [0 -a(3) a(2);
                  a(3) 0 -a(1);
                  -a(2) a(1) 0];
        end
        
        function B = calcB(p,ab)
        B = 2*[(p(1)*eye(3)+func.tilde(p(2:end)))*ab,... 
            p(2:end)*ab'-(p(1)*eye(3)+func.tilde(p(2:end)))*func.tilde(ab)];     
        end

        function E = calcE(p)
        E = [-p(2) p(1) -p(4) p(3); 
             -p(3) p(3) p(1) -p(2);
             -p(4) -p(2) p(2) p(1)];     
        end

        function A = calcA(p)
        A = 2*[p(1)^2+p(2)^2-0.5 p(2)*p(3)-p(1)*p(4) p(2)*p(4)+p(1)*p(3);
             p(2)*p(3)+p(1)*p(4) p(1)^2+p(3)^2-0.5 p(3)*p(4)-p(1)*p(2);
             p(2)*p(4)-p(1)*p(3) p(3)*p(4)+p(1)*p(2) p(1)^2+p(4)^2-0.5];     
        end
    end
end