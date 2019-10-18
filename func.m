classdef func
    methods(Static)
        % functions
        function R = calcR(theta)
            R1 = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
            R2 =[cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
            R3 = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
            R = R1*R2*R3;
        end
        
        function P = R2p(R)
            e0 = sqrt((trace(R) + 1)/4);
            e1 = (R(3,2)-R(2,3))/(4*e0);
            e2 = (R(1,3)-R(3,1))/(4*e0);
            e3 = (R(2,1)-R(1,2))/(4*e0);
            P = [e0; e1; e2; e3];
        end
        
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
        
        function DP1 = getDP1(i,j,f,ground,print)
            DP1.Phi = i.ab'*i.A'*j.A*j.ab - f.f; 
            DP1.Nu = f.fd;
            DP1.Gamma = -i.a'*j.Bpdab*j.Pd - j.a'*i.Bpdab*i.Pd - 2*i.ad'*j.ad + f.fdd;
            if isequal(ground,'false')
                DP1.Phi_r = zeros(1,6);
                DP1.Phi_p = [j.a'*i.Bpab, i.a'*j.Bpab];
            else
                DP1.Phi_r = zeros(1,3);
                DP1.Phi_p = [j.a'*i.Bpab];
            end
            if isequal(print,'true')
                fprintf('-----------DP1-----------\n')
                fprintf('Phi is %f\n', DP1.Phi)
                fprintf('Nu is %f\n', DP1.Nu)
                fprintf('Gamma is %f\n', DP1.Gamma)
                fprintf('Partial derivatives for Phi_r and Phi_p are [ ');
                fprintf('%g ',DP1.Phi_r);
                fprintf('] and [ ');
                fprintf('%g ',DP1.Phi_p);
                fprintf(']\n');
            end
        end
        
        function DP2 = getDP2(i,j,ij,f,ground,print)
            DP2.Phi = i.ab'*i.A'*ij.d - f.f;
            DP2.Nu = f.fd;
            DP2.Gamma = -i.a'*j.Bpdsbq*j.Pd + i.a'*i.Bpdsbp*i.Pd - ij.d'*i.Bpdab*i.Pd - 2*i.ad'*ij.dd + f.fdd; 
            if isequal(ground,'false')
                DP2.Phi_r = [-i.a , i.a];
                DP2.Phi_p = [-i.a'*i.Bpsbp+ij.d'*i.Bpab, i.a'*j.Bpsbq]; 
            else
               DP2.Phi_r = [-i.a];
               DP2.Phi_p = [-i.a'*i.Bpsbp+ij.d'*i.Bpab];  
            end
            if isequal(print,'true')
                fprintf('-----------DP2-----------\n')
                fprintf('Phi is %f\n', DP2.Phi)
                fprintf('Nu is %f\n', DP2.Nu)
                fprintf('Gamma %f\n', DP2.Gamma)
                fprintf('Partial derivatives for DP2 Phi_r and Phi_p are [ ');
                fprintf('%g ',DP2.Phi_r);
                fprintf('] and [ ');
                fprintf('%g ',DP2.Phi_p);
                fprintf(']\n');
            end        
        end
        
        function CD = getCD(i,j,ij,f,ground,print)
            CD.Phi = ij.c'*ij.d - f.f;
            CD.Nu = f.fd;
            CD.Gamma = ij.c'*i.Bpdsbp*i.Pd - ij.c'*j.Bpdsbq*j.Pd + f.fdd;
            if isequal(ground,'false')
                CD.Phi_r = [-ij.c', ij.c'];
                CD.Phi_p = [-ij.c'*i.Bpsbp, ij.c'*j.Bpsbq];
            else
                CD.Phi_r = [-ij.c'];
                CD.Phi_p = [-ij.c'*i.Bpsbp];
            end
            if isequal(print,'true')
                fprintf('-----------CD-----------\n')
                fprintf('Phi is %f\n', CD.Phi)
                fprintf('Nu is %f\n', f.fd)
                fprintf('Gamma %f\n', CD.Gamma)
                fprintf('Partial derivatives for CD Phi_r and Phi_p are [ ');
                fprintf('%g ',CD.Phi_r);
                fprintf('] and [ ');
                fprintf('%g ',CD.Phi_p);
                fprintf(']\n');
            end    
        end
        
        function pnorm = getpnorm(i,j,ij,f,ground,print)
            pnorm.Phi = 0.5*i.P'*i.P - 0.5;
            pnorm.Nu = 0;
            pnorm.Gamma = -2*i.Pd'*i.Pd;
            pnorm.Phi_r = zeros(1,3);
            pnorm.Phi_p = 2*i.P';
        end
        
        function D = getD(i,j,ij,f,ground,print)
            D.Phi = ij.d'*ij.d - f.f;
            D.Nu = f.fd;
            D.Gamma = 2*ij.d'*i.Bpdsbp*i.Pd - 2*ij.d'*j.Bpdsbq*j.Pd - 2*ij.dd'*ij.dd + f.fdd;
            if isequal(ground,'false')
                D.Phi_r = [-2*ij.d', 2*ij.d'];
                D.Phi_p = [-2*ij.d'*i.Bpsbp, 2*ij.d'*j.Bpsbq];
            else
                D.Phi_r = [-2*ij.d'];
                D.Phi_p = [-2*ij.d'*i.Bpsbp];
            end
            if isequal(print,'true')
                fprintf('-----------D-----------\n')
                fprintf('Phi is %f\n', D.Phi)
                fprintf('Nu is %f\n', f.fd)
                fprintf('Gamma %f\n', D.Gamma) 
                fprintf('Partial derivatives for D Phi_r and Phi_p are [ ');
                fprintf('%g ',D.Phi_r);
                fprintf('] and [ ');
                fprintf('%g ',D.Phi_p);
                fprintf(']\n');
            end    
        end
    end
end