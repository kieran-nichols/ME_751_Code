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
%             e1 = sqrt(0.25 * (-trace(R) + 1 + 2* R(1,1)));
%             e2 = sqrt(0.25 * (-trace(R) + 1 + 2* R(2,2)));
%             e3 = sqrt(0.25 * (-trace(R) + 1 + 2* R(3,3)));
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
        
        function DP1 = getDP1(i_ab,i_P,i_Pd,j_ab,j_P,j_Pd,f,ground,print)            
            i_A = func.calcA(i_P);
            i_a = i_A*i_ab;            
            j_A = func.calcA(j_P);
            j_a = j_A*j_ab; 
            i_Bpab = func.calcB(i_P,i_ab);
            j_Bpab = func.calcB(j_P,j_ab);
            i_Bpdab = func.calcB(i_Pd,i_ab); 
            j_Bpdab = func.calcB(j_Pd,j_ab);
            i_ad = i_Bpab*i_Pd;
            j_ad = j_Bpab*j_Pd;

            DP1.Phi = i_a'*j_a - f.f; 
            DP1.Nu = f.fd;
            DP1.Gamma = -i_a'*j_Bpdab*j_Pd - j_a'*i_Bpdab*i_Pd - 2*i_ad'*j_ad + f.fdd;
            if isequal(ground,'false')
                DP1.Phi_r = zeros(1,6);
                DP1.Phi_p = [j_a'*i_Bpab, i_a'*j_Bpab];
            else
                DP1.Phi_r = zeros(1,3);
                DP1.Phi_p = [j_a'*i_Bpab];
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
        
        function DP2 = getDP2(i_ab,i_P,i_Pd,i_sbp,j_ab,j_P,j_Pd, j_sbq,f,ground,print) 
            % Solve for all of the B matrices 
            ij_d = j_r + j_A*j_sbq - i_r - i_A*i_sbp;
            i_Bpdab = func.calcB(i_Pd,i_ab); 
            j_Bpdab = func.calcB(j_Pd,j_ab); 
            i_Bpdsbp = func.calcB(i_Pd,i_sbp); 
            j_Bpdsbq = func.calcB(j_Pd,j_sbq); 
            i_Bpab = func.calcB(i_P,i_ab);
            j_Bpab = func.calcB(j_P,j_ab);
            i_Bpsbp = func.calcB(i_P,i_sbp);
            j_Bpsbq = func.calcB(j_P,j_sbq);

            % Solve for a dots and ij.d
            i_ad = i_Bpab*i_Pd;
            j_ad = j_Bpab*j_Pd;
            % i.add = i.Bpdab*i.pd + i.Bpab*i.pdd;            
            DP2.Phi = i_ab'*i_A'*ij_d - f.f;
            DP2.Nu = f_fd;
            DP2.Gamma = -i_a'*j_Bpdsbq*j_Pd + i_a'*i_Bpdsbp*i_Pd - ij_d'*i_Bpdab*i_Pd - 2*i_ad'*ij_dd + f.fdd; 
            if isequal(ground,'false')
                DP2.Phi_r = [-i_a , i_a];
                DP2.Phi_p = [-i_a'*i_Bpsbp+ij_d'*i_Bpab, i_a'*j_Bpsbq]; 
            else
               DP2.Phi_r = [-i_a];
               DP2.Phi_p = [-i_a'*i_Bpsbp+ij_d'*i_Bpab];  
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
        
        function CD = getCD(ij_c,i_r,i_P,i_Pd,i_sbp,j_r,j_P,j_Pd, j_sbq,f,ground,print)
            i_A = func.calcA(i_P);          
            j_A = func.calcA(j_P); 
            % Solve for all of the B matrices 
            ij_d = j_r + j_A*j_sbq - i_r - i_A*i_sbp;
            i_Bpdsbp = func.calcB(i_Pd,i_sbp); 
            j_Bpdsbq = func.calcB(j_Pd,j_sbq); 
            i_Bpsbp = func.calcB(i_P,i_sbp);
            j_Bpsbq = func.calcB(j_P,j_sbq);          

            CD.Phi = ij_c'*ij_d - f.f;
            CD.Nu = f.fd;
            CD.Gamma = ij_c'*i_Bpdsbp*i_Pd - ij_c'*j_Bpdsbq*j_Pd + f.fdd;
            if isequal(ground,'false')
                CD.Phi_r = [-ij_c', ij_c'];
                CD.Phi_p = [-ij_c'*i_Bpsbp, ij_c'*j_Bpsbq];
            else
                CD.Phi_r = [-ij_c'];
                CD.Phi_p = [-ij_c'*i_Bpsbp];
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
        
        function pnorm = getpnorm(i_P,i_Pd)
            pnorm.Phi = 0.5*i_P'*i_P - 0.5;
            pnorm.Nu = 0;
            pnorm.Gamma = -i_Pd'*i_Pd;
            pnorm.Phi_r = zeros(1,3);
            pnorm.Phi_p = i_P';
        end
        
        function D = getD(i,j,ij,f,ground,print) 
                        % Solve for all of the B matrices
            j.sbq = [0,0,0]'; 
            ij.d = j.r + j.A*j.sbq - i.r - i.A*i.sbp;
            i.Bpdab = func.calcB(i.Pd,i.ab); 
            j.Bpdab = func.calcB(j.Pd,j.ab); 
            i.Bpdsbp = func.calcB(i.Pd,i.sbp); 
            j.Bpdsbq = func.calcB(j.Pd,j.sbq); 
            i.Bpab = func.calcB(i.P,i.ab);
            j.Bpab = func.calcB(j.P,j.ab);
            i.Bpsbp = func.calcB(i.P,i.sbp);
            j.Bpsbq = func.calcB(j.P,j.sbq);
            i.Bpdab = func.calcB(i.Pd,i.ab); 
            j.Bpdab = func.calcB(j.Pd,j.ab); 
            i.Bpab = func.calcB(i.P,i.ab);
            j.Bpab = func.calcB(j.P,j.ab);

            % Solve for a dots and ij.d
            i.ad = i.Bpab*i.Pd;
            j.ad = j.Bpab*j.Pd;
            % i.add = i.Bpdab*i.pd + i.Bpab*i.pdd;            

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