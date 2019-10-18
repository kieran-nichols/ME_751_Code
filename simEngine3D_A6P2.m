%%
%     % system3D.m
%     % Kieran Nichols
clc; clear all;
%%
i.r = [0, 0, -2]'; 
offset = pi/2;
t = 0;
itr = 1000*10;
t_end = 10;
time = 0:t_end/itr:t_end;
itr = itr+1; %plus one to get the matrices to be the same amount of columns
posp = zeros(3,itr);
velp = zeros(3,itr); % p represents prime reference frame aka point O'
accp = zeros(3,itr);
posq = zeros(3,itr); %q represents point q
velq = zeros(3,itr);
accq = zeros(3,itr);
new_itr = 10;
posp_m = zeros(7,new_itr);
velp_m = zeros(7,new_itr);
accelp_m = zeros(7,new_itr);
L = 2;

Case = '1';

for k = 1:1
    t = time(k);
    theta = -pi/4*cos(2*t)+offset;
    theta_d = pi/2*sin(2*t);
    theta_dd = pi*cos(2*t);
    i.R = func.calcR(theta);
    i.P = func.R2p(i.R);
    i.A = func.calcA(i.P);
    i.ab = [0, 1, 0]'; 
    i.a = i.A*i.ab;
    i.bb = [0, 0, 1]'; 
    i.b = i.A*i.ab;
    i.Pd = [0, 0, 0, 0]';
    i.Pdd =[0;0;0;0];
    i.sbp = [0.2, -1.0, 1.5]'; 
    j.r = [0, 0, 0]'; 
    j.P = [0, 0, 0, 0]';
    j.A = func.calcA(j.P);
    j.ab = [0, 1, 0]'; 
    j.a = j.A*j.ab;
    j.Pd = [0, 0, 0, 0]';
    j.Pdd =[0;0;0;0];
    % ij.dd = j.rd + j.Bpsbq*j.Pd - i.rd - i.Bpsbp*i.Pd;

    % Solve for all of the B matrices
    j.sbq = [-2,0,0]'; 
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
    f.f = cos(theta); 
    f.fd = ((pi*sin(2*t)*sin((pi*cos(2*t))/4 - pi/2))/2); 
    f.fdd = (pi*cos(2*t)*sin((pi*cos(2*t))/4 - pi/2) - (pi^2*sin(2*t)^2*cos((pi*cos(2*t))/4 - pi/2))/4); 

    % Test case 1
    if isequal(Case,'1')
        % Solve for DP1 and Rev constraints
        Drive_DP1 = func.getDP1(i,j,f,'true','false');
        Geo_DP1_1 = func.getDP1(i,j,f,'true','false');
        Geo_DP1_2 = func.getDP1(i,j,f,'true','false');
        ij.c = [1, 0, 0]';
        Geo_CD1 = func.getCD(i,j,ij,f,'true','false');
        ij.c = [0, 1, 0]';
        Geo_CD2 = func.getCD(i,j,ij,f,'true','false');
        ij.c = [0, 0, 1]';
        Geo_CD3 = func.getCD(i,j,ij,f,'true','false');
        Pnorm = func.getpnorm(i,j,ij,f,'true','false');

        J = [Drive_DP1.Phi_r Drive_DP1.Phi_p;
           Geo_DP1_1.Phi_r Geo_DP1_1.Phi_p;
           Geo_DP1_2.Phi_r Geo_DP1_2.Phi_p;
           Geo_CD1.Phi_r Geo_CD1.Phi_p;
           Geo_CD2.Phi_r Geo_CD2.Phi_p;
           Geo_CD3.Phi_r Geo_CD3.Phi_p;
           Pnorm.Phi_r Pnorm.Phi_p];
       
       J_Gamma = [Drive_DP1.Gamma;
           Geo_DP1_1.Gamma;
           Geo_DP1_2.Gamma;
           Geo_CD1.Gamma;
           Geo_CD2.Gamma;
           Geo_CD3.Gamma;
           Pnorm.Gamma];

       J_Nu = [Drive_DP1.Nu;
           Geo_DP1_1.Nu;
           Geo_DP1_2.Nu;
           Geo_CD1.Nu;
           Geo_CD2.Nu;
           Geo_CD3.Nu;
           Pnorm.Nu];

          J_Phi = [Drive_DP1.Phi;
           Geo_DP1_1.Phi;
           Geo_DP1_2.Phi;
           Geo_CD1.Phi;
           Geo_CD2.Phi;
           Geo_CD3.Phi;
           Pnorm.Phi];

            posq(:,k) = [0;0;0];
            velq(:,k) = [0;0;0];
            accq(:,k) = [0;0;0];

        %using Newton Rhapson
        % initial guess
        posp_m(1:3,1) = L*[0; sin(f.f); -cos(f.f)];
        velp_m(1:3,1) = L*[0.01;0.01;0.01]; %[0; sin(f.fd); -cos(f.fd)];
        accelp_m(1:3,1) = L*[0.01;0.01;0.01]; %L*[0; sin(f.fdd); -cos(f.fdd)];
        for m=1:new_itr-1
            posp_m(:,m+1) = [posp_m(1:3,m);i.P] - J\J_Phi; 
            i.P = posp_m(4:7,m);
            velp_m(:,m+1) = J\J_Nu;
%         velp_m(:,m+1) = [velp_m(1:3,m);i.Pd] - [DP1.Phi_r,DP1.Phi_p]\DP1.Nu;
            i.Pd = velp_m(4:7,m);
%         accelp_m(:,m+1) = [accelp_m(1:3,m);i.Pdd] - J\DP1.Gamma;
        accelp_m(:,m+1) = J\J_Gamma;
        i.Pdd = accelp_m(4:7,m);
        posp(:,k) = posp_m(1:3,new_itr);
        velp(:,k) = velp_m(1:3,new_itr);
        posp(:,k) = posp_m(1:3,new_itr);
        end
    else
        f.fx = 0;
        f.fy = cos(theta); 
        f.fz = -sin(theta);
        % function solve
        f.fdx = 0;
        f.fdy = -theta_d*sin(theta); 
        f.fdz = -theta_d*cos(theta);
        f.fddx = 0;
        f.fddy = -(theta_dd*sin(theta)+theta_d^2*cos(theta)); 
        f.fddz = -(theta_dd*cos(theta)-theta_d^2*sin(theta)); 
        % using funcitons
        posp(1:3,k) = L*[f.fx;f.fy;f.fz];
        velp(1:3,k) = L*[f.fdx;f.fdy;f.fdz];
        accp(1:3,k) = L*[f.fddx;f.fddy;f.fddz]; 
    end
    if k==1
        J_Phi
        J
        J_Nu
        J_Gamma
    end
end


% disp('complete');

