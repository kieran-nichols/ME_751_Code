%% Dynamics Engine
% Kieran Nichols
clc; clear all;

Case = 1;
if Case == 0
%%
% Initilizing variables
% using test case from some students on piazza
L = 2;
i.r = [0, 0, -L]'; 
offset = pi/2;
t = 0;
itr = 1000;
t_end = 10;
time = 0:t_end/itr:t_end;
itr = itr+1; %plus one to get the matrices to be the same amount of columns
torque = zeros(3,itr);
posp = zeros(3,itr);
velp = zeros(3,itr); % p represents prime reference frame aka point O'
accp = zeros(3,itr);

for k = 1:itr
    t = time(k);
    theta = -pi/4*cos(2*t)+offset;
    theta_d = pi/2*sin(2*t);
    theta_dd = pi*cos(2*t);

    f.f = cos(theta); 
    f.fd = ((pi*sin(2*t)*sin((pi*cos(2*t))/4 - pi/2))/2); 
    f.fdd = (pi*cos(2*t)*sin((pi*cos(2*t))/4 - pi/2) - (pi^2*sin(2*t)^2*cos((pi*cos(2*t))/4 - pi/2))/4); 

    ca = 0.05^2; % cross sectional area
    sa = 0.05*2*L; % side area
    m = 7800 * ca * 2 * L; % mass = density*volume
    I = [1/3*m*(0.05^2+4^2), 0,                 0;
        0,                  1/3*m*(0.05^2+4^2), 0;
        0,                  0,                 1/3*m*(0.05^2+0.05^2)];

    Mass = [m, 0, 0;
            0, m, 0;
            0, 0, m];
        
    Phi_q = [1 0 0 0;
             0 1 0 -L*cos(theta);
             0 0 1 L*cos(theta)];
         
    Qa = [0; 0; -m*9.81];
    
    qdd = f.fdd;
    
    Lamda = Phi_q'\(Qa - Mass*qdd); 
    torque(:,i) = Lambda;
    % Newton Euler equations
    % M*qdd + Phi_q'*Lamda = Qa
    % Lamda = Phi_q'\(Qa - M*qdd)
    % Lamda is the reaction force and torque that ensures that the contraints are fulfilled
    % qdd is the generalized accelerations
    % Phi_q is the Jacobian
    % Qa is the applied external force that do not have geometric contraints which will be gravity?
    % M is the Mass Matrix
    


figure
subplot(3,1,1)
plot(time,posp(1,:),time,posp(2,:),time,posp(3,:))
title('Position of point O-prime')
ylabel('position (m)')
legend('X','Y','Z')
subplot(3,1,2)
plot(time,velp(1,:),time,velp(2,:),time,velp(3,:))
title('Velocity of point O-prime')
ylabel('velocity (m/s)')
subplot(3,1,3)
plot(time,accp(1,:),time,accp(2,:),time,accp(3,:))
title('Acceleration of point O-prime')
ylabel('acceleration (m/s/s)')
xlabel('time(s)')
% 
figure
subplot(3,1,1)
plot(time,torque(1,:))
title('Torquex')
ylabel('Tx (Nm)')
legend('X','Y','Z')
subplot(3,1,2)
plot(time,torque(2,:))
title('Torque y')
ylabel('Ty (Nm)')
subplot(3,1,3)
plot(time,torque(3,:))
title('Torque_z')
ylabel('Tz (Nm)')
xlabel('time(s)')
end
else
%%
% Position, Vel, Acceleration Analysis
clc; clear all;

%%
% Initilizing variables
% using test case from some students on piazza

% i.r = [0, 0, -2]'; 
offset = pi/2;
t = 0;
itr = 1000; %1000
t_end = 10; %10
time = 0:t_end/itr:t_end;
itr = itr+1; %plus one to get the matrices to be the same amount of columns
posp = zeros(3,itr);
velp = zeros(3,itr); % p represents prime reference frame aka point O'
accelp = zeros(3,itr);
posq = zeros(3,itr); %q represents point q
velq = zeros(3,itr);
accelq = zeros(3,itr);
new_itr = 10;
posp_m = zeros(7,new_itr);
velp_m = zeros(7,new_itr);
accelp_m = zeros(7,new_itr);
L = 2;

% intial vel and accel
j.P = [1, 0, 0, 0]';
i.Pd = [0, 0, 0, 0]';
i.Pdd =[0;0;0;0];
j.Pd = [0, 0, 0, 0]';
j.Pdd =[0;0;0;0];
epi = 10^-1;
% i.r = L*[0 sin(pi/4) -cos(pi/4)]';
theta = -pi/4*cos(2*0)+offset;
%     theta_d = pi/2*sin(2*t);
%     theta_dd = pi*cos(2*t);
i.r = L*[0 cos(theta) -sin(theta)]';
%     i.R = func.calcR(theta);
i.ab = [1, 0, 0]';
i.sbp = [-L, 0, 0]'; 
j.r = [0, 0, 0]'; 
j.A = func.calcA(j.P);
j.sbq = [0, 0, 0]';
i.A = [0 0 1; cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0];
i.P = func.R2p(i.A);
% i.A = func.calcA(i.P);

for k = 1:itr
    m = 1; diff_norm = 10; % dummy norm
    while m<new_itr && diff_norm>epi
    t = time(k);
    theta = -pi/4*cos(2*t)+offset;
%     theta_d = pi/2*sin(2*t);
%     theta_dd = pi*cos(2*t);
    i.r = L*[0 cos(theta) -sin(theta)]';
%     i.R = func.calcR(theta);
    i.ab = [1, 0, 0]';
    i.sbp = [-L, 0, 0]'; 
    j.r = [0, 0, 0]'; 
    j.A = func.calcA(j.P);
    j.sbq = [0, 0, 0]';
%     j.a = j.A*j.ab;
    % ij.dd = j.rd + j.Bpsbq*j.Pd - i.rd - i.Bpsbp*i.Pd;

    f.f = cos(theta); 
    f.fd = ((pi*sin(2*t)*sin((pi*cos(2*t))/4 - pi/2))/2); 
    f.fdd = (pi*cos(2*t)*sin((pi*cos(2*t))/4 - pi/2) - (pi^2*sin(2*t)^2*cos((pi*cos(2*t))/4 - pi/2))/4); 

%     f.f = sin(pi/4. * cos(2* t)); 
%     f.fd = -pi/2. * sin(2* t) * cos(pi/4.* cos(2* t)); 
%     f.fdd = -pi/4. * (4* cos(2* t) * cos(pi/4.* cos(2* t)) + pi *sin(2.* t)*sin(2* t) *sin(pi/4.*cos(2* t))); 
    % Solve for DP1 and Rev constraints
    j.ab = [0, 1, 0]'; 
    Drive_DP1 = func.getDP1(i.ab,i.P,i.Pd,j.ab,j.P,j.Pd,f,'true','false');
    % Geo constraints are not affected by driving function
    f.f = 0; f.fd = 0; f.fdd = 0;
    i.ab = [0, 1, 0]';
    Geo_DP1_1 = func.getDP1(i.ab,i.P,i.Pd,j.ab,j.P,j.Pd,f,'true','false');
    j.ab = [0, 0, 1]'; % rotate about z
%     j.a = j.A*j.ab;
    Geo_DP1_2 = func.getDP1(i.ab,i.P,i.Pd,j.ab,j.P,j.Pd,f,'true','false');
%     j.ab = [1, 0, 0]'; % rotate about x
    ij.c = [1, 0, 0]';
    Geo_CD1 = func.getCD(ij.c,i.r,i.P,i.Pd,i.sbp,j.r,j.P,j.Pd, j.sbq,f,'true','false');
    ij.c = [0, 1, 0]';
    Geo_CD2 = func.getCD(ij.c,i.r,i.P,i.Pd,i.sbp,j.r,j.P,j.Pd, j.sbq,f,'true','false');
    ij.c = [0, 0, 1]';
    Geo_CD3 = func.getCD(ij.c,i.r,i.P,i.Pd,i.sbp,j.r,j.P,j.Pd, j.sbq,f,'true','false');
    Pnorm = func.getpnorm(i.P,i.Pd);

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

%         posq(:,k) = [0;0;0];
%         velq(:,k) = [0;0;0];
%         accq(:,k) = [0;0;0];

    %using Newton Rhapson
    % initial guess
%     posp_m(1:3,1) = L*[0; sin(f.f); -cos(f.f)];
%     velp_m(1:3,1) = L*[0.01;0.01;0.01]; %[0; sin(f.fd); -cos(f.fd)];
%     accelp_m(1:3,1) = L*[0.01;0.01;0.01]; %L*[0; sin(f.fdd); -cos(f.fdd)];
%     for m=1:new_itr-1
        posp_m(:,m+1) = [posp_m(1:3,m);i.P] - J\J_Phi; 
        i.P = posp_m(4:7,m+1);
        velp_m(:,m+1) = J\J_Nu;
        i.Pd = velp_m(4:7,m+1);
        accelp_m(:,m+1) = J\J_Gamma;
        i.Pdd = accelp_m(4:7,m+1);
        diff = posp_m(:,m+1) - posp_m(:,m);
        diff_norm = norm(diff);
        m = m+1;
    end
    posp(:,k) = posp_m(1:3,m-1);
    velp(:,k) = velp_m(1:3,m-1);
    accelp(:,k) = accelp_m(1:3,m-1);

end

figure
subplot(3,1,1)
plot(time,posp(1,:),time,posp(2,:),time,posp(3,:))
title('Position of point O-prime')
ylabel('position (m)')
legend('X','Y','Z')
subplot(3,1,2)
plot(time,velp(1,:),time,velp(2,:),time,velp(3,:))
title('Velocity of point O-prime')
ylabel('velocity (m/s)')
subplot(3,1,3)
plot(time,accelp(1,:),time,accelp(2,:),time,accelp(3,:))
title('Acceleration of point O-prime')
ylabel('acceleration (m/s/s)')
xlabel('time(s)')
% 
% figure
% subplot(3,1,1)
% plot(time,posq(1,:),time,posq(2,:),time,posq(3,:))
% title('Position of point Q')
% ylabel('position (m)')
% legend('X','Y','Z')
% subplot(3,1,2)
% plot(time,velq(1,:),time,velq(2,:),time,velq(3,:))
% title('Velocity of point Q')
% ylabel('velocity (m/s)')
% subplot(3,1,3)
% plot(time,accq(1,:),time,accq(2,:),time,accq(3,:))
% title('Acceleration of point Q')
% ylabel('acceleration (m/s/s)')
% xlabel('time(s)')

end



