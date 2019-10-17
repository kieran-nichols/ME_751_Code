%%
%     % system3D.m
%     % Kieran Nichols
clc; clear all;

%%
% Initilizing variables
% using test case from some students on piazza
% Test case 1
% nbodies = 2;

% body i
i.r = [8, 6, -3]'; 
i.p = [0;0;0];
i.q = [-1;0;0];
i.P = [4, 3, -5, 1]';
i.P = i.P/norm(i.P);
i.A = func.calcA(i.P);
i.ab = [-1.2, 1 ,0.3]'; 
i.a = i.A*i.ab;
i.Pd = [-0.2, 1.3, 3.4, 0]';
i.Pd(4) = -dot(i.Pd,i.P)/i.P(4);
i.Pd = i.Pd/norm(i.Pd);
i.Pdd =[0;0;0;0];
i.rd = [7;8;9];
i.rdd = [0;0;0];
i.sbp = [0.1, -0.3, 6.0]'; % ?

% body j
j.r = [-0.5, 1.6, -6.3]'; 
j.P = [3.3, -4, 5.1, 6]';
j.P = j.P/norm(j.P,2);
j.A = func.calcA(j.P);
j.ab = [1.2, 4.5, 3.1]'; 
j.a = j.A*j.ab;
j.Pd = [0.6, -3.7, 5.1, 0]';
j.Pd(4) = -dot(j.Pd,j.P)/j.P(4);
j.Pd = j.Pd/norm(j.Pd,2);
j.Pdd =[0;0;0;0];
j.rd = [11;12;13];
j.rdd = [0;0;0];
j.sbq = [0.2, -1.0, 1.5]'; 

ij.c = [0.3, 0.4, -6]'; % diff between p coordinates on i and q coordinates on j

%%
% % % Test case 2
% % % nbodies = 2;
% % % body i
% i.r = [2, 0, 0]'; 
% % i.p = [0;0;0];
% % i.q = [-1;0;0];
% i.P = [0, 0, 0, 1]';
% % i.P = i.P/norm(i.P);
% i.A = func.calcA(i.P);
% i.ab = [2, 0, 0]'; 
% i.a = i.A*i.ab;
% i.Pd = [0, 0, 0, 0]';
% % i.Pd(4) = -dot(i.Pd,i.P)/i.P(4);
% % i.Pd = i.Pd/norm(i.Pd);
% i.Pdd =[0;0;0;0];
% i.rd = [0;0;0];
% i.rdd = [0;0;0];
% i.sbp = [2, 0, 0]'; 
% 
% % body j
% j.r = [0, 0, 0]'; 
% j.P = [0, 0, 0, 0]';
% % j.P = j.P/norm(j.P,2);
% j.A = func.calcA(j.P);
% j.ab = [0, 0, 0]'; 
% j.a = j.A*j.ab;
% j.Pd = [0, 0, 0, 0]';
% % j.Pd(4) = -dot(j.Pd,j.P)/j.P(4);
% % j.Pd = j.Pd/norm(j.Pd,2);
% j.Pdd =[0;0;0;0];
% j.rd = [0;0;0];
% j.rdd = [0;0;0];
% j.sbq = [0, 0, 0]'; 
% 
% % ij.c = [0.3, 0.4, -6]'; % diff between p coordinates on i and q coordinates on j
%%
t = 0;
itr = 1000;
t_end = 10;
time = 0.001:t_end/itr:t_end;
pos = zeros(3,itr);
vel = zeros(3,itr);
acc = zeros(3,itr);
L = 2;
offset = pi/2;

for k = 1:1 %itr
    t = time(k);
    theta = pi/4*cos(2*t);
    theta_d = -pi/2*sin(2*t);
    theta_dd = -pi*cos(2*t);
    f.f = 1.2;%cos(theta+offset); 
    f.fd = 2.5;%-theta_d*sin(theta+offset); 
    f.fdd = 0.2;%-theta_dd*sin(theta+offset)-theta_d^2*cos(theta+offset); 

    % Solve for all of the B matrices
    i.Bpdab = func.calcB(i.Pd,i.ab); 
    j.Bpdab = func.calcB(j.Pd,j.ab); 
    i.Bpdsbp = func.calcB(i.Pd,i.sbp); 
    j.Bpdsbq = func.calcB(j.Pd,j.sbq); 
    i.Bpab = func.calcB(i.P,i.ab);
    j.Bpab = func.calcB(j.P,j.ab);
    i.Bpsbp = func.calcB(i.P,i.sbp);
    j.Bpsbq = func.calcB(j.P,j.sbq);

    % Solve for a dots and ij.d
    i.ad = i.Bpab*i.Pd;
    j.ad = j.Bpab*j.Pd;
    % i.add = i.Bpdab*i.pd + i.Bpab*i.pdd;
    ij.d = j.r + j.A*j.sbq - i.r - i.A*i.sbp;
    ij.dd = j.rd + j.Bpsbq*j.Pd - i.rd - i.Bpsbp*i.Pd;

    % Solve for Phi and Gamma for DP1 and CD constraints
    DP1 = func.getDP1(i,j,f,'true','true');
    DP2 = func.getDP2(i,j,ij,f,'true','true');
    CD = func.getCD(i,j,ij,f,'true','true');
    D = func.getD(i,j,ij,f,'true','true');
    
    nu = f.fd; % RHS velocity is equal to diffrential of function ?
    
    pos(:,k) = L*[-cos(theta);sin(theta);0];
    vel(:,k) = L*[cos(f.fd);sin(f.fd);0];
    acc(:,k) = L*[cos(f.fdd);sin(f.fdd);0];     
end
% figure
% subplot(3,1,1)
% plot(time,pos(1,:),time,pos(2,:),time,pos(3,:))
% subplot(3,1,2)
% plot(time,vel(1,:),time,vel(2,:),time,vel(3,:))
% subplot(3,1,3)
% plot(time,acc(1,:),time,acc(2,:),time,acc(3,:))

% fprintf('The values of Phi.DP1 and Phi.CD are %f and %f\n', Phi.DP1,Phi.CD)
% fprintf('The values of Phi.DP1 is %f\n', Phi.DP1)
% fprintf('The right hand side of velocity equation is %f\n', nu)
% fprintf('The right hand side of acceleration equation contains coefficients of %f\n', Gamma.DP1)
% fprintf('The right hand side of acceleration equation contains coefficients of %f and %f\n', Gamma.DP1, Gamma.CD)
% fprintf('The values of Phi.DP2 and Phi.D are %f and %f\n', Phi.DP2,Phi.D)
% fprintf('The right hand side of velocity equation is %f\n', nu)
% fprintf('The right hand side of acceleration equation contains coefficients of %f and %f\n', Gamma.DP2, Gamma.D)

% disp('complete');

