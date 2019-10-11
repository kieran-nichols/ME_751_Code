%%
%     % system3D.m
%     % Kieran Nichols
clc; clear all;

%%
% % Initilizing variables
% nbodies = 2;
% % body i
% i.r = [0;0;0]; 
% i.p = [0;0;0];
% i.q = [-1;0;0];
% i.P = [1;0;0;0];
% i.A = func.calcA(i.P);
% i.ab = i.q - i.p; 
% i.a = i.A*i.ab;
% i.Pd = [0;0;0;0];
% i.Pdd =[0;0;0;0];
% i.rd = [0;0;0];
% i.rdd = [0;0;0];
% i.sbp = i.p; % ?
% 
% % body j
% j.r = [-1;0;0]; 
% j.p = [0;0;0];
% j.q = [-1;0;0];
% j.P = [1;0;0;0];
% j.A = func.calcA(j.P);
% j.ab = j.q - j.p; 
% j.a = j.A*j.ab;
% j.Pd = [0;0;0;0];
% j.Pdd =[0;0;0;0];
% j.rd = [0;0;0];
% j.rdd = [0;0;0];
% j.sbq = j.q; % ?
% 
% c = i.p - j.q; % diff between p coordinates on i and q coordinates on j
% 
% f.f = [0]; % ?
% f.fd = [0]; % ?
% f.fdd = [0]; % ?

%%
% Initilizing variables
% using test case from some students on piazza
nbodies = 2;
% body i
i.r = [8, 6, -3]'; 
% i.p = [0;0;0];
% i.q = [-1;0;0];
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

c = [0.3, 0.4, -6]'; % diff between p coordinates on i and q coordinates on j

f.f = 1.2;
f.fd = 2.5; 
f.fdd = 0.2; 
%%
nu = f.fd; % RHS velocity is equal to diffrential of function

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
Phi.DP1 = i.ab'*i.A'*j.A*j.ab - f.f; 
Gamma.DP1= -i.a'*j.Bpdab*j.Pd - j.a'*i.Bpdab*i.Pd - 2*i.ad'*j.ad + f.fdd; 

Phi.DP2 = i.ab'*i.A'*ij.d - f.f; 
Gamma.DP2= -i.a'*j.Bpdsbq*j.Pd + i.a'*i.Bpdsbp*i.Pd - ij.d'*i.Bpdab*i.Pd - 2*i.ad'*ij.dd + f.fdd; 

Phi.CD = c'*ij.d - f.f;
Gamma.CD = c'*i.Bpdsbp*i.Pd - c'*j.Bpdsbq*j.Pd + f.fdd;

Phi.D = ij.d'*ij.d - f.f;
Gamma.D = 2*ij.d'*i.Bpdsbp*i.Pd - 2*ij.d'*j.Bpdsbq*j.Pd - 2*ij.dd'*ij.dd + f.fdd;

% fprintf('The values of Phi.DP1 and Phi.CD are %f and %f\n', Phi.DP1,Phi.CD)
% fprintf('The right hand side of velocity equation is %f\n', nu)
% fprintf('The right hand side of acceleration equation contains coefficients of %f and %f\n', Gamma.DP1, Gamma.CD)

fprintf('The values of Phi.DP2 and Phi.D are %f and %f\n', Phi.DP2,Phi.D)
fprintf('The right hand side of velocity equation is %f\n', nu)
fprintf('The right hand side of acceleration equation contains coefficients of %f and %f\n', Gamma.DP2, Gamma.D)
%%
% Phi_r are the jacobian wrt partial derivatives of r
% Phi_p are the jacobian wrt partial derivatives of p
DP1.Phi_r = zeros(1,6);
DP1.Phi_p = [j.a'*i.Bpab, i.a'*j.Bpab];

DP2.Phi_r = [-i.a , i.a];
DP2.Phi_p = [-i.a'*i.Bpsbp+ij.d'*i.Bpab, i.a'*j.Bpsbq]; %probably wrong

CD.Phi_r = [-c' c'];
CD.Phi_p = [-c'*i.Bpsbp, c'*j.Bpsbq];

D.Phi_r = [-2*ij.d', 2*ij.d'];
D.Phi_p = [-2*ij.d'*i.Bpsbp, 2*ij.d'*j.Bpsbq];

% fprintf('\n');
% fprintf('The partial derivatives for DP1 Phi_r are [ ');
% fprintf('%g ',DP1.Phi_r);
% fprintf(']\n');
% 
% fprintf('The partial derivatives for DP1 Phi_p are [ ');
% fprintf('%g ',DP1.Phi_p);
% fprintf(']\n');
% 
% fprintf('The partial derivatives for CD Phi_r are [ ');
% fprintf('%g ',CD.Phi_r);
% fprintf(']\n');
% 
% fprintf('The partial derivatives for CD Phi_p are [ ');
% fprintf('%g ',CD.Phi_p);
% fprintf(']\n');

fprintf('\n');
fprintf('The partial derivatives for DP2 Phi_r are [ ');
fprintf('%g ',DP2.Phi_r);
fprintf(']\n');

fprintf('The partial derivatives for DP2 Phi_p are [ ');
fprintf('%g ',DP2.Phi_p);
fprintf(']\n');

fprintf('The partial derivatives for D Phi_r are [ ');
fprintf('%g ',D.Phi_r);
fprintf(']\n');

fprintf('The partial derivatives for D Phi_p are [ ');
fprintf('%g ',D.Phi_p);
fprintf(']\n');

% disp('complete');

