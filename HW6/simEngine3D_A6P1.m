%%
%     % system3D.m
%     % Kieran Nichols
clc; clear all;

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
t = 0;
f.f = 1.2;
f.fd = 2.5; 
f.fdd = 0.2;

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

DP1 = func.getDP1(i,j,f,'false','true');
DP2 = func.getDP2(i,j,ij,f,'false','true');
CD = func.getCD(i,j,ij,f,'false','true');
D = func.getD(i,j,ij,f,'false','true');