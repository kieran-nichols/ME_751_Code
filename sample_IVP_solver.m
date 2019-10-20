% Solves the IVP y_dot= -100*y + sin(t) and y(0)=0
% yExactis the analytical solution
% yFEis the solution obtained with Forward Euler
% yBEis the solution obtained with Backward Euler
%
% Input: the integration step-size
h = input('Input step size h:');
tend = 8;
tm=0:h:tend;
yExact= 1/10001*(100*sin(tm)-cos(tm)+exp(-100*tm));
yFE= zeros(size(tm'));
yBE= zeros(size(tm'));
for i=2:1:length(yFE)
yFE(i) = yFE(i-1)*(1-100*h) + h*sin(tm(i-1));
end
dummyINV= 1/(1+100*h);
for i=2:1:length(yBE)
yBE(i) = yBE(i-1)*dummyINV+ h*dummyINV*sin(tm(i));
end

figure
subplot(3,1,1)
plot(tm,yFE)
subplot(3,1,2)
plot(tm,yBE)