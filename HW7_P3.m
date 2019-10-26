clc; clear all;
% Solves the IVP for HW7 prob 3

% ybe is the solution obtained with Backward Euler
%

h = .001; % Step size
tend = 10;
t=1:h:tend;
itr = length(t);
y = zeros(1,itr);
yE = zeros(1,itr);

% yd = -y^2 - 1/t^4; % ydot

% where y(1) = 1; t from 1 to 10s
% Initial vaules
y(1) = 1;

% x and y variable for Backward Euler Method
% ybe = 2.1; % zeros(1,itr);

% Backward Euler equation
% y(n) = y(n-1) + h*yd(n);

% Based of this IVP
% y(n) - y(n-1) + h * (-y^2 - 1/t^4) = 0;

% which is equivalent to g(y(n),t(n)) = [0]; ?

% BDF
yb4 = 0; yb3 = 0; yb2 = 0; yb1 = 0; yb = 0;
ybdf = zeros(1,itr);

% Iterate through time using Backward Euler and Jacobian
epi = 10^-3;
norm_corr = 10; %dummy high number
norm_corr1 = 10;

for n = 2:itr
    % yE is the exact analytical solution
    yE(n) = 1/t(n) + (1/t(n)^2)*tan(1/t(n)+pi-1);
    tbe = t(n-1);
    ybe = y(n-1);    
    m = 1;
    % BDF method
    if n > 4
        yb = ybdf(n-4); yb1 = ybdf(n-3); yb2 = ybdf(n-2); yb3 = ybdf(n-1); yb4 = ybdf(n);
    end
    while norm_corr > epi && m < 100 && norm_corr1 > epi
        % Non linear system
        G = ybe - y(n-1) - h * (-ybe^2 - 1/tbe^4);
        % Jacobian
        Jxy = [-4*h/tbe^5, 1+2*h*ybe]; % [1-2*h*ybe]; ? -4*h/t(n)^5
        % correction
        corr = -Jxy\G; % Which order? what should Jxy look like?
%         tbe =  tbe + corr(1);
        ybe =  ybe + corr(2);
        norm_corr = norm(corr);
        % continue itr if condition is not met
        if n > 4
            G1 = yb4 - 1/25*(48*yb3 - 36*yb2 + 16*yb1 -3*yb + 12*h*(-yb4^2 - 1/tbe^4)); 
            Jxy1 = [1 - 24/25*h*yb];
            corr1 = -Jxy\G;
            yb4 = yb4 + corr(2);
            norm_corr1 = norm(corr1);
        end
        m = m + 1;
    end
%     m
    % redifined xn and yn with the last itr of while loop info
    y(n) = ybe;
    if n<5
        ybdf(n)= yE(n);
    else
    ybdf(n) = yb4;
    end
    
end
 
% Plot calculated xn and yn vs time
figure
% subplot(3,1,1)
plot(t,y,t,yE, t,ybdf)
ylim(2*[-2,2])
legend('ybdf','yExact','ybdf')
% subplot(3,1,2)
% plot(t,y)
% ylim(.1*[-2,2])