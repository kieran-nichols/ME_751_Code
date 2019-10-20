% Solves the IVP for HW7 prob 2
% yExact is the analytical solution
% yExact = 1/t + 1/t^2*tan(1/t+pi-1)

% ybe is the solution obtained with Backward Euler
%

h = 0.01; % Step size
tend = 20;
tm=0:h:tend;
itr = length(tm);
x = zeros(1,itr);
y = zeros(1,itr);
% yd = -y^2 - 1/t^4; % ydot

% where y(1) = 1; t from 1 to 10s
% Initial vaules
y(1) = 1;

% x and y variable for Backward Euler Method
% xbe = 1.1; % zeros(1,itr); 
% ybe = 2.1; % zeros(1,itr);

% Backward Euler equation
% x(n) = x(n-1) + h*xd(n);
% y(n) = y(n-1) + h*yd(n);

% Based of this IVP
% y(n) = x(1 - y/(1 + x^2));

% x(n) = x(n-1) + h(-x(n) - 4*x(n)*y(n)/(1 + x(n)^2));
% y(n) = y(n-1) + h*x(n)*(1 - y/(1 + x(n)^2));

% which is equivalent to g(x(n),y(n)) = [0;0];

% Iterate through time using Backward Euler and Jacobian
epi = 10^-3;
norm_corr = 10; %dummy high number

for n = 2:itr
%     xbe = x(n-1);
    ybe = y(n-1);
    m = 1;
    while norm_corr > epi  
        % Non linear system
        G = [xbe*(1+h) + 4*h*xbe*ybe/(1+xbe^2) - x(n-1);
            -h*xbe + ybe + h*xbe*ybe/(1+xbe^2) - y(n-1)];
        % Jacobian
        Jxy = [1 + h + 4*h*ybe*(1-xbe^2)/(1+xbe^2)^2, 4*h*xbe/(1+xbe^2);
          -h + h*(1-xbe^2)/(1+xbe^2)^2, 1 + h*xbe/(1+xbe^2)];
        % correction
        corr = Jxy\-G;
        xbe =  xbe + corr(1);
        ybe =  ybe + corr(2);
        norm_corr = norm(corr);
        % continue itr if condition is not met
        m = m + 1;
    end
    % redifined xn and yn with the last itr of while loop info
    x(n) = xbe;
    y(n) = ybe;  
end
 
% Plot calculated xn and yn vs time
figure
subplot(3,1,1)
plot(tm,x)
subplot(3,1,2)
plot(tm,y)