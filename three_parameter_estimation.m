function results = three_parameter_estimation(initParams) 
% Estimation of the parameters s, km4 and km5 by non-linear least squares optimization
data = load("PosFeed_Expdata");
tspan = data.tspan;
exp = data.exp;
 
opts = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxIter', 150, 'Diagnostics', 'off', 'Display', 'iter');
initParams=[0.8;0.05;0.05];

loBound = [0.1 0.05 0.05]; 
upBound = [0.8 0.1 0.1];

lsqnonlin(@residual,  log10 (initParams), log10 (loBound), log10 (upBound), opts);

function R = residual(b)       
a = 10.^b;
[T,Y]= reactionsolve(a);
residual = zeros([length(tspan), 2]);
for i = 1:length(T)
    for ii = 1:length(tspan)
        if T(i) == tspan(ii)
            residual(ii, 1) = exp(ii, 1)-Y(i, 1);
            residual(ii, 2) = exp(ii, 2)-Y(i, 2);
        end
    end
end
R = residual(:);
SSE = sum(R(:).^2)
results=a;
   
subplot(2,1,1);
plot(T,Y(:,1));
hold on;
grid minor
plot(tspan(:,1),exp(:,1),'*');
xlabel('Time');
ylabel('Act');
hold off;
str = "Parameters: " + sprintf('%g | ', a);
title(str); 
    
subplot(2,1,2);
plot(T,Y(:,2));
hold on;
grid minor
plot(tspan(:,1),exp(:,2),'*');
xlabel('Time');
ylabel('yP');
hold off;
      
hold off;
drawnow;
end

function [T,Y] = reactionsolve(a)
s=a(1);
km4=a(2);
km5=a(3);

% exp(tspan = 0)
x0 = [0.8; 0.6];
t = 0:0.01:12;
[T,Y] = ode45(@reaction, t, x0, []);

function dx = reaction(t,x)
Act=x(1);
yp=x(2);

par.ytot=1;
par.E=0.5;

% Values of fixed parameters
par.k1 = 1;
par.k2 = 0.8;
par.k3 = 1.2;
par.k4 = 1;
par.k5 = 1;

Act_dot = par.k1*s+par.k2*yp-par.k3*Act;
yp_dot = par.k4*Act*(par.ytot-yp)/(km4+par.ytot-yp)-par.k5*par.E*yp/(km5+yp);

dx = zeros(2,1);
dx(1)=Act_dot;
dx(2)=yp_dot;   
end
end
end


