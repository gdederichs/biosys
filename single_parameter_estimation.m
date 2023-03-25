function results =single_parameter_estimation(initParams) 
% Estimation of the parameters by non-linear least squares optimization

data = load("PosFeed_Expdata");
tspan = data.tspan;
exp = data.exp;
 
opts = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxIter', 150, 'Diagnostics', 'off', 'Display', 'iter');
initParams=[0.8];

% single parameter estimated s
loBound = [0.05]; 
upBound = [0.8];

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
% size(R)
% SSE_Act = R(1:4).*R(1:4);
% SSE_ACT = sum(SSE_Act(:))
% 
% SSE_yp = R(5:8).*R(5:8);
% SSE_yp = sum(SSE_Act(:))
results=a;
   
subplot(2,1,1);
plot(T,Y(:,1));
hold on;
grid minor
plot(tspan(:,1),exp(:,1),'*');
xlabel('Time');
ylabel('Act');
hold off;
str = "Parameter: "+string(a);
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

% exp(tspan = 0)
x0 = [0.8; 0.6];
t = 0:0.01:12;
[T,Y] = ode45(@reaction, t, x0, []);

function dx = reaction(t,x)

Act=x(1);
yp=x(2);

par.ytot=1;
par.E=0.5;


par.k1=1;
par.k2=0.8;
par.k3=1.2;
par.k4=1;
par.k5=1;
par.km4=0.05;
par.km5=0.05;


Act_dot = par.k1*s+par.k2*yp-par.k3*Act;
yp_dot = par.k4*Act*(par.ytot-yp)/(par.km4+par.ytot-yp)-par.k5*par.E*yp/(par.km5+yp);

dx = zeros(2,1);
dx(1)=Act_dot;
dx(2)=yp_dot;
    
end

end

end


