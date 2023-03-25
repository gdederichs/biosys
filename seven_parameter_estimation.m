function results = seven_parameter_estimation(initParams) 

% Estimation of the parameters by non-linear least squares optimization
data = load("PosFeed_Expdata");
tspan = data.tspan;
exp = data.exp;

figure()
hold on
grid minor
plot(tspan,exp,'*')
title('True Data')
legend('Act','yp')
xlabel('time')
ylabel('Concentration')
set(gcf,'Position',[100 100 1000 600])
saveas(gcf,'Results/paramest_data.png')
hold off

 
opts = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxIter', 150, 'Diagnostics', 'off', 'Display', 'iter');
initParams=[0.1;0.1;0.1;0.1;0.1; 0.1;0.05;0.05];

% 8 parameters estimated s, k1, k2, k3, k4, k5, km4, km5
loBound = [0.1 0.1 0.1 0.1 0.1 0.1 0.05 0.05]; 
upBound = [0.8 1.5 1.5 1.5 1.5 1.5 0.1 0.1];

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
legend('Estimation', 'True Data')
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
legend('Estimation', 'True Data')
xlabel('Time');
ylabel('yP');
hold off;
      
hold off;
drawnow;
end

function [T,Y] = reactionsolve(a)

s=a(1);
k1=a(2);
k2=a(3);
k3=a(4);
k4=a(5);
k5=a(6);
km4=a(7);
km5=a(8);

% exp(tspan = 0)
x0 = [0.8; 0.6];
t = 0:0.01:12;
[T,Y] = ode45(@reaction, t, x0, []);

function dx = reaction(t,x)

Act=x(1);
yp=x(2);

par.ytot=1;
par.E=0.5;

Act_dot = k1*s+k2*yp-k3*Act;
yp_dot = k4*Act*(par.ytot-yp)/(km4+par.ytot-yp)-k5*par.E*yp/(km5+yp);

dx = zeros(2,1);
dx(1)=Act_dot;
dx(2)=yp_dot;
    
end

end

end


