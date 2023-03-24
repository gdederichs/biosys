function results = seven_parameter_estimation(initParams) 
% Estimation of the parameters by non-linear least squares optimization

data = load("PosFeed_Expdata");
tspan = data.tspan;
exp = data.exp;

% Adding noise (var = 0.01, sd = 0.1)
for i = 1:length(tspan)
    tspan(i) = tspan(i)+randn(1,1)*0.04; 
    exp(i) = exp(i)+randn(1,1)*0.04;
end

opts = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxIter', 150, 'Diagnostics', 'off', 'Display', 'iter');
initParams=[0.1;0.1;0.1;0.1;0.1; 0.1;0.05;0.05];

% 8 parameters estimated s, k1, k2, k3, k4, k5, km4, km5
loBound = [0.1 0.1 0.1 0.1 0.1 0.1 0.05 0.05]; 
upBound = [0.8 1.5 1.5 1.5 1.5 1.5 0.1 0.1];

lsqnonlin(@residual,  log10 (initParams), log10 (loBound), log10 (upBound), opts);

function R = residual(b)
        
   a = 10.^b;
         
   [T,Y]= reactionsolve(a);
       
   residual = exp - Y;
   R = residual(:);
       
   results=a;
       
   subplot(2,1,1);
   plot(T,Y(:,1));
   hold on;
   plot(tspan(:,1),exp(:,1),'o');
   xlabel('Time');
   ylabel('Act');
   hold off;
   str = sprintf('%g | ', a);
   title(str);

         
   subplot(2,1,2);
   plot(T,Y(:,2));
   hold on;
   plot(tspan(:,1),exp(:,2),'o');
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
    x0 = [exp(1,1); exp(1,2)];

    [T,Y] = ode45(@reaction, tspan, x0, []);

function dx = reaction(t,x)

    Act=x(1);
    yp=x(2);

    par.ytot=1;
    par.E=0.2;

    Act_dot = k1*s+k2*yp-k3*Act;
    yp_dot = k4*Act*(par.ytot-yp)/(km4+par.ytot-yp)-k5*par.E*yp/(km5+yp);

    dx = zeros(2,1);
    dx(1)=Act_dot;
    dx(2)=yp_dot;
    
end

end

end

