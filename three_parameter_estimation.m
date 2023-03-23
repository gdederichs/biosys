function results = three_parameter_estimation(initParams) 
% Estimation of the parameters by non-linear least squares optimization

% Course TKP4195, Biosystems Engineering

salut = load("PosFeed_Expdata");
tspan = salut.tspan;
exp = salut.exp;
 
opts = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxIter', 150, 'Diagnostics', 'off', 'Display', 'iter');
initParams=[0.8;0.05;0.05];

% 3 parameters estimated s, km4, km5
loBound = [0.05 0.05 0.05]; 
upBound = [0.8 0.1 0.1];


res = lsqnonlin(@residual,  log10 (initParams), log10 (loBound), log10 (upBound), opts);

   function R = residual(b)
        
       a = 10.^b;
         
       [T ,Y ]= reactionsolve (a);
       
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

function [ T, Y ] = reactionsolve( a )


s=a(1);
km4=a(2);
km5=a(3);

% exp(tspan = 0)
x0 = [0.8; 0.6];

[T,Y] = ode45(@reaction, tspan, x0, []);

function dx = reaction(t, x)

Act=x(1);
yp=x(2);

par.ytot=1;
par.E=0.2;

par.k1 = 0.1;
par.k2 = 1.1324;
par.k3 = 1.5;
par.k4 = 0.1;
par.k5 = 0.502;

Act_dot = par.k1*s+par.k2*yp-par.k3*Act;
yp_dot = par.k4*Act*(par.ytot-yp)/(km4)-par.k5*par.E*yp/(km5+yp);

dx = zeros(2,1);
dx(1)=Act_dot;
dx(2)=yp_dot;
    
end

end

end


