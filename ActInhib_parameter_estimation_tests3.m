function results =ActInhib_parameter_estimation_tests3(initParams) 
% Estimation of the parameters by non-linear least squares optimization
% problem. 
% we will estimate the Activator inhibitor and find parameter s that causes oscillations. 
% Chapter 6.
% we only "know" from experiments we conducted how Act and yP behave, but
% we do not know the parameters. We want to estimate this, such as the
% residuals are minimum (seeedit lecture notes).
%
% run for example: model_parameter_estimation_tests3([0.1;0.1;0.1]) for 3 parameters. 
% run for example: model_parameter_estimation_tests6([0.1;0.1;0.1; 0.5;0.5;0.05 ]) for 6 parameters. 

% Course TKP4195, Biosystems Engineering

sol=[0 0 0]; tspan=0;
load actInhib_Expdata

exp=sol;
   
opts = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxIter', 150, 'Diagnostics', 'off', 'Display', 'iter');
initParams=[0.1;0.1;0.1; 0.5;0.5;0.05]

% positive feedback reaction from 
% we will estimate only 3 parameters, than 5
% first km5, km6, k5, k6, s
loBound = [0.05 0.05 0.1 0.1 0.1]; 
 upBound = [0.1 0.1 1.5 1.5 0.8];
% s Km5 Km6
%loBound = [0 0 0]; 
%upBound = [1 1 1];
% loBound = [0.01]; 
 %upBound = [0.4];

% loBound = [0.0001 .0452 .05 .0444 .0001 250 .007]; 
% upBound = [0.48 1.7 1.7 1.9 .02 350 .008];

res = lsqnonlin(@residual,  log10 (initParams), log10 (loBound), log10 (upBound), opts);

   function R = residual(b)
        
       a = 10.^b;
         
       [T ,Y ]= reactionsolve (a);
       
      % r1= exp(:,1)-Y(:,1);
      % r2= exp(:,2)-Y(:,2);
      % r3= exp(:,3)-Y(:,3);
      % r= r1+r2+r3;
       
       
       residual = exp - Y;
       R = residual(:);
       
       results=a;
       
          subplot(3,1,1);
          plot(T,Y(:,1));
          hold on;
          plot(tspan(:,1),exp(:,1),'o');
          xlabel('Time');
          ylabel('Act');
          hold off;
          str = sprintf('%g | ', a);
          title(str);

         
          subplot(3,1,2);
          plot(T,Y(:,2));
          hold on;
          plot(tspan(:,1),exp(:,2),'o');
          xlabel('Time');
          ylabel('yP');
          hold off;
          
          subplot(3,1,3);
          plot(T,Y(:,3));
          hold on;
          plot(tspan(:,1),exp(:,3),'o');
          xlabel('Time');
          ylabel('yP');
          hold off;
          
          hold off;
          drawnow;
   end

function [ T, Y ] = reactionsolve( a )


s=a(1);
km5=a(2);
km6=a(3);
k5=a(4);
k6=a(5);

x0 = [0.6; 0.8; 0];

[T,Y] = ode45(@reaction, tspan, x0, []);

function dx = reaction(t, x)

R=x(1);
X=x(2);
Ep=x(3);

% Fixed parameters:
par.k0=1;
par.k1=1;
par.k2=1;
par.k22=0;
par.k3=1;
par.k4=0.5;
%par.k5=1;
%par.k6=1;

%par.s=0.5;
%km5=0.08;
%km6=0.09;

par.Etot=1;
par.E=0;

%Diff lign
R_dot=par.k0*Ep+par.k1*s-par.k22*R-par.k2*X*R;
X_dot=par.k3*Ep-par.k4*X;
Ep_dot=(k5*R.*(par.Etot-Ep))./(km5+par.Etot-Ep) - k6*Ep./(km6+Ep);

dx = zeros(3,1);
dx(1)=R_dot;
dx(2)=X_dot;
dx(3)=Ep_dot;
    
end

end

end


