clc
clear

repression_ass7_g

function repression_ass7_g()
% Nadi 
% this is a skeleton program that simulates repression. You can change the
% differential equations with anything you want, to visualize other 2
% variables system. 
% 2014
close all
%colors = repmat('krgbmc',1,300) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 1:  Define constants 

beta = 1 ;
gamma = 1 ;
delta = 0.2 ;
l_0 = 4 ;

p = 4 ;
sigma = 1 ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2:  Define parameters that will be varied 
   
lics = [8,3.2,3,2] ;%%initial conditions
LacYics=[3,1.3,1.2,1];
tests = length(lics) ;
l_ext = [2.5] ; %%diff parameters
trials = length(l_ext) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 3:  Define time step, simulation time, initialize matrices 
dt    = 1e-2 ; % s 
tlast = 20.000 ;  % s
iterations = fix(tlast/dt) ;
time = dt*(0:iterations-1) ;

l_all = zeros(tests,trials,iterations) ;
l_last = zeros(tests,trials) ;
LacY_all = zeros(tests,trials,iterations) ;
LacY_last = zeros(tests,trials) ;

figure(1)
hold on 
figure(2)
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 4:  Run for different initial conditions and values of l_ext 
for iii=1:tests % different initial conditions
    
  for ii=1:trials % different values of l_ext
     
    l = lics(iii) ;
    LacY = LacYics(iii) ;
    x=zeros(2,iterations);
    dx=zeros(2,iterations);
    x(:,1) = [l,LacY];
  
    for i = 1:iterations

      l = x(1,i);
      LacY = x(2,i);
      %R_all(iii,ii,i) = R ;
      
      % Insert the diff equations here:
      dl_Euler = beta*l_ext*LacY-gamma*l;
      dlacY_Euler = delta+p*l^4/(l^4+l_0^4)-sigma*LacY;
      
      % integrate (Euler)
     
      [t,x] = ode45(dl_Euler);
      % save to memory: 
      
    end % of this time step
     

    % plot here (and save to memory) different l_ext values
    

 end

    % plot here (and save to memory) different initial conditions...
    
    
end
end

