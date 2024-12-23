function []=PosFeed_Exp_generator()
% main function. 
% you run it by simply runing positiveFeedback in the command line. Sjekk at du ligger i riktig mappe (samme som alle filene. opprett ny mappe!)
% jeg delte her til funksjoner: 
% ode15 løser diff ligning fra tiden t0 til tend, gitt initielle betingelser x0 (x01 og x02);
% Tine, prøv å debug. for dette kan du markere med musen en linje (på venstre siden) slik at en rød markering åpnes, og når du kjøre den, den skal stoppe der og vente på ordre fra deg. da kan du bruke F5 å fortsette, F10 å fortsette en linje, osv. slik ser du hva som skjer line for linje. 
% PosFeed_ODE er diff ligning fil. bruk den som mal! 
close all
color='b';
par=PosFeedParameters_test();
% par.s=0.05; S line above koshland-golbeter (one cutting point)
% par.s=0.4; S line below KG function (one cutting point)
options=[];
[t,x]=ode15s(@PosFeed_ODE,[0:4:15],[0.8 0.6],options,par);

PlotTime(t,x);
tspan=t; exp=x;
save PosFeed_Expdata  tspan exp
end

%% Plot:

function []=PlotTraj(t,x,color)
% just plot the trajoctries. Straight forward!
% -----------------
% By Nadi S Bar, 
% For the course TKP4195
% January 2014

if nargin<3
    color='b'; % in case color is not sent to the function, default. 
end

figure(1)
plot(x(:,1),x(:,2),'-','color',color);
hold on
plot(x(1,1),x(1,2),'*r');
plot(x(end,1),x(end,2),'ob','MarkerFaceColor','k');
xlabel('Activator');ylabel('Y_p');

end
%% phase plot

function []=PhasePlot(par,fig)
% plots the phase diagram using quiver
% do not forget the .* and ./ operator in the ODE Act_dot and yP_dot, 
% since you are feeding matrices, and
% you actualy need to calculate point wise.
% -----------------
% By Nadi S Bar, 
% For the course TKP4195
% January 2014

figure(fig);
hold on

[Act, yP] = meshgrid(0:0.05:1, 0:.05:1);

Act_dot=par.k1*par.s+par.k2*yP-par.k3*Act;
yP_dot=(par.k4*Act.*(par.yT-yP))./(par.km4+par.yT-yP) - par.k5*par.E*yP./(par.km5+yP);

quiver(Act,yP,Act_dot, yP_dot)


end
%% plot time 
function []=PlotTime(t,x)
% -----------------
% By Nadi S Bar, 
% For the course TKP4195
% January 2014


figure(2)
subplot(2,1,1);
plot(t,x(:,1),'b-');
ylabel('Activator')
xlabel('Time (s)')
subplot(2,1,2);
plot(t,x(:,2),'r--','linewidth',2);
xlabel('Time (s)')
ylabel('Y_p');

end
%% Nullclines

function PlotNullCline(par,fig,color)
% Plot the nullclines
% Uses symbolic toolbox
% write the ODE funciton f1...n as strings first
% then evaluate them useing eval
% then plot them using ezplot (function plot)
% -----------------
% By Nadi S Bar, 
% For the course TKP4195
% January 2014
figure(fig)
hold on;

str1='par.k1*par.s+par.k2*yP-par.k3*Act'
str2='(par.k4*Act.*(par.yT-yP))./(par.km4+par.yT-yP) - par.k5*par.E*yP./(par.km5+yP)'
syms Act yP

Null1=eval(str1);
Null2=eval(str2);
h=ezplot(Null1,[0,1,0,1]);set(h,'color',color,'linewidth',2);grid;hold on
h=ezplot(Null2,[0,1,0,1]);set(h,'color','b','linewidth',2);grid;hold on
end

function par=PosFeedParameters_test()

par.k1=1;
par.k2=0.8;
par.k3=1.2;
par.k4=1;
par.k5=1;
par.s=0.1;
par.km4=0.05;
par.km5=0.05;
par.yT=1;
par.E=0.5;
par.nypar=1;
end

function x_dot=PosFeed_ODE(t,x,par)
% Initielle betingelser 
Act=x(1);
yP=x(2);
% Diff lign
Act_dot=par.k1*par.s+par.k2*yP-par.k3*Act;
yP_dot=(par.k4*Act.*(par.yT-yP))./(par.km4+par.yT-yP) - par.k5*par.E*yP./(par.km5+yP);
x_dot=[Act_dot;yP_dot];
end
