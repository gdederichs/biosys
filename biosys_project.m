clc
clear

%% Lactose input at lext = 0; arbitrary initial conditions
% Initial cdts and params
x0 = init_cond();
par = param();
% Timespan for simulation
tspan = [0 50]; 
options = [];

[t1,x] = ode45(@diff_eq,tspan,x0,options,par);
l = x(:, 1);
LacI = x(:, 2);
LacY = x(:, 3);
LacAZ = x(:, 4);

% Simulation of the Lac operon without lactose input

% figure(1)
% plot(t1, l, t1, LacI, t1, LacY, t1, LacAZ)
% legend('Intracellular Lactose', 'LacI', 'LacY', 'LacA and LacZ')
% xlabel('time')
% ylabel('concentration')

%% Lactose input at lext = 1; initial cdts are steady state of previous
x0 = [l(end), LacI(end), LacY(end), LacAZ(end)]; %x0 = [lin; LacI0; LacY0; LacAZ0];
par.lext = 1;

[t2,x] = ode45(@diff_eq,tspan,x0,options,par);
l_1 = x(:, 1);
LacI_1 = x(:, 2);
LacY_1 = x(:, 3);
LacAZ_1 = x(:, 4);

% Simulation of the lac operon with lactose input

% figure(2)
% plot(t2, l_1, t2, LacI_1, t2, LacY_1, t2, LacAZ_1)
% legend('Intracellular Lactose', 'LacI', 'LacY', 'LacA and LacZ')
% xlabel('time')
% ylabel('concentration')

%% Merge graphs

l_merge = cat(1,l,l_1);
LacI_merge = cat(1,LacI,LacI_1);
LacY_merge = cat(1,LacY,LacY_1);
LacAZ_merge = cat(1,LacAZ,LacAZ_1);
t_merge = cat(1,t1,t1(end)+t2);

% Plotting entire simulation
figure(3);
hold on
grid minor
plot(t_merge, l_merge, t_merge, LacI_merge, t_merge, LacY_merge, t_merge, LacAZ_merge)
xline(t1(end),Label = "Extracellular Lactose added",LineStyle="--")
legend('Intracellular Lactose', 'LacI', 'LacY', 'LacA and LacZ','',Location = 'northwest')
xlabel('time')
ylabel('concentration')
ylim([0,inf])
set(gcf,'Position',[100 100 1000 600])
saveas(gcf,'Results/complete_evo.png')
hold off

%% Functions
function dxdt = diff_eq(t,x,par)
% Variables 
l = x(1);
LacI = x(2);
LacY = x(3);
LacAZ = x(4);

% Differential equations
l_dot = par.k1*par.lext*LacY - l*par.kd*LacAZ;
LacI_dot = par.k2/(1+(l/par.Kl)^par.n)-par.kdI*LacI;
LacY_dot = par.k3/(1+(LacI/par.KI)^par.m)-par.kdY*LacY;
LacAZ_dot = 2*par.k3/(1+(LacI/par.KI)^par.m)-par.kdAZ*LacAZ;

dxdt = [l_dot;LacI_dot;LacY_dot;LacAZ_dot];
end


function par = param()
par.k1 = 0.9;       % Absorption rate of lext to l
par.k2 = 1;         % Formation of LacI
par.k3 = 3;         % Formation of LacY, LacA and LacZ
par.kd = 0.5;       % Degradation of l
par.kdI = 0.4;      % Degradation of LacI
par.kdY = 0.2;      % Degradation of LacY
par.kdAZ = 0.2;     % Degradation of LacA and LacZ
par.Kl = 0.03;      % Inhibition of LacI (smaller value -> bigger inhinition)
par.KI = 0.05;      % Inhibition of LacY, LacA and LacZ (smaller value -> bigger inhinition)
par.lext = 0.0;     % lext constant
par.n = 1;          % LacI_dot Hill coefficient
par.m = 1;          % LacY_dot and LacAZ_dot Hill coefficient
end

function x0 = init_cond()
lin = 6;
LacI0 = 0.01;
LacY0 = 1;
LacAZ0 = 2*LacY0;

x0 = [lin; LacI0; LacY0; LacAZ0];
end