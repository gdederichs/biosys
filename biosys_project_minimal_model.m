clc
clear

%% Minimal Model Simulation
% Initial cdts and params
x0 = init_cond();
par = param();
l_list = [8 3.2 3 2];
LacY_list = [3 1.3 1.2 1];
% Timespan for simulation
tspan = [0 20]; 
options = [];

for i = 1:length(l_list)
    plot_t_evo(tspan, set_x0(l_list(i), LacY_list(i)), options, par)
end

%% Jacobian, augmentation
% syms l LacY beta lext gamma delta sigma p l0
% 
% l_dot = beta*lext*LacY-gamma*l;
% LacY_dot = delta+p*l^4/(l^4+l0^4)-sigma*LacY;
% 
% A = jacobian([l_dot,LacY_dot],[l LacY]);
% B = jacobian([l_dot,LacY_dot],[beta lext gamma delta sigma p l0]);
% 
% % S_dot = A*S + B

%% Nullclines
plot_NullCline(6,12)
vectorField(par,6, 12)
title('Nullclines')
xlabel('[LacY]')
ylabel('[Lactose]')
legend('dLacY/dt = 0','dl/dt = 0')
view([90 -90])



%% Functions
function nl = plot_NullCline(x_lim, y_lim)
par = param();

syms l LacY 
LacY_dot = par.delta+par.p*l^4/(l^4+par.l0^4)-par.sigma*LacY;
l_dot = par.beta*par.lext*LacY-par.gamma*l;


figure()
hold on
grid minor
nl1=ezplot(LacY_dot,[0,x_lim,0,y_lim]);
set(nl1,color = 'r');
nl2=ezplot(l_dot,[0,x_lim,0,y_lim]);
set(nl2,color = 'b');
nl=[nl1,nl2];
end

function vectorField(par, x_lim, y_lim)
[LacY ,l] = meshgrid(0:0.5:x_lim, 0:0.5:y_lim);

LacY_dot = par.delta+par.p.*l.^4./(l.^4+par.l0.^4)-par.sigma.*LacY;
l_dot = par.beta.*par.lext.*LacY-par.gamma.*l;

quiver(LacY, l, LacY_dot, l_dot)
xlim([0 x_lim])
ylim([0 y_lim])

end



function dxdt = diff_eq(t,x,par)
% Variables 
l = x(1);
LacY = x(2);

% Differential equations
l_dot = par.beta*par.lext*LacY-par.gamma*l;
LacY_dot = par.delta+par.p*l^4/(l^4+par.l0^4)-par.sigma*LacY;


dxdt = [l_dot;LacY_dot];
end


function par = param()
par.beta = 1;
par.gamma = 1;
par.delta = 0.2;
par.sigma = 1;
par.l0 = 4;
par.p = 4;
par.lext = 2.5;
end


function x0 = init_cond()
l_init = 4;
LacY0 = 1;

x0 = [l_init; LacY0];
end

function x0 = set_x0(l_init, LacY0)

x0 = [l_init, LacY0];
end

function plot_t_evo(tspan, x0, options, par)
[t,x] = ode45(@diff_eq,tspan,x0,options,par);
l = x(:, 1);
LacY = x(:, 2);

figure()
hold on
grid minor
plot(t, l, t, LacY)
legend('Intracellular Lactose', 'LacY')
xlabel('time')
ylabel('concentration')
title("Time evolution with l_0 = "+x0(1)+" and LacY_0 = "+x0(2))
end