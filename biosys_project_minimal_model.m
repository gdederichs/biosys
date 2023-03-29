clc
clear

%% Minimal Model Simulation
% Initial cdts and params
x0 = init_cond();
l_list = [8 3.2 3 2];
LacY_list = [3 1.3 1.2 1];
lext_list = 0:0.01:7;
par = param();

% Timespan for simulation
dt    = 1e-2 ;  
tlast = 20.000 ;
iterations = fix(tlast/dt) ;
tspan = dt*(0:iterations-1) ; 
options = [];

%% Nullclines
plot_NullCline(6,12);
vectorField(par,6, 12);
title('Nullclines')
xlabel('[LacY]')
ylabel('[Lactose]')
legend('dLacY/dt = 0','dl/dt = 0')
view([90 -90])
set(gcf,'Position',[100 100 1000 600])
saveas(gcf,'Results/nullclines.png')
hold off

%% Temporal evolution with varying init cdt
figure()
hold on
grid minor
for ii = 1:min(length(l_list), length(LacY_list))
    subplot(min(length(l_list), length(LacY_list))/2,2,ii)
    [last_l, last_LacY] = plot_t_evo(tspan, set_x0(l_list(ii), LacY_list(ii)), options, par,1);
    ylim([0 max(max(last_l, last_LacY), max(l_list(ii), LacY_list(ii)))+1])
end

set(gcf,'Position',[10 10 1500 900])
saveas(gcf,'Results/temporal_evo.png')
hold off

%% Bifurcation
final_LacY83 = [];
final_LacY21 = [];

for i = 1:length(lext_list)
    par.lext = lext_list(i);
    
    [final_l_83, final_LacY_83] = plot_t_evo(tspan, set_x0(8, 3), options, par,0);
    [final_l_21, final_LacY_21] = plot_t_evo(tspan, set_x0(2, 1), options, par,0);

    final_LacY83 = [final_LacY83 final_LacY_83];
    final_LacY21 = [final_LacY21 final_LacY_21];
    
end

% Values of lext for which bistability is no longer observed
for iii = 1:length(lext_list)
    if iii~=1 && abs(final_LacY83(iii)-final_LacY21(iii)) < 0.001 && abs(final_LacY83(iii-1)-final_LacY21(iii-1)) > 0.001
        largest_lext = lext_list(iii)
    elseif iii~=length(lext_list) && abs(final_LacY83(iii)-final_LacY21(iii)) < 0.001 && abs(final_LacY83(iii+1)-final_LacY21(iii+1)) > 0.001
        smallest_lext = lext_list(iii)
    end
end

% Plotting of bifurcation diagram
figure()
hold on
grid minor
plot(lext_list, final_LacY21, lext_list, final_LacY83)
legend("l_0 = 2 and LacY_0 = 1", "l_0 = 8 and LacY_0 = 3", Location = "southeast")
title("Bifurcation diagram")
ylabel("[LacY]_s_s")
xlabel("[lext]")
set(gcf,'Position',[100 100 1000 600])
saveas(gcf,'Results/bifurcation.png')

%% Sensitivity analysis
[t,S] = ode45(@S_dot,tspan,zeros([12,1]));
S = reshape(S.',2,6,[]);
S_normalized = S(:, :, end)/max(abs(S(:, :, end)), [], 'all')


%% Functions
function dsdt = S_dot(t, S)
beta = 1;
gamma = 1;
delta = 0.2;
sigma = 1;
l0 = 4;
p= 4;
lext = 1;

syms l LacY

l_dotA = beta*lext*LacY-gamma*l;
LacY_dotA = delta+p*l^4/(l^4+l0^4)-sigma*LacY;

A = matlabFunction(jacobian([l_dotA,LacY_dotA],[l LacY]));

syms betaa gammaa deltaa sigmaa pp l00
ll=1;
LacYY=1;
l_dotB = betaa*lext*LacYY-gammaa*ll;
LacY_dotB = deltaa+pp*ll^4/(ll^4+l00^4)-sigmaa*LacYY;

B = matlabFunction(jacobian([l_dotB,LacY_dotB],[betaa gammaa deltaa sigmaa pp l00]));

ds = A(ll)*reshape(S, [2, 6]) + B(l0, p);
dsdt = reshape(ds, [12, 1]);
end

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

function [last_l, last_LacY] = plot_t_evo(tspan, x0, options, par, want_plot)
[t,x] = ode45(@diff_eq,tspan,x0,options,par);
l = x(:, 1);
LacY = x(:, 2);

if want_plot == 1
    hold on
    grid minor
    plot(t, l, t, LacY)
    legend('Intracellular Lactose', 'LacY')
    xlabel('time')
    ylabel('concentration')
    title("Time evolution with l_0 = "+x0(1)+" and LacY_0 = "+x0(2))
end

last_LacY = LacY(end);
last_l = l(end);
end