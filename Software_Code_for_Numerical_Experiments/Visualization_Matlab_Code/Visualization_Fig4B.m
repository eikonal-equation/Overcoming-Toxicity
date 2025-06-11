function Visualization_Fig4B()
% Visualization_Fig4B.m
% This function generates the Fig. 4B in the paper, which visualizes
% the population dynamics of a system with 
% constitutive toxin production (Eq.[2] in the paper with "delta == 0" and "a == 1".)
% The visualization includes the population dynamics of the system under
% different initial conditions, and the fraction of the killer population over time.
% The initial conditions are chosen such that the killers win in one case
% and the sensitives win in the other case.
%
% Author: MingYi Wang, Cornell University
% Last Modified: 06/2025
%
clear; clc;
%% Parameters
n = 800; % number of mesh points along on side
h = 1/n; % mesh size

epsilon = 0.2; %The cost of constitutive toxin production
gamma = 1.0; %The rescaled cost of toxin-production rate
rks = 0.85; %The rescaled ratio of the growth rate of the killer to that of the sensitive
rho = 0.65; % rescaled survival rate
lamb = 1; % rescaled arrival rate of the population

% Policy matrix of constitutive toxin production
Policy_mat = ones(n+1);
Policy_mat(1,1) = 0; % For visualization purpose, we set the first element to 0

num_dilution = 100; % number of dilutions
dt = h; % time step for the ODE solver

%% Killers win
% Initial condition for the killers win
xc = 0.5;
yc = 0.7;
% Generate paths for 100 dilutions
choice = 'conservative'
[xlist_full,ylist_full,policy_list,xstart, ystart, xend, yend] ...
    = Optimal_path_DilutionIndexed(num_dilution,h,Policy_mat,lamb,rks,gamma,epsilon,dt,xc,yc,choice,rho);
% Get the actual normalized populations
nk_win = xlist_full.*ylist_full;
ns_win = (1-xlist_full).*ylist_full;
f_win = xlist_full;
N_win = ylist_full;

save('killer win.mat','nk_win','ns_win','f_win','N_win');

%% Sensitives win
% Initial condition for the sensitives win
xc = 0.35;
yc = 0.7;
% Generate paths for 100 dilutions
choice = 'conservative'
[xlist_full,ylist_full,policy_list,xstart, ystart, xend, yend] ...
    = Optimal_path_DilutionIndexed(num_dilution,h,Policy_mat,lamb,rks,gamma,epsilon,dt,xc,yc,choice,rho);
% Get the actual normalized populations
nk_lose = xlist_full.*ylist_full;
ns_lose = (1-xlist_full).*ylist_full;

f_lose = xlist_full;
N_lose = ylist_full;

save('killer lose.mat','nk_lose','ns_lose','f_lose','N_lose');

%% Plotting them together
mygreen = [50,205,50]/255;
nn = length(nk_win);
t_list = dt.*(1:nn);

tiledlayout(2,1)
nexttile
hold on
plot(t_list,nk_win,'r-','linewidth',1.5);
plot(t_list,ns_win,'b-','linewidth',1.5);
hold off
yyaxis right
plot(t_list,f_win,'-.','linewidth',2.5,'color',mygreen);

legend('Constitutive killer ($n_K$)','Sensitive ($n_S$)','Fraction of the killer ($f$)',...
    'fontsize',12,'location','best','interpreter','latex');
title('Killers win','fontsize',15)
grid on
xlim([0,82]);
ylim([0.5,1.01]);
ax = gca;
ax.YAxis(2).Color = mygreen;

nexttile
hold on
plot(t_list,nk_lose,'r-','linewidth',1.5);
plot(t_list,ns_lose,'b-','linewidth',1.5);
hold off
xlabel("Rescaled time (t)",'fontsize',18);
% ylabel("Normalized population",'fontsize',18);

yyaxis right
plot(t_list,f_lose,'-.','linewidth',2.5,'color',mygreen);

% legend('Constitutive killer','Sensitive','Freq of the killer (f)','fontsize',15,'location','best');
title('Sensitives win','fontsize',15)
grid on
% grid minor
xlim([0,82]);
ylim([0,0.44]);
ax = gca;
ax.YAxis(2).Color = mygreen;

ttt = text(-9,0.1,'Normalized population','fontsize',18);
ttt2 = text(89,0.17,'Fraction of the killer','fontsize',18,'Color',mygreen);
ttt.Rotation = 90;
ttt2.Rotation = 90;

end
