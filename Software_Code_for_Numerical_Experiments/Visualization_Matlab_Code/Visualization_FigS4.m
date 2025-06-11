function Visualization_FigS4()
% Visualization_FigS4.m
% This function generates the Fig. S4 in the paper, which visualizes
% the theoretical limits of a single strain under infinitely many dilutions.
% The theoretical limits are derived from Eq. [S3.1.10] in the SI Appendix.
%
% Author: MingYi Wang, Cornell University
% Last Modified: 06/2025
%
clear; clc;
%% Define the parameters and theoretical limit functions (Eq. [S3.1.10] in the SI Appendix)  
theo_lim = @(r,rho,T) 1./(rho+(1-rho)./(1-exp(-(r.*T+log(rho)))));
Zero_or_not = @(r,rho,T) r.*T+log(rho);
T_comp = @(r, rho) -log(rho)./r;

Tf = 1e2; % Total time horizon for the simulation
tt = 0:0.1:Tf; % Time vector
rr = 0:0.01:1; % Rescaled survival rate vector
[X,Y] = meshgrid(rr,tt);
%% For sensitives only 
r = 1
phase = zeros(length(tt),length(rr));

for i = 1:length(tt)
    for j = 1:length(rr)
        if Zero_or_not(r,rr(j),tt(i)) > 0
            phase(i,j) = theo_lim(r,rr(j),tt(i));
        end
    end
end

figure
contourf(X,Y,phase,20,'EdgeColor','none');
hold on
plot(rr,T_comp(r,rr),'m--','linewidth',2);
hold otheo_lim
axis square
xlabel('$\rho$','FontSize',20,'Interpreter','latex');
ylabel('$T$','FontSize',20,'Interpreter','latex');
colorbar();
ylim([0,5]);
xlim([0,1]);
clim([0 1]);
%% For toxin producers only
epsilon = 0.2;
r = 0.85*(1-epsilon);
phase = zeros(length(tt),length(rr));

for i = 1:length(tt)
    for j = 1:length(rr)
        if Zero_or_not(r,rr(j),tt(i)) > 0
            phase(i,j) = theo_lim(r,rr(j),tt(i));
        end
    end
end

figure
contourf(X,Y,phase,20,'EdgeColor','none');
hold on
plot(rr,T_comp(r,rr),'m--','linewidth',2);
hold otheo_lim
axis square
xlabel('$\rho$','FontSize',20,'Interpreter','latex');
ylabel('$T$','FontSize',20,'Interpreter','latex');
colorbar();
ylim([0,5]);
xlim([0,1]);
clim([0 1]);

end