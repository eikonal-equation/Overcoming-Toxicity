function Visualization_FigS6()
% Visualization_FigS6.m
% This function generates the Fig. S6 in the paper, which visualizes
% the population dynamics of a single strain under random dilutions.
% The visualization includes the contour plots of the mean population size
% after a series of random dilutions for both sensitive and toxin-producing strains.
%
% Author: MingYi Wang, Cornell University
% Last Modified: 06/2025
%
clear; clc;
%%  Parameter values
% Define the parameters and theoretical limit functions (Eq. [S3.1.10] in the SI Appendix)  
theo_lim = @(r,rho,T) 1./(rho+(1-rho)./(1-exp(-(r.*T+log(rho)))));
Zero_or_not = @(r,rho,T) r.*T+log(rho);
T_comp = @(r, rho) -log(rho)./r;
% Theoretial result after a single dilution (Eq. [S3.1.4] in the SI Appendix)
pop_single = @(t,x,r) 1./(1 + ((1-x)./x)*exp(-r.*t));

Tf = 1e2; % Total time horizon for the simulation
num_dilution = 200; % Number of dilutions for the simulation
xinit = 0.5; % Initial fraction of the killer for all samples
Num_Samples = 1e5; % Number of samples for the Monte Carlo simulation
% For comparison purposes when plotting
tt = 0:0.1:Tf; % Time vector
rr = 0:0.01:1; % Rescaled survival rate vector
%% Sensitives only; r = 1
r = 1
inv_lamb_list = [0.01:0.05:5.01]; % The inverse of the rescaled survival rate
lamb_list = 1./inv_lamb_list; % The rescaled survival rate
rho_list = 0:0.01:1; % Rescaled survival rate vector
mean_table_sensitive = zeros(length(lamb_list),length(rho_list)); % Initialize the mean table for sensitive population
tic
for i = 1:length(lamb_list)
    lamb = lamb_list(i)
    for j = 1:length(rho_list)
        rho = rho_list(j)
        pop_list = zeros(Num_Samples,1); % Initialize the population list for each sample
        parfor ii = 1:Num_Samples
            % Call the single population random dilution function 
            % to get the population after random dilutions
            single_pop_dilutionIndex = singlePop_randDilution(num_dilution,r,xinit,rho,lamb);
            pop_list(ii) = single_pop_dilutionIndex(end);

        end
        mean_table_sensitive(i,j) = mean(pop_list);
    end
end
t_mc = toc;
fprintf('The elapsed time of completing the Monte Carlo simulation (single strain;sensitive) is %.3f seconds.\n',t_mc)

%% Fig. S6A; sensitives only
[X,Y] = meshgrid(rho_list,1./lamb_list);
figure
contourf(X,Y,mean_table_sensitive,10,'EdgeColor','none');
hold on
plot(rr,T_comp(r,rr),'m--','linewidth',2);
hold off
xlabel('$\rho$','FontSize',18,'Interpreter','latex');
ylabel('$1/\lambda$','FontSize',18,'Interpreter','latex');
colorbar();
ylim([0,5]);
xlim([0,1]);
clim([0 1]);
axis square
%% Toxin-producing strain only
r = 0.68; % Rescaled growth rate of the toxin-producing strain
inv_lamb_list = [0.01:0.05:5.01];
lamb_list = 1./inv_lamb_list;
rho_list = 0:0.01:1;
mean_table_killer = zeros(length(lamb_list),length(rho_list));
tic
for i = 1:length(lamb_list)
    lamb = lamb_list(i)
    for j = 1:length(rho_list)
        rho = rho_list(j)
        pop_list = zeros(Num_Samples,1);
        parfor ii = 1:Num_Samples
            single_pop_dilutionIndex = singlePop_randDilution(num_dilution,r,xinit,rho,lamb);
            pop_list(ii) = single_pop_dilutionIndex(end);
        end
        mean_table_killer(i,j) = mean(pop_list);
    end
end
t_mc = toc
fprintf('The elapsed time of completing the Monte Carlo simulation (single strain; toxin-producers) is %.3f seconds.\n',t_mc)
%% Fig. S6B; toxin-producers only
r = 0.68
[X,Y] = meshgrid(rho_list,1./lamb_list);
figure
contourf(X,Y,mean_table_killer,15,'Edgecolor','none');
hold on
plot(rr,T_comp(r,rr),'m--','linewidth',2);
hold off
xlabel('$\rho$','FontSize',18,'Interpreter','latex');
ylabel('$1/\lambda$','FontSize',18,'Interpreter','latex');
colorbar();
ylim([0,5]);
xlim([0,1]);
clim([0 1]);
axis square

end