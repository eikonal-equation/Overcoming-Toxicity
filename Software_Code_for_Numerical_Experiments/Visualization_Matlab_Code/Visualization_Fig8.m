function Visualization_Fig8(xc, yc,N,lamb_list,rho_list)
% Visualization_Fig8.m
% This function generates the Fig. 8 in the paper, which visualizes
% the probability performance for the strategic, tactical, and constitutive cases
% for different values of the dilution frequency (lambda) and survival fraction (rho).
% The visualization includes the heat maps of the probability performance
% for the strategic, tactical, and constitutive cases at a specific point (xc, yc).
%
% Input:
%   xc: x-coordinate of the point of interest
%   yc: y-coordinate of the point of interest
%   N: number of points in each direction for the contour plots
%   lamb_list: list of arrival rates (dilution frequencies)
%   rho_list: list of survival rates (survival fractions)
%
% Author: MingYi Wang, Cornell University
% Last Modified: 06/2025
%
clear; clc;
%% Create tables to store the value function data
strategic_table = zeros(length(lamb_list),length(rho_list));
tactical_table = zeros(length(lamb_list),length(rho_list));
const_table = zeros(length(lamb_list),length(rho_list));
xx = linspace(0,1,N); % sample points in x-direction
xindx = find(xx == xc,1); % Find the index of xc in xx
yindx = find(xx == yc,1); % Find the index of yc in xx
data_precision = 'double';
% Loop over the lambda and rho values
for ii = 1:length(lamb_list)

    lamb = lamb_list(ii)

    for jj = 1:length(rho_list)

        rho = rho_list(jj)

        % Value function file for the strategically-optimal case
        wFile = fopen(["Sto_strategic_vertical_valuefn_rho",num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat']);
        perf_strategic = fread(wFile, data_precision);
        perf_strategic=  reshape(perf_strategic,[N,N]);
        fclose(wFile);

    
        % Value function file for the constitutive case
        wFile = fopen(["Perf_vertical_const_valuefn_rho",num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat']);
        perf_const = fread(wFile, data_precision);
        perf_const=  reshape(perf_const,[N,N]);
        fclose(wFile);

    
        % Value function file for the tactically-optimal case
        wFile = fopen(["Perf_vertical_tactic_valuefn_rho",num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat']);
        perf_tactic = fread(wFile, data_precision);
        perf_tactic=  reshape(perf_tactic,[N,N]);
        fclose(wFile);

        %  Store the performance metrics in the respective tables
        strategic_table(ii, jj) = perf_strategic(yindx, xindx);
        tactical_table(ii, jj) = perf_tactic(yindx, xindx);
        const_table(ii, jj) = perf_const(yindx, xindx);

    end
end

%% Constitutive
figure
imagesc(rho_list,lamb_list,const_table);
view([0 0 1]);
xlabel('Survival fraction ($\rho$)','FontSize',20,'Interpreter','latex');
ylabel('Dilution frequency ($\lambda$)','FontSize',20,'Interpreter','latex');
colorbar();
yticks(lamb_list);
clim([0 1]);
colormap("hot")
axis square
set(gca,'YDir','normal');

%% Tactically-optimal
figure
imagesc(rho_list,lamb_list,tactical_table);
view([0 0 1]);
xlabel('Survival fraction ($\rho$)','FontSize',20,'Interpreter','latex');
ylabel('Dilution frequency ($\lambda$)','FontSize',20,'Interpreter','latex');
colorbar();
yticks(lamb_list);
clim([0 1]);
colormap("hot")
axis square
set(gca,'YDir','normal');

%% Strategic vs. Tactical
figure
imagesc(rho_list,lamb_list,strategic_table - tactical_table);
view([0 0 1]);
xlabel('Survival fraction ($\rho$)','FontSize',20,'Interpreter','latex');
ylabel('Dilution frequency ($\lambda$)','FontSize',20,'Interpreter','latex');
colorbar();
yticks(lamb_list);
colormap("cool")
axis square
set(gca,'YDir','normal');

end
