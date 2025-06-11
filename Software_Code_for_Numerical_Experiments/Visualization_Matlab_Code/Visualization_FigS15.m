function Visualization_FigS15()
% Visualization_FigS15.m
% This function generates the Fig. S15 in the paper, which visualizes
% the Monte Carlo simulations for constitutive killers but with random outcomes after each dilution.
% （Specifically for epsilon == 0.1 and r_ks == 1)
% The visualization includes the mean population size after a series of random dilutions
% for different initial conditions (ICs) and the comparison between two sampling methods:
% "Binomial sampling" and "Binomial sampling + Picking IC uniformly random in a cell."
% 
%
% Author: MingYi Wang, Cornell University
% Last Modified: 06/2025
%
clear; clc;
%% Initialization
n = 800; % number of grid points in each direction
h = 1/n; % grid spacing
dt = h; % time step size

epsilon = 0.1; % The cost of constitutive toxin production
gamma = 1.0; % The rescaled cost of toxin-production rate
rks = 1.0; % The rescaled ratio of the growth rate of the killer to that of the sensitive
lamb = 1.0; % The rescaled arrival rate
rho = 0.65; % The rescaled survival rate

Kc = 1e5; % shared carrying capacity
Num_Samples = 1e5; % Number of samples for the Monte Carlo simulation
num_dilution = 200; % Number of dilutions for the simulation
Tp = 1; % Inter-dilution time

xc_list = linspace(0.05,0.15,5); % Initial fraction of the killer (f_0)
yc_list = linspace(0.1,0.9,9); % Initial total population (N_0)
diff_x = (xc_list(2) - xc_list(1))/2 % Difference between two points in xc_list
diff_y = (yc_list(2) - yc_list(1))/2 % Difference between two points in yc_list
%% Fig. S12A: MC simulations: Binomial sampling
tic
% Initialize lists to store the results
pdf_f_list = cell(length(xc_list),length(yc_list));
pdf_N_list = cell(length(xc_list),length(yc_list));

for i = 1:length(xc_list)
    xc = xc_list(i);
    for j = 1:length(yc_list)
        yc = yc_list(j);
        pdf_f = zeros(1,Num_Samples);
        pdf_N = zeros(1,Num_Samples);
        parfor ii = 1:Num_Samples
            xc_sample = xc;
            yc_sample = yc;
            % Call the function to obtain the population dynamics for 200 dilutions
            [xlist_full,ylist_full,xstart, ystart, xend, yend] ...
                = Const_path_DilutionBinomial(Tp,num_dilution,h,rks,gamma,epsilon,dt,xc_sample,yc_sample,rho,Kc);
            % Store the final fraction of the killer and the final total population
            pdf_f(ii) = xlist_full(end);
            pdf_N(ii) = ylist_full(end);

            ii
        end
        pdf_f_list{i,j} = pdf_f;
        pdf_N_list{i,j} = pdf_N;
    end
    
end
t_total = toc
fprintf('The elapsed time of MC simulations with Binomial sampling is %.3f seconds.\n',t_total)
%% Fig. S12A: plotting the mean table

% The mean table is a 2D matrix where each entry corresponds to the mean fraction
mean_table = zeros(length(yc_list),length(xc_list));
[X,Y] = meshgrid(xc_list,yc_list);
for i = 1:length(xc_list)
    for j = 1:length(yc_list)
        pdf_f =  pdf_f_list{i,j} ; 
        mean_table(j,i) = mean(pdf_f);
    end
end
load("const_bd.mat","const_bd");
figure
imagesc(xc_list,yc_list,mean_table);
hold on
contour(X,Y,const_bd,'k--','linewidth',2);
hold off
xlabel('Initial fraction of the killer (x = f_0)','FontSize',14);
ylabel('Initial total population (y = N_0)','FontSize',14);
xticks(xc_list);
yticks(yc_list);
colorbar();
colormap("jet");
clim([0 1]);
pbaspect([1 2 1])
set(gca,'YDir','normal');
%% Fig. S12B: MC simulations: Binomial sampling + Picking IC uniformly random in a cell
tic
pdf_f_list_unif = cell(length(xc_list),length(yc_list));
pdf_N_list_unif = cell(length(xc_list),length(yc_list));

for i = 1:length(xc_list)
    xc = xc_list(i);
    for j = 1:length(yc_list)
        yc = yc_list(j);
        pdf_f = zeros(1,Num_Samples);
        pdf_N = zeros(1,Num_Samples);
        parfor ii = 1:Num_Samples

            xc_low = xc - diff_x; % Lower bound for xc "cell"
            xc_up = xc + diff_x;  % Upper bound for xc "cell"
            yc_low = yc - diff_y; % Lower bound for yc "cell"
            yc_up = yc + diff_y;  % Upper bound for yc "cell"

            xc_sample = xc_low + (xc_up - xc_low)*rand; % Sample xc uniformly in the cell
            yc_sample = yc_low + (yc_up - yc_low)*rand; % Sample yc uniformly in the cell
            
            % Call the function to obtain the population dynamics for 200 dilutions
            [xlist_full,ylist_full,xstart, ystart, xend, yend] ...
                = Const_path_DilutionBinomial(Tp,num_dilution,h,rks,gamma,epsilon,dt,xc_sample,yc_sample,rho,Kc);

            pdf_f(ii) = xlist_full(end);
            pdf_N(ii) = ylist_full(end);

            ii
        end
        pdf_f_list_unif{i,j} = pdf_f;
        pdf_N_list_unif{i,j} = pdf_N;
    end
    
end
t_total = toc
fprintf('The elapsed time of MC simulations with uniform IC is %.3f seconds.\n',t_total)
%% Fig. S12B: plotting the mean table
mean_table = zeros(length(yc_list),length(xc_list));
[X,Y] = meshgrid(xc_list,yc_list);
for i = 1:length(xc_list)
    for j = 1:length(yc_list)
        pdf_f =  pdf_f_list_unif{i,j} ; 
        mean_table(j,i) = mean(pdf_f);
    end
end

figure
imagesc(xc_list,yc_list,mean_table);
hold on
contour(X,Y,const_bd,'k--','linewidth',2);
hold off
xlabel('Mean initial fraction of the killer (x = f_0)','FontSize',14);
ylabel('Mean initial total population (y = N_0)','FontSize',14);
xticks(xc_list);
yticks(yc_list);
colorbar();
colormap("jet");
clim([0 1]);
pbaspect([1 2 1])
set(gca,'YDir','normal');

%% Fig. S12C: Plotting two sample paths
load("killer_win_FigS15.mat","nk_win","ns_win","f_win");
load("killer_lose_FigS15.mat","nk_lose","ns_lose","f_lose");

mygreen = [50,205,50]/255;
nn = length(nk_win);
t_list = dt.*(1:nn);
% figure 
tiledlayout(2,1)
nexttile
hold on 
plot(t_list,nk_win,'r-','linewidth',1.5);
plot(t_list,ns_win,'b-','linewidth',1.5);
hold off

yyaxis right
plot(t_list,f_win,'-.','linewidth',2.5,'color',mygreen);

legend('Constitutive killer ($n_K$)','Sensitive ($n_S$)','Fraction of the killer ($f$)',...
    'fontsize',12,'location','southoutside','interpreter','latex');
title('Killers win','fontsize',15)
grid on
% grid minor
xlim([0,205]);
ylim([0.05,1.01]);
ax = gca;
ax.YAxis(2).Color = mygreen;


nexttile
hold on 
plot(t_list,nk_lose,'r-','linewidth',1.5);
plot(t_list,ns_lose,'b-','linewidth',1.5);
hold off
xlabel("Rescaled time (t)",'fontsize',18);


yyaxis right
plot(t_list,f_lose,'-.','linewidth',2.5,'color',mygreen);


title('Sensitives win','fontsize',15)
grid on
xlim([0,205]);
ax = gca;
ax.YAxis(2).Color = mygreen;

ttt = text(-20,0.03,'Normalized population','fontsize',18);
ttt2 = text(225,0.06,'Fraction of the killer','fontsize',18,'Color',mygreen);
ttt.Rotation = 90;
ttt2.Rotation = 90;

end