function Visualization_Fig7(Filename_list, File_parameters)
% Visualization_Fig7.m
% This function generates the Fig. 7 in the paper, which visualizes
% the probability performance metric for the constitutive toxin production case
% and the tactically-optimal toxin production case.
% The visualization includes the contour plots of the probability performance metric
% and the tactically-optimal policy (with superimposed boundary of the strategically-optimal Toxin-On region).
% Input:
%   Filename_list: n array containing the paths to the data files with the following order:
%                  1: performance metric of the constitutive case;
%                  2: performance metric of the tactically-optimal case;
%                  3: optimal policy for the tactically-optimal case;
%                  4: optimal policy for the strategically-optimal case
%   File_parameters: A string representing the path to the parameter file
%
% Author: MingYi Wang, Cornell University
% Last Modified: 06/2025
%
clear; clc;
%% Read in the data files
data_precision1 = 'double';
data_precision2 = 'uint8';
% Read in parameter data
uFile_params = fopen(File_parameters);
u_params = fread(uFile_params,11, data_precision1);
Nx = u_params(1) + 1 % Number of grid points in x-direction
Ny = u_params(2) + 1 % Number of grid points in y-direction
dx = u_params(3) % Grid spacing in x-direction
dy = u_params(4) % Grid spacing in y-direction
epsilon = u_params(5) % Cost of constitutive toxin production
rks = u_params(6) % Rescaled ratio of the growth rate of the killer to that of the sensitive
gamma = u_params(7) % Rescaled cost of toxin-production rate
lamb = 1 % Rescaled arrival rate
rho = 0.65 % Rescaled survival rate
vic_thres = u_params(10) % Victory threshold
defeat_thres = u_params(11) % Defeat threshold
fclose(uFile_params);

% Drift part of the population growth dynamics (Eq.[2] in the paper with "delta == 0")
fx = @(x,y,a)  x.*(1-x).*((1-y).*(rks.*(1-epsilon.*a) - 1) + a.*gamma.*x.*y);
fy = @(x,y,a)  y.*(1-y).*(1 + (rks.*(1-epsilon.*a) - 1).*x) - a.*gamma.*y.^2.*x.*(1-x);


% For constitutive case
wFile = fopen(Filename_list(1));
perf_const = fread(wFile, data_precision1);
perf_const =  reshape(perf_const,[Ny,Nx]);
fclose(wFile);
% For tactically-optimal case
wFile = fopen(Filename_list(2));
perf_tactic = fread(wFile, data_precision1);
perf_tactic = reshape(perf_tactic,[Ny,Nx]);
fclose(wFile);

% Optimal policy file for the tactically-optimal case
aFile = fopen(Filename_list(3));
policy_tactic = fread(aFile, data_precision2);
policy_tactic =  reshape(policy_tactic,[1601,1601]);
fclose(aFile);

% Optimal policy file for the strategically-optimal case
aFile = fopen(Filename_list(4));
policy_strategic = fread(aFile, data_precision2);
policy_strategic =  reshape(policy_tstrategic,[1601,1601]);
fclose(aFile);
%% Plotting initalization
xx = linspace(0,1,Nx); % sample points in x-direction
yy = linspace(0,1,Ny); % sample points in y-direction
[X,Y] = meshgrid(xx,yy); % create a meshgrid for plotting
load("const_bd.mat","const_bd");
load("myopic_bd.mat","mypoic_bd");
%% Fig. 7A: Plotting the probability performance metric for the constitutive case
figure;
contourf(X,Y,perf_const,12,'Edgecolor','none');
hold on
contour(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),const_bd(2:end-1,2:end-1), [1], 'k--','LineWidth',2);
xline(defeat_thres,'m:','linewidth',2);
xline(vic_thres,'m:','linewidth',2);
hold off
xlabel('Initial fraction of the killer (x = f_0)','FontSize',14);
ylabel('Initial total population (y = N_0)','FontSize',14);
axis equal
hcb = colorbar();
clim([0 1]);
cc = linspace(0,1,10);
set(hcb,'Ticks',cc);
myCticks = arrayfun(@(x) sprintf('%.1f',x),cc,'un',0);
set(hcb,'Ticklabels',myCticks);

%% Fig. 7B: Plotting the probability performance metric for the tactically-optimal case
figure;
contourf(X,Y,perf_tactic,12,'Edgecolor','none');
hold on
contour(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),myopic_bd_bd(2:end-1,2:end-1), [1], 'k--','LineWidth',2);
xline(defeat_thres,'m:','linewidth',2);
xline(vic_thres,'m:','linewidth',2);
hold off
xlabel('Initial fraction of the killer (x = f_0)','FontSize',14);
ylabel('Initial total population (y = N_0)','FontSize',14);
axis equal
hcb = colorbar();
clim([0 1]);
cc = linspace(0,1,10);
set(hcb,'Ticks',cc);
myCticks = arrayfun(@(x) sprintf('%.1f',x),cc,'un',0);
set(hcb,'Ticklabels',myCticks);

%% Fig. 7C: Plotting the tactically-optimal policy
% Compute the vector field
N=20; % Number of grid points we intended to plot the vector field
NN = 1600; % Number of grid points we used to compute the vector field
eps = 1e-6; % Numerical tolerance for zero
x_plot=linspace(0,1,N+1);
y_plot=linspace(0,1,N+1);
[x,y]=meshgrid(x_plot',y_plot');
policy_st_downsample = policy_strategic(1:NN/N:end,1:NN/N:end);
policy_tac_downsample = policy_tactic(1:NN/N:end,1:NN/N:end);
% Compute the rate functions for different "a" values
% a = 0; not producing toxin
df_a0=fx(x,y,0); dN_a0=fy(x,y,0);
dN_a0(abs(dN_a0-0) < eps) = 0;
df_a0(abs(df_a0-0) < eps) = 0;
rr_a0=sqrt(df_a0.^2+dN_a0.^2);
u_a0=df_a0./rr_a0; v_a0=dN_a0./rr_a0;
% a = 1; producing toxin at the maximum rate
df_a1=fx(x,y,1); dN_a1=fy(x,y,1); 
dN_a1(abs(dN_a1-0) < eps) = 0;
df_a1(abs(df_a1-0) < eps) = 0;
rr_a1=sqrt(df_a1.^2+dN_a1.^2);
u_a1=df_a1./rr_a1; v_a1=dN_a1./rr_a1;

figure
contourf(X,Y,policy_tactic, [0, 1],'EdgeColor','none');
hold on
for ii = 2:N
    for ij = [2:1:N+1]
        if isnan(u_a0(ii,ij)) || isnan(v_a0(ii,ij)) || isnan(u_a1(ii,ij)) || isnan(v_a1(ii,ij))
            continue
        else
            ah = annotation('arrow',...
                'Color', [169/255,169/255,169/255],...
                'headStyle','cback1','HeadLength',3.5,'HeadWidth',3.5);
            set(ah,'parent',gca);
            if policy_tac_downsample(ii,ij) == 0
                set(ah,'position',[x(ii,ij) y(ii,ij) 0.03*u_a0(ii,ij) 0.03*v_a0(ii,ij)]);
            else
                set(ah,'position',[x(ii,ij) y(ii,ij) 0.03*u_a1(ii,ij) 0.03*v_a1(ii,ij)]);
            end
        end
    end
end
contour(X(:,1:1:1584),Y(:,1:1:1584),policy_strategic(:,1:1:1584),[1],'r--','LineWidth',2.5);
hold off
xlabel('Initial fraction of the killer (x = f_0)','FontSize',14);
ylabel('Initial total population (y = N_0)','FontSize',14);
axis equal
colormap("copper")

end
