function Visualization_Fig4CtoF()
% Visualization_Fig4CtoF.m
% This function generates the Figs. 4C to 4F in the paper, which visualizes
% the time untill competitive exclusion with constitutive toxin production.
% The visualization includes the contour plots of the time till competitive exclusion, 
% (red and blue indicate f=1 and f=0 limits,respectively).
% Figs. 4C and Fig. 4D assume the same survival rate (rho = 0.65) but different inter-dilution times (T = 1 and T = 2, respectively).
% Figs. 4E and Fig. 4F  assume the same basal death rate (delta = 0.43) as Fig. 4A
% while different inter-dilution times (T = 4 and T = 8, respectively).
%
% Author: MingYi Wang, Cornell University
% Last Modified: 06/2025
%
clear; clc;
%% Read in data from files
data_precision1 = 'double';
data_precision2 = 'int32';
% Read in parameter data
filename_parameters = 'Regular_dilutions_constitutive_DomainParameters.dat';
% Paramter Values
uFile_params = fopen(filename_parameters);
u_params = fread(uFile_params,11, data_precision1);

Nx = u_params(1) + 1 % Number of grid points in x-direction
Ny = u_params(2) + 1 % Number of grid points in y-direction
dx = u_params(3) % Grid spacing in x-direction
dy = u_params(4) % Grid spacing in y-direction 
epsilon = u_params(5) % Cost of constitutive toxin production
rks = u_params(6) % Rescaled ratio of the growth rate of the killer to that of the sensitive
gamma = u_params(7) % Rescaled cost of toxin-production rate
Tf = u_params(8) % Regular inter-dilution time
fM = u_params(9) + 1 % Number of time slices
dt = u_params(10) % Time step size
num_dilution = u_params(11) %  Number of dilutions

fclose(uFile_params);

% Read in iteration map data
filename_4C = 'limit_iter_rho650_T1.dat';
filename_4D = 'limit_iter_rho650_T2.dat';
filename_4E = 'limit_iter_rho178500_T4.dat';
filename_4F = 'limit_iter_rho3186_T8.dat';

wFile = fopen(filename_4C);
IterMap_4C = fread(wFile, data_precision2);
IterMap_4C =  reshape(IterMap_4C,[Ny,Nx]);
fclose(wFile);

wFile = fopen(filename_4D);
IterMap_4D = fread(wFile, data_precision2); 
IterMap_4D =  reshape(IterMap_4D,[Ny,Nx]);
fclose(wFile);

wFile = fopen(filename_4E);
IterMap_4E = fread(wFile, data_precision2);
IterMap_4E =  reshape(IterMap_4E,[Ny,Nx]);
fclose(wFile);

wFile = fopen(filename_4F);
IterMap_4F = fread(wFile, data_precision2);
IterMap_4F =  reshape(IterMap_4F,[Ny,Nx]);
fclose(wFile);

%% Plotting
yy = linspace(0,1,Ny);
xx = linspace(0,1,Nx);

[X,Y] = meshgrid(xx,yy);
%% 4C
% Load the necessary data for Fig. 4C from the files generated in Fig. 4A and Fig. 3B
load("Fslice_4A.mat","uu");
load('separatrix_4.mat','y_unstable_backward',"y_stable_backward",'y_unstable_forward',"y_stable_forward");
figure;
contourf(X,Y,IterMap_4C,60,'Edgecolor','none');
hold on
contour(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),uu(2:end-1,2:end-1), [1], 'k--','LineWidth',2);
plot(0.35,0.7,'diamond','markersize',5.5,'markerfacecolor','k','markeredgecolor','k');
plot(0.5,0.7,'^','markersize',5.5,'markerfacecolor','k','markeredgecolor','k');
% Covert the stable and unstable separatrices in the KS-coordinates
%  to fractions of the killer and the total population
nn1 = y_stable_forward(:,1) + y_stable_forward(:,2);
ff1 = y_stable_forward(:,1)./nn1;
ff1(ff1>1) = 1;
nn2 = y_stable_backward(:,1) + y_stable_backward(:,2);
ff2 = y_stable_backward(:,1)./nn2;
ff2(ff2>1) = 1;
indx = find(nn2 < 0,1) - 1;
ff2(ff2<1e-3) = nan;
nn2(nn2<1e-3) = nan;
% Plot the separatrix
plot(ff1,nn1,'-.','LineWidth',2.7,'Color','m');
plot(ff2(1:indx),nn2(1:indx),'-.','LineWidth',2.7,'Color','m');

hold off

xlabel('Initial fraction of the killer (x = f_0)','FontSize',14);
ylabel('Initial total population (y = N_0)','FontSize',14);
axis equal
h = colorbar();
clim([-430,270]);
colormap(centered);
% Set custom ticks for the colorbar
set(h, 'Ticks', [-400, -300, -200, -100, 0, 100, 200]);
% Set corresponding custom labels, converting negative values to positive for the blue region
set(h, 'TickLabels', {'400', '300', '200', '100', '0', '100', '200'});

const_bd = IterMap_4C;
const_bd(IterMap_4C < 0) = 0; % Set the cost boundary
const_bd(IterMap_4C >= 0) = 1; % Set the cost boundary

save("const_bd.mat","const_bd");
%% 4D
figure
contourf(X,Y,IterMap_4D*2,60,'Edgecolor','none');
hold on
plot(0.35,0.7,'diamond','markersize',5.5,'markerfacecolor','k','markeredgecolor','k');
plot(0.5,0.7,'^','markersize',5.5,'markerfacecolor','k','markeredgecolor','k');
hold off

xlabel('Initial fraction of the killer (x = f_0)','FontSize',14);
ylabel('Initial total population (y = N_0)','FontSize',14);
axis equal
h = colorbar();
clim([-430,270]);
colormap(centered);
% Set custom ticks for the colorbar
set(h, 'Ticks', [-400, -300, -200, -100, 0, 100, 200]);
% Set corresponding custom labels, converting negative values to positive for the blue region
set(h, 'TickLabels', {'400', '300', '200', '100', '0', '100', '200'});
%% 4E 
figure;
contourf(X,Y,IterMap_4E*4,60,'Edgecolor','none');
hold on
plot(ff1,nn1,'-.','LineWidth',2.7,'Color','m');
plot(ff2(1:indx),nn2(1:indx),'-.','LineWidth',2.7,'Color','m');
hold off
xlabel('Initial fraction of the killer (x = f_0)','FontSize',14);
ylabel('Initial total population (y = N_0)','FontSize',14);
axis equal
h = colorbar();
clim([-430,270]);
colormap(centered);
% Set custom ticks for the colorbar
set(h, 'Ticks', [-400, -300, -200, -100, 0, 100, 200]);
% Set corresponding custom labels, converting negative values to positive for the blue region
set(h, 'TickLabels', {'400', '300', '200', '100', '0', '100', '200'});

%% 4F
figure
contourf(X,Y,IterMap_4F*8,60,'Edgecolor','none');
hold on
plot(ff1,nn1,'-.','LineWidth',2.7,'Color','m');
plot(ff2(1:indx),nn2(1:indx),'-.','LineWidth',2.7,'Color','m');
hold off

xlabel('Initial fraction of the killer (x = f_0)','FontSize',14);
ylabel('Initial total population (y = N_0)','FontSize',14);
axis equal
h = colorbar();
clim([-430,270]);
colormap(centered);
% Set custom ticks for the colorbar
set(h, 'Ticks', [-400, -300, -200, -100, 0, 100, 200]);
% Set corresponding custom labels, converting negative values to positive for the blue region
set(h, 'TickLabels', {'400', '300', '200', '100', '0', '100', '200'});

end
