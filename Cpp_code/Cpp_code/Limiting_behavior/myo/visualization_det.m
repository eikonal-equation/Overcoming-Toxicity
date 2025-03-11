%% Read in data from files
clear all; clc;
%close all;
precision1 = 'double';
precision2 = 'int32';
precision3 = 'uint8';
filename1 = 'limit_iter_rho650_T1.dat';
filename2 = 'limit_f_rho650_T1.dat';
filename3 = 'limit_check_rho650_T1.dat';

filename11 = 'limit_iter_rho650_T2.dat';
% Read in parameter data
filename_parameters = 'Max_freq_finite_horizon_const_DomainParameters.dat';


% Paramter Values
uFile_params = fopen(filename_parameters);
u_params = fread(uFile_params,10, precision1);

Nx = u_params(1) + 1
Ny = u_params(2) + 1
dx = u_params(3)
dy = u_params(4)
epsilon = u_params(5)
rks = u_params(6)
gamma = u_params(7)
Tf = u_params(8)
fM = u_params(9) + 1
dt = u_params(10)
fclose(uFile_params);


% Value function file
wFile = fopen(filename1);
IterMap = fread(wFile, precision2);
IterMap =  reshape(IterMap,[Ny,Nx]);
fclose(wFile);
% % Value function file
% wFile = fopen(filename2);
% w = fread(wFile, precision1);
% w =  reshape(w,[Ny,Nx]);
% fclose(wFile);
% % Value function file
% wFile = fopen(filename3);
% cMat = fread(wFile, precision3);
% cMat =  reshape(cMat,[Ny,Nx]);
% fclose(wFile);

% Value function file
wFile = fopen(filename11);
IterMap2 = fread(wFile, precision2);
IterMap2 =  reshape(IterMap2,[Ny,Nx]);
fclose(wFile);

load('fslice.mat');
load('const slice.mat');
%% Plotting
yy = linspace(0,1,Ny);
xx = linspace(0,1,Nx);

[X,Y] = meshgrid(xx,yy);

figure;
% tiledlayout(1,2)
% ax1 = nexttile;
contourf(X,Y,IterMap,60,'Edgecolor','none');
hold on
contour(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),ff(2:end-1,2:end-1), [1], 'k--','LineWidth',2);
hold off

xlabel('Initial fraction of the killer (x = f_0)','FontSize',14);
ylabel('Initial total population (y = N_0)','FontSize',14);
axis equal
% title('$T=1$','fontsize',20,'Interpreter','latex');
h = colorbar();
colormap(centered);

% Set custom ticks for the colorbar
set(h, 'Ticks', [-300, -200, -100, 0, 100, 200]);

% Set corresponding custom labels, converting negative values to positive for the blue region
set(h, 'TickLabels', {'300', '200', '100', '0', '100', '200'});

% ax2 = nexttile;
% contourf(X,Y,IterMap2*2,50,'Edgecolor','none');
% xlabel('Initial fraction of the killer (x)','FontSize',14);
% ylabel('Initial total population (y)','FontSize',14);
% axis equal
% title('$T=2$','fontsize',20,'Interpreter','latex');
% colorbar();
% colormap(ax2,centered);

% figure
% contourf(X,Y,cMat,'Edgecolor','none');
% xlabel('Initial fraction of the killer (x)','FontSize',14);
% ylabel('Initial total population (y)','FontSize',14);
% axis equal
% % title('Value function','fontsize',16,'Interpreter','latex');
% colorbar();

