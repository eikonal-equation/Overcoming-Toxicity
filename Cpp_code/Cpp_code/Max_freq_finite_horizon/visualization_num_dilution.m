%% Read in data from files
clear all; clc;
%close all;
precision1 = 'double';
precision2 = 'uint8';
filename3 = 'Max_freq_dilution_f.dat';
filename4 = 'Max_freq_dilution_N.dat';
% Read in parameter data
filename_parameters = 'Max_freq_finite_horizon_single_DomainParameters.dat';


% Paramter Values
uFile_params = fopen(filename_parameters);
u_params = fread(uFile_params,11, precision1);

Nx = u_params(1) + 1
Ny = u_params(2) + 1
dx = u_params(3)
dy = u_params(4)
epsilon = u_params(5)
rks = u_params(6)
gamma = u_params(7)
Tf = u_params(8)
fM = u_params(9)+1
dt = u_params(10)
num_dilution = u_params(11)
fclose(uFile_params);
stroage = (num_dilution - mod(num_dilution,10))/10+2;

wFile = fopen(filename3);
f = fread(wFile, precision1);
f =  reshape(f,[Ny,Nx, stroage]);
fclose(wFile);

wFile = fopen(filename4);
N = fread(wFile, precision1);
N =  reshape(N,[Ny,Nx, stroage]);
fclose(wFile);
%% Plotting
yy = linspace(0,1,Ny);
xx = linspace(0,1,Nx);

[X,Y] = meshgrid(xx,yy);
uu = zeros(Ny,Nx);

slice_num = 10;

% w_slice = w(:,:,slice_num);
% policy_slice = policy(:,:,slice_num);
% uu(w_slice < X) = 1;

f_slice1 = f(:,:,slice_num);
N_slice1 = N(:,:,slice_num);


figure;
subplot(1,2,1)
contourf(X,Y,f_slice1,'Edgecolor','none');
xlabel('Fraction of the Killer (f)','FontSize',14);
ylabel('Total population (N)','FontSize',14);
axis equal
title('Relative abundance of killer in the limit','fontsize',16,'Interpreter','latex');
colorbar();

subplot(1,2,2)
contourf(X,Y,(1-f_slice1).*N_slice1,'Edgecolor','none');
xlabel('Fraction of the Killer (f)','FontSize',14);
ylabel('Total population (N)','FontSize',14);
axis equal
title('Sensitive population in the limit','fontsize',16,'Interpreter','latex');
colorbar();
