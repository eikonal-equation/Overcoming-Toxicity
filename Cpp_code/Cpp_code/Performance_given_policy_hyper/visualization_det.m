%% Read in data from files
clear all; clc;
%close all;
precision1 = 'double';
precision2 = 'uint8';

% Read in parameter data
filename_parameters = 'output/Perf_hyper_DomainParameters.dat';


% Paramter Values
uFile_params = fopen(filename_parameters);
u_params = fread(uFile_params, 11, precision1);

Nx = u_params(1) + 1
Ny = u_params(2) + 1
dx = u_params(3)
dy = u_params(4)
epsilon = u_params(5)
rks = u_params(6)
gamma = u_params(7)
% lamb = u_params(8)
% rho = u_params(9)
vic_thres = u_params(10)
defeat_thres = u_params(11)
fclose(uFile_params);

defeat_indx = floor(defeat_thres/dx)+1;
vic_indx = floor(vic_thres/dx)+1;

lamb = 0.75;
rho = 0.6;

filename2 = ['output/Perf_hyper_const_valuefn_rho',num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat'];

% Value function file
wFile = fopen(filename2);
w_ult = fread(wFile, precision1);
w_ult=  reshape(w_ult,[1601,1601]);
fclose(wFile);
%% Plotting
yy = linspace(0,1,Ny);
xx = linspace(0,1,Nx);

[X,Y] = meshgrid(xx,yy);

ff = @(x) 1./x;

figure;
% subplot(1,2,1)
contourf(X,Y,w_ult,15,'Edgecolor','none');
hold on
plot(xx, vic_thres.*ff(xx),'m:','linewidth',2)
plot(xx, defeat_thres.*ff(xx),'m:','linewidth',2)
hold off
xlabel('Initial fraction of the killer (x = f_0)','FontSize',14);
ylabel('Initial total population (y = N_0)','FontSize',14);
axis equal
% title('Value function','fontsize',16,'Interpreter','latex');
colorbar();
xlim([0,1]);
ylim([0,1]);

% subplot(1,2,2)
% contourf(X,Y,policy,[0 1]);
% hold on
% plot(xx, vic_thres.*ff(xx),'m:','linewidth',2)
% plot(xx, defeat_thres.*ff(xx),'m:','linewidth',2)
% hold off
% xlabel('Frequency (f)','FontSize',14);
% ylabel('Total population (N)','FontSize',14);
% axis equal
% title('Optimal Toxin-On/Off Region','fontsize',16,'Interpreter','latex');
% colorbar();
% xlim([0,1]);
% ylim([0,1]);

