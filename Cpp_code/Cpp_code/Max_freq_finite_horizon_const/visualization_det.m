%% Read in data from files
clear all; clc;
%close all;
precision1 = 'double';
precision2 = 'uint8';
filename1 = 'Max_freq_finite_horizon_const_valuefn.dat';
% filename2 = 'Max_freq_finite_horizon_policy.dat';
filename3 = 'Max_freq_finite_horizon_inf_limit_f_const.dat';
filename4 = 'Max_freq_finite_horizon_inf_limit_N_const.dat';
% Read in parameter data
filename_parameters = 'Max_freq_finite_horizon_const_DomainParameters.dat';


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
num_dilution = u_params(11)+1
fclose(uFile_params);

% fM = 2

% Value function file
wFile = fopen(filename1);
w = fread(wFile, precision1);
w =  reshape(w,[Ny,Nx, fM]);
fclose(wFile);
% % Optimal policy file
% aFile = fopen(filename2);
% policy = fread(aFile, precision2);
% policy = reshape(policy,[Ny,Nx, fM]);
% fclose(aFile);
% Value function file
wFile = fopen(filename3);
f = fread(wFile, precision1);
f =  reshape(f,[Ny,Nx, fM]);
fclose(wFile);

wFile = fopen(filename4);
N = fread(wFile, precision1);
N =  reshape(N,[Ny,Nx, fM]);
fclose(wFile);
%% Plotting
yy = linspace(0,1,Ny);
xx = linspace(0,1,Nx);

[X,Y] = meshgrid(xx,yy);
uu = zeros(Ny,Nx);

slice_num = 161;

w_slice = w(:,:,slice_num);
% policy_slice = policy(:,:,slice_num);
uu(w_slice < X) = 1;

f_slice1 = f(:,:,slice_num);
N_slice1 = N(:,:,slice_num);

figure;
% subplot(1,2,1)
contourf(X,Y,w_slice,12,'Edgecolor','none');
hold on
contour(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),uu(2:end-1,2:end-1), [1], 'k--','LineWidth',2);
% plot(0.35,0.7,'o','markersize',5.5,'markerfacecolor','c','markeredgecolor','k');
% plot(0.5,0.7,'o','markersize',5.5,'markerfacecolor','m','markeredgecolor','k');
hold off

xlabel('Initial fraction of the killer (x = f_0)','FontSize',14);
ylabel('Initial total population (y = N_0)','FontSize',14);
axis equal
% title('Value function','fontsize',16,'Interpreter','latex');
hcb = colorbar();
clim([0 1]);
cc = linspace(0,1,10);
set(hcb,'Ticks',cc);
myCticks = arrayfun(@(x) sprintf('%.1f',x),cc,'un',0);
set(hcb,'Ticklabels',myCticks);
% 
% subplot(1,2,2)
% contourf(X,Y,policy_slice,[0 1]);
% xlabel('Fraction of the Killer (f)','FontSize',14);
% ylabel('Total population (N)','FontSize',14);
% axis equal
% title('Optimal Toxin-On/Off Region','fontsize',16,'Interpreter','latex');
% colorbar();
% f_slice2 = f_slice1;
% f_slice2(f_slice2 > 0.999999999) = 1;
% figure;
% subplot(1,2,1)
% contourf(X,Y,f_slice2,'Edgecolor','none');
% hold on
% contour(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),uu(2:end-1,2:end-1), [1], 'k--','LineWidth',2);
% plot(0.4,0.7,'o','markersize',5.5,'markerfacecolor','c','markeredgecolor','k');
% plot(0.5,0.7,'o','markersize',5.5,'markerfacecolor','m','markeredgecolor','k');
% hold off
% xlabel('Initial fraction of the killer (x)','FontSize',14);
% ylabel('Initial total population (y)','FontSize',14);
% axis equal
% title('Relative abundance of killer in the limit','fontsize',16,'Interpreter','latex');
% hcb = colorbar();
% clim([0,1]);
% cc = linspace(0,1,10);
% set(hcb,'Ticks',cc);
% myCticks = arrayfun(@(x) sprintf('%.1f',x),cc,'un',0);
% set(hcb,'Ticklabels',myCticks);
%%
figure
% subplot(1,2,2)
% contourf(X,Y,(1-f_slice1).*N_slice1,12,'Edgecolor','none');
contourf(X,Y,N_slice1,12,'Edgecolor','none');
xlabel('Fraction of the killer (f)','FontSize',14);
ylabel('Total population (N)','FontSize',14);
axis equal
% title('Sensitive population in the limit','fontsize',16,'Interpreter','latex');
hcb = colorbar();
clim([0.3 0.8]);
% cc = linspace(0,1,10);
% set(hcb,'Ticks',cc);
% myCticks = arrayfun(@(x) sprintf('%.1f',x),cc,'un',0);
% set(hcb,'Ticklabels',myCticks);

%%
%close all;
precision1 = 'double';
precision2 = 'uint8';
filename3 = 'Max_freq_finite_horizon_const_valuefn_lastslice.dat';
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

wFile = fopen(filename3);
w = fread(wFile, precision1);
w =  reshape(w,[Ny,Nx]);
fclose(wFile);

yy = linspace(0,1,Ny);
xx = linspace(0,1,Nx);

[X,Y] = meshgrid(xx,yy);
uu = zeros(Ny,Nx);


% policy_slice = policy(:,:,slice_num);
uu(w < X) = 1;


figure;
% subplot(1,2,1)
contourf(X,Y,w,12,'Edgecolor','none');
hold on
contour(X,Y,uu, [1], 'r--','LineWidth',2);
hold off
xlabel('Fraction of the killer (f)','FontSize',14);
ylabel('Total population (N)','FontSize',14);
axis equal
% title('Value function','fontsize',16,'Interpreter','latex');
hcb = colorbar();
clim([0 1]);
cc = linspace(0,1,10);
set(hcb,'Ticks',cc);
myCticks = arrayfun(@(x) sprintf('%.1f',x),cc,'un',0);
set(hcb,'Ticklabels',myCticks);


%%
% Extract the contour lines
C = contourc(X(1,:), Y(:,1), uu, [0.5 0.5]); % Use 0.5 as the level to extract the boundary

% Extract the x and y coordinates of the boundary
x_boundary = C(1, 2:end);
y_boundary = C(2, 2:end);

% Plot the boundary
figure
hold on;
plot(x_boundary, y_boundary, 'r', 'LineWidth', 2);
hold off;


% Extract the contour lines
C2 = contourc(X(1,:), Y(:,1), f_slice1, [0.5 0.5]); % Use 0.5 as the level to extract the boundary

% Extract the x and y coordinates of the boundary
x_boundary2 = C2(1, 2:end);
y_boundary2 = C2(2, 2:end);

% Plot the boundary
figure
hold on;
plot(x_boundary2, y_boundary2, 'r', 'LineWidth', 2);
hold off;

%%
aa = (x_boundary < 1e-3) | (x_boundary > 0.999);
x_boundary(aa) = [];
y_boundary(aa) = [];


figure
hold on
plot(x_boundary2, y_boundary2, 'LineWidth', 2);
plot(x_boundary, y_boundary, 'LineWidth', 2);
hold off



% Define the curves (x_boundary, y_boundary) and (x_boundary2, y_boundary2)
curve1 = [x_boundary; y_boundary];
curve2 = [x_boundary2; y_boundary2];

% Find the intersection points
[x0,y0,iout,jout] = intersections(x_boundary,y_boundary,x_boundary2,y_boundary2);

% Plot the curves and intersection points
plot(x_boundary, y_boundary, 'r', 'LineWidth', 2);
hold on;
plot(x_boundary2, y_boundary2, 'b', 'LineWidth', 2);
plot(x0, y0, 'go', 'MarkerSize', 10);
hold off;


x0 = x0(1);
y0 = y0(1);

kx = find(x0 <= xx',1);
ky = find(y0 <= xx',1);

xleft = xx(kx-1); xright = xx(kx);
ybot = xx(ky-1); ytop = xx(ky);
u11 = w(ky-1,kx-1);
u12 = w(ky-1,kx);
u21 = w(ky,kx-1);
u22 = w(ky,kx);
%bilinear interpolation
w_interp = BilinearInterp(u11,u12,u21,u22,xleft,xright,...
    ytop,ybot,x0,y0,dx,dy);


%%
TT = 1;
AA = ones(Ny,Nx,161);

choice = "conservative"
[xlist,ylist,policy_list] = path_finiteT(x0,y0,Nx-1,dt,TT,AA,choice,gamma,epsilon,rks);

rho = 0.65;
nafter = ylist(end)*rho

%%
figure
contourf(X,Y,f_slice1,12,'EdgeColor','none');
hold on
plot(x_boundary, y_boundary, 'c--','linewidth',2);
plot(x_boundary2, y_boundary2, 'c--','linewidth',2);
plot(xlist,ylist,'r-','LineWidth',2);
hold off
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
axis equal