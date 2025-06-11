function Visualization_Fig5(File_parameters)
% Visualization_Fig5.m
% This function generates the Fig. 5 in the paper, which visualizes
% the improvement of the value function
% from the constitutive toxin production case (Fig. 4A) to the myopic
% toxin production case (Fig. 5A).
% The visualization includes the contour plots of the value function 
% improvement, the optimal policy, and the time until competitive exclusion.
% Input:
%   File_parameters: A string representing the path to the parameter file
%
% Author: MingYi Wang, Cornell University
% Last Modified: 06/2025
%
clear; clc;
%% Read in the data files
data_precision1 = 'double';
data_precision2 = 'uint8';
filename_valfn = 'Regular_Myopic_opt_valuefn.dat';
filename_policy = 'Regular_Myopic_opt_policy.dat';
% Read in parameter data
uFile_params = fopen(File_parameters);
u_params = fread(uFile_params,10, data_precision1);
Nx = u_params(1) + 1 % Number of grid points in x-direction
Ny = u_params(2) + 1 % Number of grid points in y-direction
dx = u_params(3) % Grid spacing in x-direction
dy = u_params(4) % Grid spacing in y-direction
epsilon = u_params(5) % Cost of constitutive toxin production
rks = u_params(6) % Rescaled ratio of the growth rate of the killer to that of the sensitive
gamma = u_params(7) % Rescaled cost of toxin-production rate
Tf = u_params(8) % Regular inter-dilution time
fM = u_params(9)+1 % Number of time slices
dt = u_params(10) % Time step size
fclose(uFile_params);

% Drift part of the population growth dynamics (Eq.[2] in the paper with "delta == 0")
fx = @(x,y,a)  x.*(1-x).*((1-y).*(rks.*(1-epsilon.*a) - 1) + a.*gamma.*x.*y);
fy = @(x,y,a)  y.*(1-y).*(1 + (rks.*(1-epsilon.*a) - 1).*x) - a.*gamma.*y.^2.*x.*(1-x);

% Value function file
wFile = fopen(filename_valfn);
valfn = fread(wFile, data_precision1);
valfn =  reshape(valfn,[Ny,Nx, fM]);
fclose(wFile);
% Optimal policy file
aFile = fopen(filename_policy);
policy = fread(aFile, data_precision2);
policy = reshape(policy,[Ny,Nx, fM]);
fclose(aFile);
% Load the slice of the value function from the constitutive case from Fig. 4A
load("Fslice_valfn4A.mat","valfn_slice");
valfn_slice_4A = valfn_slice;
%% Plotting initalization
xx = linspace(0,1,Nx); % sample points in x-direction
yy = linspace(0,1,Ny); % sample points in y-direction
[X,Y] = meshgrid(xx,yy); % create a meshgrid for plotting
slice_num = floor(Tf/dt) + 1; % Determine the slice number based on the total time and time step size
% Extract the slice of the value function and policy
valfn_slice = valfn(:,:,slice_num);
policy_slice = policy(:,:,slice_num);

ww = zeros(Ny,Nx); % Initialize a matrix for the region that f(T^-) < f_0
ww(valfn_slice < X) = 1; % Set the region where f(T^-) < f_0 to 1

save("Fslice_5A.mat","ww");
%% Fig. 5A: Plotting improvement from constitutive case
figure;
contourf(X,Y, valfn_slice - valfn_slice_4A,'Edgecolor','none');
xlabel('Initial fraction of the killer (x = f_0)','FontSize',14);
ylabel('Initial total population (y = N_0)','FontSize',14);
axis equal
colorbar();

%% Fig. 5C: Plotting the myopic-optimal policy at the 0-th slice
N=20; % Number of grid points we intended to plot the vector field
NN = 1600; % Number of grid points we used to compute the vector field
eps = 1e-6; % Numerical tolerance for zero
x_plot=linspace(0,1,N+1);
y_plot=linspace(0,1,N+1);
[x,y]=meshgrid(x_plot',y_plot');
policy_slice_downsample = policy_slice(1:NN/N:end,1:NN/N:end);
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
contourf(X,Y,policy_slice, [0, 1],'EdgeColor','none');
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
            if policy_slice_downsample(ii,ij) == 0
                set(ah,'position',[x(ii,ij) y(ii,ij) 0.03*u_a0(ii,ij) 0.03*v_a0(ii,ij)]);
            else
                set(ah,'position',[x(ii,ij) y(ii,ij) 0.03*u_a1(ii,ij) 0.03*v_a1(ii,ij)]);
            end
        end
    end
end
hold off
xlabel('Initial fraction of the killer (x = f_0)','FontSize',14);
ylabel('Initial total population (y = N_0)','FontSize',14);
axis equal
colormap("copper")

%% Fig. 5B: Plotting the time until competitive exclusion
% Read in the data files
load("const_bd.mat","const_bd");
filename_itermap = 'Regular_myopic_opt_limit_iter_rho650_T1.dat';
data_precision3 = 'int32';
% Itermap file
itFile = fopen(filename_itermap);
Itermap = fread(itFile, data_precision3);
Itermap =  reshape(Itermap,[Ny,Nx, fM]);
fclose(itFile);

myopic_bd = Itermap;
myopic_bd(Itermap < 0) = 0; % Set the cost boundary
myopic_bd(Itermap >= 0) = 1; % Set the cost boundary
save("mypoic_bd.mat","myopic_bd");

% Plotting the time until competitive exclusion
figure;
contourf(X,Y,Itermap,60,'Edgecolor','none');
hold on
contour(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),const_bd(2:end-1,2:end-1), [1], 'k--','LineWidth',2);
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

end
