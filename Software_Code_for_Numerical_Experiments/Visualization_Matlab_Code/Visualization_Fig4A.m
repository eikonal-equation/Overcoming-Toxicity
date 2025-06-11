function Visualization_Fig4A(File_valfn,File_parameters)
% Visualization_Fig4A.m
% This function generates the Fig. 4A in the paper, which visualizes
% the value function of the optimal control problem for a system with
% constitutive toxin production (Eq.[2] in the paper with "delta == 0" and "a == 1".)
% The visualization includes the contour plot of the value function and
% the region where the value function is less than the initial fraction of the killer.
% Input:
%   File_valfn: A string representing the path to the value function file
%   File_parameters: A string representing the path to the parameter file
%
% Author: MingYi Wang, Cornell University
% Last Modified: 06/2025
%
clear; clc;
%% Read in the data files
data_precision = 'double';
% Paramter Values
uFile_params = fopen(File_parameters);
u_params = fread(uFile_params,10,data_precision);

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

% Value function file
wFile = fopen(File_valfn);
valfn = fread(wFile, data_precision);
valfn  = reshape(valfn ,[Ny,Nx, fM]);
fclose(wFile);

%% Plotting
xx = linspace(0,1,Nx); % sample points in x-direction
yy = linspace(0,1,Ny); % sample points in y-direction
[X,Y] = meshgrid(xx,yy); % create a meshgrid for plotting

slice_num = floor(Tf/dt) + 1; % Determine the slice number based on the total time and time step size
valfn_slice = valfn(:,:,slice_num);

uu = zeros(Ny,Nx); % Initialize a matrix for the region that f(T^-) < f_0
uu(valfn_slice < X) = 1; % Set the region where f(T^-) < f_0 to 1

save("Fslice_4A.mat","uu");
save("Fslice_valfn4A.mat","valfn_slice");
%% Plotting the value function
figure
contourf(X,Y,valfn_slice,12,'Edgecolor','none');
hold on
contour(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),uu(2:end-1,2:end-1), [1], 'k--','LineWidth',2);
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

end
