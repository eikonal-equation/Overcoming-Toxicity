function Visualization_Fig1C()
% Visualization_Fig1C.m
% This function generates the Fig. 1C in the paper, which visualizes
% the population dynamics of a system with
% constitutive toxin production (Eq.[2] in the paper with "delta == 0" and "a == 1".)
% The visualization includes the direction field of the population dynamics
% and the trajectories of the system under different initial conditions.
%
% Author: MingYi Wang, Cornell University
% Last Modified: 06/2025
%
%% Parameters
epsilon = 0.2; %The cost of constitutive toxin production
gamma = 1.0; %The rescaled cost of toxin-production rate
rks = 0.85; %The rescaled ratio of the growth rate of the killer to that of the sensitive
d = 0; % rescaled basal death rate 

% Drift part of the population growth dynamics (Eq.[2] in the paper with "delta == 0")
fx = @(x,y,a)  x.*(1-x).*((1-y).*(rks.*(1-epsilon.*a) - 1) + a.*gamma.*x.*y);
fy = @(x,y,a)  y.*(1-y).*(1 + (rks.*(1-epsilon.*a) - 1).*x) - a.*gamma.*y.^2.*x.*(1-x) - d.*y;;

%% First, we generate the trajectories
a=1; % Fix the toxin production rate to 1
% Exact duplication of the initial conditions we used in the paper
traj_1 = traj2(0.3,0.7,a,fk,fs,epsilon,gamma,rks);
traj_2 = traj2(0.2,0.1,a,fk,fs,epsilon,gamma,rks);
traj_3 = traj2(0.5,0.05,a,fk,fs,epsilon,gamma,rks);
traj_4 = traj2(0.05,0.05,a,fk,fs,epsilon,gamma,rks);
traj_5 = traj2(0.05,0.4,a,fk,fs,epsilon,gamma,rks);
traj_6 = traj2(0.6,0.4,a,fk,fs,epsilon,gamma,rks);
traj_7 = traj2(0.7,0.08,a,fx,fy,epsilon,gamma,rks);

% Extract the x and y coordinates of the points in the trajectories
x_traj_1 = traj_1(1,:);
y_traj_1 = traj_1(2,:);

x_traj_2 = traj_2(1,:);
y_traj_2 = traj_2(2,:);

x_traj_3 = traj_3(1,:);
y_traj_3 = traj_3(2,:);

x_traj_4 = traj_4(1,:);
y_traj_4 = traj_4(2,:);

x_traj_5 = traj_5(1,:);
y_traj_5 = traj_5(2,:);

x_traj_6 = traj_6(1,:);
y_traj_6 = traj_6(2,:);

x_traj_7 = traj_7(1,:);
y_traj_7 = traj_7(2,:);

% Compute the direction field
N=20; % Number of grid points we intended to plot the direction field
NN = 800; % Number of grid points we used to compute the direction field
eps = 1e-6; % Numerical tolerance for zero
% Discretization
x_plot=linspace(0,1,N+1);
y_plot=linspace(0,1,N+1);
x_vecf=linspace(0,1,NN+1);
y_vecf=linspace(0,1,NN+1);
[x,y]=meshgrid(x_plot',y_plot');
[xx,yy]=meshgrid(x_vecf',y_vecf');
% Compute the drift
df=fx(x,y,a); dN=fy(x,y,a);
% Normalize the drift
dN(abs(dN-0) < eps) = 0;
df(abs(df-0) < eps) = 0;
rr=sqrt(df.^2+dN.^2);
u=df./rr; v=dN./rr;

% Background matrix for visualization
Zmat = ones(NN+1);
Zmat(1,1) = 0; % For visualization purpose, we set the first element to 0

%% Plotting 
figure
contourf(xx,yy,Zmat);
hold on
for ii = 2:20
    for ij = [1:1:N+1]
        if isnan(u(ii,ij)) || isnan(v(ii,ij))
            continue
        else
        ah = annotation('arrow',...
            'Color', [169/255,169/255,169/255],...
            'headStyle','cback1','HeadLength',3.5,'HeadWidth',3.5);
        set(ah,'parent',gca);
        set(ah,'position',[x(ii,ij) y(ii,ij) 0.04*u(ii,ij) 0.04*v(ii,ij)]);
        end
    end
end
hold on 
plot(x_traj_1,y_traj_2,'-','LineWidth',2.7,'Color','k');
plot(x_traj_2,y_traj_2,'-','LineWidth',2.7,'Color','k');
plot(x_traj_3,y_traj_3,'-','LineWidth',2.7,'Color','k');
plot(x_traj_4,y_traj_4,'-','LineWidth',2.7,'Color','k');
plot(x_traj_5,y_traj_5,'-','LineWidth',2.7,'Color','k');
plot(x_traj_6,y_traj_6,'-','LineWidth',2.7,'Color','k');
plot(x_traj_7,y_traj_7,'-','LineWidth',2.7,'Color','k');
plot(1,1,'o','MarkerEdgeColor','k','MarkerFaceColor','r','markersize',6.5)
% Select points along the trajectory for arrows
arrow_index1 = 12; 
arrow_index2 = 18; 
addingArrow2(arrow_index1,arrow_index2,x_traj_1,y_traj_1,'k');
arrow_index1 = 10; 
arrow_index2 = 16; 
addingArrow2(arrow_index1,arrow_index2,x_traj_2,y_traj_2,'k');
arrow_index1 = 8; 
arrow_index2 = 16; 
arrow_index3 = 31; 
addingArrow3(arrow_index1,arrow_index2,arrow_index3,x_traj_3,y_traj_3,'k');
arrow_index1 = 12; 
arrow_index2 = 18; 
arrow_index3 = 108; 
addingArrow3(arrow_index1,arrow_index2,arrow_index3,x_traj_4,y_traj_4,'k');
arrow_index1 = 5; 
arrow_index2 = 15; 
addingArrow2(arrow_index1,arrow_index2,x_traj_5,y_traj_5,'k');
arrow_index1 = 3; 
arrow_index2 = 8; 
addingArrow2(arrow_index1,arrow_index2,x_traj_6,y_traj_6,'k');
arrow_index1 = 10; 
arrow_index2 = 16; 
addingArrow2(arrow_index1,arrow_index2,x_traj_7,y_traj_7,'k');

hold off
axis([0 1 0 1])
axis equal
shading flat
xlabel('Fraction of the killer (f)','FontSize',14);
ylabel('Total population (N)','FontSize',14);
set(findall(gca,'type','text'),'visible','on');
colormap("copper")

%% subfunctions
function full_traj = traj2(f0,N0,a,fx,fy,epsilon,gamma,rks)
    % This function computes the trajectory of the system
    % using the ODE45 solver.
    % Inputs:
    % f0: initial fraction of the killer
    % N0: initial population size
    % a: the toxin production rate
    % fx: the "f" drift part of the population growth model
    % fy: the "N" drift part of the population growth model
    % epsilon: the cost of constitutive toxin production
    % gamma: the rescaled cost of toxin-production rate
    % rks: the rescaled ratio of the growth rate of the killer to that of the sensitive
    % Outputs:
    % full_traj: the trajectory of the system
    tspan = [0 500];
    [T,Z] = ode45((@(t,y)[fx(y(1),y(2),a);fy(y(1),y(2),a)]),...
        tspan,[f0 N0]);
    fval=Z(:,1);
    Nval=Z(:,2);
    full_traj  = [fval';Nval'];
end

function addingArrow2(arrow_index1,arrow_index2,x,y,mycolor)
    % This function adds TWO arrows to the plot at specified indices
    % Inputs:
    % arrow_index1: index of the first arrow
    % arrow_index2: index of the second arrow
    % x: x-coordinates of the trajectory
    % y: y-coordinates of the trajectory
    % mycolor: color of the arrows
    
    % Calculate the direction of the arrows
    dx1 = x(arrow_index1 + 1) - x(arrow_index1);
    dy1 = y(arrow_index1 + 1) - y(arrow_index1);

    dx2 = x(arrow_index2 + 1) - x(arrow_index2);
    dy2 = y(arrow_index2 + 1) - y(arrow_index2);


    % Plot the arrows using quiver
    ah = annotation('arrow',...
        'Color', mycolor,...
        'headStyle','vback2','HeadLength',7,'HeadWidth',11,'LineStyle','none');
    set(ah,'parent',gca);
    set(ah,'position',[x(arrow_index1) y(arrow_index1) 0.04*dx1 0.04*dy1]);
    ah = annotation('arrow',...
        'Color', mycolor,...
        'headStyle','vback2','HeadLength',7,'HeadWidth',11,'LineStyle','none');
    set(ah,'parent',gca);
    set(ah,'position',[x(arrow_index2) y(arrow_index2) 0.04*dx2 0.04*dy2]);
end


function addingArrow3(arrow_index1,arrow_index2,arrow_index3,x,y,mycolor)
    % This function adds THREE arrows to the plot at specified indices
    % Inputs:
    % arrow_index1: index of the first arrow
    % arrow_index2: index of the second arrow
    % arrow_index3: index of the third arrow
    % x: x-coordinates of the trajectory
    % y: y-coordinates of the trajectory
    % mycolor: color of the arrows

    % Calculate the direction of the arrows
    dx1 = x(arrow_index1 + 1) - x(arrow_index1);
    dy1 = y(arrow_index1 + 1) - y(arrow_index1);

    dx2 = x(arrow_index2 + 1) - x(arrow_index2);
    dy2 = y(arrow_index2 + 1) - y(arrow_index2);

    dx3 = x(arrow_index3 + 1) - x(arrow_index3);
    dy3 = y(arrow_index3 + 1) - y(arrow_index3);
    % Plot the arrows using quiver
    ah = annotation('arrow',...
        'Color', mycolor,...
        'headStyle','vback2','HeadLength',7,'HeadWidth',11,'LineStyle','none');
    set(ah,'parent',gca);
    set(ah,'position',[x(arrow_index1) y(arrow_index1) 0.04*dx1 0.04*dy1]);
    ah = annotation('arrow',...
        'Color', mycolor,...
        'headStyle','vback2','HeadLength',7,'HeadWidth',11,'LineStyle','none');
    set(ah,'parent',gca);
    set(ah,'position',[x(arrow_index2) y(arrow_index2) 0.04*dx2 0.04*dy2]);
    ah = annotation('arrow',...
        'Color', mycolor,...
        'headStyle','vback2','HeadLength',7,'HeadWidth',11,'LineStyle','none');
    set(ah,'parent',gca);
    set(ah,'position',[x(arrow_index3) y(arrow_index3) 0.04*dx3 0.04*dy3]);
end

end
