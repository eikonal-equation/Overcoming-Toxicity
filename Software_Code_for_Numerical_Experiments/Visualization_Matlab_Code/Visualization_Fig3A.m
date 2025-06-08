function Visualization_Fig3A()
% Visualization_Fig3A.m
% This function generates the Fig. 3A in the paper, which visualizes
% the population dynamics of a system with
% constitutive toxin production (Eq.[2] in the paper with "a == 1".)
% The visualization includes the vector field of the population dynamics
% and the trajectories of the system under different initial conditions.
%
% Author: MingYi Wang, Cornell University
% Last Modified: 06/2025
%
%% Parameters
epsilon = 0.2; %The cost of constitutive toxin production
gamma = 1.0; %The rescaled cost of toxin-production rate
rks = 0.85; %The rescaled ratio of the growth rate of the killer to that of the sensitive
d = 0.4308; % rescaled basal death rate (== -log(rho))

% Drift part of the population growth dynamics (Eq.[2] in the paper with "delta == 0")
fx = @(x,y,a)  x.*(1-x).*((1-y).*(rks.*(1-epsilon.*a) - 1) + a.*gamma.*x.*y);
fy = @(x,y,a)  y.*(1-y).*(1 + (rks.*(1-epsilon.*a) - 1).*x) - a.*gamma.*y.^2.*x.*(1-x) - d.*y;

%% plot on qp-plane
a=1; % Fix the toxin production rate to 1
% Exact duplication of the initial conditions we used in the paper
traj_1 = traj2(0.95,0.05,a,fx,fy,epsilon,gamma,rks);
traj_2 = traj2(0.3,0.1,a,fx,fy,epsilon,gamma,rks);
traj_3 = traj2(0.5,0.1,a,fx,fy,epsilon,gamma,rks);
traj_4 = traj2(0.12,0.05,a,fx,fy,epsilon,gamma,rks);
traj_5 = traj(0.3,0.97,a,fx,fy,epsilon,gamma,rks);
traj_6 = traj(0.6,0.96,a,fx,fy,epsilon,gamma,rks);
traj_7 = traj(0.7,0.08,a,fx,fy,epsilon,gamma,rks);
traj_8 = traj(0.15,0.97,a,fx,fy,epsilon,gamma,rks);
traj_9 = traj(0.85,0.96,a,fx,fy,epsilon,gamma,rks);

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

x_traj_8 = traj_8(1,:);
y_traj_8 = traj_8(2,:);

x_traj_9 = traj_9(1,:);
y_traj_9 = traj_9(2,:);

% Compute the vector field
N=20; % Number of grid points we intended to plot the vector field
NN = 800; % Number of grid points we used to compute the vector field
eps = 1e-6; % Numerical tolerance for zero
% Discretization
x_plot=linspace(0,1,N+1);
y_plot=linspace(0,1,N+1);
x_vecf=linspace(0,1,NN+1);
y_vecf=linspace(0,1,NN+1);
[x,y]=meshgrid(x_plot',y_plot');
[xx,yy]=meshgrid(x_vecf',y_vecf');
% Compute the drift
df=fs(x,y,a); dN=fk(x,y,a);
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
plot(x_traj_1,y_traj_1,'-','LineWidth',2.7,'Color','k');
plot(x_traj_2,y_traj_2,'-','LineWidth',2.7,'Color','k');
plot(x_traj_3,y_traj_3,'-','LineWidth',2.7,'Color','k');
plot(x_traj_4,y_traj_4,'-','LineWidth',2.7,'Color','k');
plot(x_traj_5,y_traj_5,'-','LineWidth',2.7,'Color','k');
plot(x_traj_6,y_traj_6,'-','LineWidth',2.7,'Color','k');
plot(x_traj_7,y_traj_7,'-','LineWidth',2.7,'Color','k');
plot(x_traj_8,y_traj_8,'-','LineWidth',2.7,'Color','k');
plot(x_traj_9,y_traj_9,'-','LineWidth',2.7,'Color','k');


load('separatrix_4.mat','y_unstable_backward',"y_stable_backward",'y_unstable_forward',"y_stable_forward");
% Plot the stable and unstable manifolds
nn1 = y_stable_forward(:,1) + y_stable_forward(:,2);
ff1 = y_stable_forward(:,1)./nn1;
ff1(ff1>1) = 1; 
nn2 = y_stable_backward(:,1) + y_stable_backward(:,2);
ff2 = y_stable_backward(:,1)./nn2;
ff2(ff2>1) = 1;
ff2(ff2<1e-3) = 1;
nn2(nn2<1e-3) = 0;
plot(ff1,nn1,'-.','LineWidth',2.7,'Color','m');
plot(ff2(1:59),nn2(1:59),'-.','LineWidth',2.7,'Color','m');
plot(1,pt1,'o','MarkerEdgeColor','k','MarkerFaceColor','r','markersize',6.5)
plot(0,pt2,'o','MarkerEdgeColor','k','MarkerFaceColor','b','markersize',6.5)
% Select points along the trajectory for arrows
arrow_index1 = 12; 
addingArrow(arrow_index1,x_11,y_11);
arrow_index1 = 10; 
addingArrow(arrow_index1,x_22,y_22);
arrow_index1 = 10; 
addingArrow(arrow_index1,x_33,y_33);
arrow_index1 = 12; 
addingArrow(arrow_index1,x_44,y_44);
arrow_index1 = 12; 
addingArrow(arrow_index1,x_55,y_55);
arrow_index1 = 12; 
addingArrow(arrow_index1,x_66,y_66);
arrow_index1 = 12; 
addingArrow(arrow_index1,x_77,y_77);
arrow_index1 = 12; 
addingArrow(arrow_index1,x_88,y_88);
arrow_index1 = 12; 
addingArrow(arrow_index1,x_99,y_99);
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

function addingArrow(arrow_index1,x,y)
    % This function adds an arrow to the plot at the specified index
    % Inputs:
    % arrow_index1: index of the arrow
    % x: x-coordinates of the trajectory    
    % y: y-coordinates of the trajectory


    % Calculate the direction of the arrows
    dx1 = x(arrow_index1 + 1) - x(arrow_index1);
    dy1 = y(arrow_index1 + 1) - y(arrow_index1);
    % Plot the arrows using quiver
    ah = annotation('arrow',...
        'Color', 'k',...
        'headStyle','vback2','HeadLength',7,'HeadWidth',11,'LineStyle','none');
    set(ah,'parent',gca);
    set(ah,'position',[x(arrow_index1) y(arrow_index1) 0.04*dx1 0.04*dy1]);

end

end
