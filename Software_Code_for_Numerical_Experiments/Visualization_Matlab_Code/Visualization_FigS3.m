function Visualization_FigS3()
% Visualization_FigS3.m
% This function generates the Fig. S3 in the paper, which visualizes
% the population dynamics of a system with constitutive toxin production
% (Eq.[S2.1] in the SI Appendix.)
% The visualization includes the direction field of the population dynamics
% and the trajectories of the system under different initial conditions
%
% Author: MingYi Wang, Cornell University
% Last Modified: 06/2025
%
clear; clc;
%%  parameters for delta == 0.2
epsilon = 0.2; % The cost of constitutive toxin production
gamma = 1.0; % The rescaled cost of toxin-production rate
rks = 0.85; % The rescaled ratio of the growth rate of the killer to that of the sensitive
d1 = 0.2; % rescaled basal death rate of the killer
d2 = 0.2; % rescaled basal death rate of the sensitive

% Drift part of the population growth dynamics (Eq.[S2.1] in the SI Appendix)
fk = @(nk,ns,a) rks.*(1-epsilon.*a).*nk.*(1-nk-ns) -d1*nk;
fs = @(nk,ns,a) ns.*(1-nk-ns)-a.*gamma.*nk.*ns - d2*ns;
% Fixed points of the system
pt1 = 1 - d1/(rks*(1-epsilon));
pt2 = 1 - d2;
pt3 = d1/(gamma*rks*(1-epsilon)) - d2/gamma;
pt4 = (d2 + gamma)/gamma - (d1*(1+gamma))/(gamma*rks*(1-epsilon));

% Generate the direction field on the (K,S)-triangle
a=1; % Fix the toxin production rate to 1
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
dNs=fs(x,y,a); dNk=fk(x,y,a);
% Normalize the drift
dNs(abs(dNs-0) < eps) = 0;
dNk(abs(dNk-0) < eps) = 0;
rr=sqrt(dNk.^2+dNs.^2);
u=dNk./rr; v=dNs./rr;
% Generate the direction field
U = ones(N+1)*NaN;
V = ones(N+1)*NaN;
% Background matrix for visualization
Zmat = ones(NN+1);
Zmat(1,1) = 0; % For visualization purpose, we set the first element to 0

for i = 1:N+1
    for j = 1:N+1
        if x_plot(i) + y_plot(j) <= 1
            U(i,j) = u(i,j);
            V(i,j) = v(i,j);
        else
            U(i,j) = NaN;
            V(i,j) = NaN;

        end
    end
end

for i = 1:NN+1
    for j = 1:NN+1
        if x_vecf(i) + y_vecf(j) > 1
            Zmat(i,j) = NaN;
        end
    end
end

% Plotting
figure
contourf(xx,yy,Zmat);
hold on
for ii = 2:20
    for ij = [1:1:N+1]
        if isnan(U(ii,ij)) || isnan(V(ii,ij))
            continue
        else
            ah = annotation('arrow',...
                'Color', [169/255,169/255,169/255],...
                'headStyle','cback1','HeadLength',3.5,'HeadWidth',3.5);
            set(ah,'parent',gca);
            set(ah,'position',[x(ii,ij) y(ii,ij) 0.04*U(ii,ij) 0.04*V(ii,ij)]);
        end
    end
end
hold on

for ii = 2:9
    % Generate a sequence of initial conditions
    xc = 0.12*ii;
    yc =  1 - xc;
    arrow_index1 = 8;
    Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks)
end
for jj = 1:4
    % Generate another sequence of initial conditions
    xc = 0.011*jj;
    yc =  0.01;
    arrow_index1 = 22 -jj;
    Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);
end
xc = 0.005; yc = 0.02;
arrow_index1 = 18;
Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);

xc = 0.06; yc = 0.006;
arrow_index1 = 17;
Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);

xc = 0.08; yc = 0.92;
arrow_index1 = 7;
Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);

xc = 0.03; yc = 0.97;
arrow_index1 = 6;
Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);

xc = 0.16; yc = 0.84;
arrow_index1 = 8;
Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);

% Plot the trajectories passing through the saddle point
hyperbolic_traj(pt3,pt4, a,fk,fs,epsilon,gamma,rks,0.2);
% Plot the fixed points
plot(0,0,'o','MarkerEdgeColor','k','MarkerFaceColor','b','markersize',6)
plot(pt1,0,'o','MarkerEdgeColor','k','MarkerFaceColor','b','markersize',6)
plot(0,pt2,'o','MarkerEdgeColor','k','MarkerFaceColor','b','markersize',6)
plot(pt3,pt4,'o','MarkerEdgeColor','k','MarkerFaceColor','r','markersize',6)
hold off

axis([0 1 0 1])
axis equal
shading flat

text(-0.1,0.5,'$n_S$','FontSize',15,'interpreter','latex');
text(0.45,-0.05,'$n_K$','FontSize',15,'Interpreter','latex');
set(gca,'visible','off');

set(findall(gca,'type','text'),'visible','on');
colormap("copper")

%%  parameters for delta == 0
epsilon = 0.2; % The cost of constitutive toxin production
gamma = 1.0; % The rescaled cost of toxin-production rate
rks = 0.85; % The rescaled ratio of the growth rate of the killer to that of the sensitive
d1 = 0; % rescaled basal death rate of the killer
d2 = 0; % rescaled basal death rate of the sensitive

% Drift part of the population growth dynamics (Eq.[S2.1] in the SI)
fk = @(nk,ns,a) rks.*(1-epsilon.*a).*nk.*(1-nk-ns) -d1*nk;
fs = @(nk,ns,a) ns.*(1-nk-ns)-a.*gamma.*nk.*ns - d2*ns;
% Fixed points of the system
pt1 = 1 - d1/(rks*(1-epsilon));
pt2 = 1 - d2;
pt3 = d1/(gamma*rks*(1-epsilon)) - d2/gamma;
pt4 = (d2 + gamma)/gamma - (d1*(1+gamma))/(gamma*rks*(1-epsilon));

% Generate the direction field on the (K,S)-triangle
a=1; % Fix the toxin production rate to 1
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
dNs=fs(x,y,a); dNk=fk(x,y,a);
% Normalize the drift
dNs(abs(dNs-0) < eps) = 0;
dNk(abs(dNk-0) < eps) = 0;
rr=sqrt(dNk.^2+dNs.^2);
u=dNk./rr; v=dNs./rr;
% Generate the direction field
U = ones(N+1)*NaN;
V = ones(N+1)*NaN;
% Background matrix for visualization
Zmat = ones(NN+1);
Zmat(1,1) = 0; % For visualization purpose, we set the first element to 0

for i = 1:N+1
    for j = 1:N+1
        if x_plot(i) + y_plot(j) <= 1
            U(i,j) = u(i,j);
            V(i,j) = v(i,j);
        else
            U(i,j) = NaN;
            V(i,j) = NaN;

        end
    end
end

for i = 1:NN+1
    for j = 1:NN+1
        if x_vecf(i) + y_vecf(j) > 1
            Zmat(i,j) = NaN;
        end
    end
end
%
figure
contourf(xx,yy,Zmat);
hold on
for ii = 2:20
    for ij = [1:1:N+1]
        if isnan(U(ii,ij)) || isnan(V(ii,ij))
            continue
        else
            ah = annotation('arrow',...
                'Color', [169/255,169/255,169/255],...
                'headStyle','cback1','HeadLength',3.5,'HeadWidth',3.5);
            set(ah,'parent',gca);
            set(ah,'position',[x(ii,ij) y(ii,ij) 0.04*U(ii,ij) 0.04*V(ii,ij)]);
        end
    end
end
hold on
for ii = 1:9
    xc = 0.1*ii;
    yc =  1 - xc;
    arrow_index1 = 4 + ii-1;
    Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);
end
for jj = 1:4
    xc = 0.005*jj;
    yc =  0.01;
    arrow_index1 = 22 -jj;
    Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);
end

for kk = 2:3
    xc = 0.025*kk;
    yc =  0.02;
    arrow_index1 = 18;
   Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);

end
xc = 0.002; yc = 0.02;
arrow_index1 = 22 -jj;
arrow_index2 = 24;
Add_trajectory2(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);

% Plot the fixed points
plot(0,0,'o','MarkerEdgeColor','k','MarkerFaceColor','b','markersize',6)
plot(pt1,0,'o','MarkerEdgeColor','k','MarkerFaceColor','b','markersize',6)
plot(0,pt2,'o','MarkerEdgeColor','k','MarkerFaceColor','b','markersize',6)
plot(pt3,pt4,'o','MarkerEdgeColor','k','MarkerFaceColor','b','markersize',6)

hold off
axis([0 1 0 1])
axis equal
shading flat

text(-0.1,0.5,'$n_S$','FontSize',15,'interpreter','latex');
text(0.45,-0.05,'$n_K$','FontSize',15,'Interpreter','latex');
set(gca,'visible','off');

set(findall(gca,'type','text'),'visible','on');
colormap("copper")


%%  parameters for delta == 0.0257
epsilon = 0.2; % The cost of constitutive toxin production
gamma = 1.0; % The rescaled cost of toxin-production rate
rks = 0.85; % The rescaled ratio of the growth rate of the killer to that of the sensitive
d1 = 0.0257; % rescaled basal death rate of the killer
d2 = 0.0257; % rescaled basal death rate of the sensitive

% Drift part of the population growth dynamics (Eq.[S2.1] in the SI)
fk = @(nk,ns,a) rks.*(1-epsilon.*a).*nk.*(1-nk-ns) -d1*nk;
fs = @(nk,ns,a) ns.*(1-nk-ns)-a.*gamma.*nk.*ns - d2*ns;
% Fixed points of the system
pt1 = 1 - d1/(rks*(1-epsilon));
pt2 = 1 - d2;
pt3 = d1/(gamma*rks*(1-epsilon)) - d2/gamma;
pt4 = (d2 + gamma)/gamma - (d1*(1+gamma))/(gamma*rks*(1-epsilon));

% Generate the direction field on the (K,S)-triangle
a=1; % Fix the toxin production rate to 1
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
dNs=fs(x,y,a); dNk=fk(x,y,a);
% Normalize the drift
dNs(abs(dNs-0) < eps) = 0;
dNk(abs(dNk-0) < eps) = 0;
rr=sqrt(dNk.^2+dNs.^2);
u=dNk./rr; v=dNs./rr;
% Generate the direction field
U = ones(N+1)*NaN;
V = ones(N+1)*NaN;
% Background matrix for visualization
Zmat = ones(NN+1);
Zmat(1,1) = 0; % For visualization purpose, we set the first element to 0

for i = 1:N+1
    for j = 1:N+1
        if x_plot(i) + y_plot(j) <= 1
            U(i,j) = u(i,j);
            V(i,j) = v(i,j);
        else
            U(i,j) = NaN;
            V(i,j) = NaN;

        end
    end
end

for i = 1:NN+1
    for j = 1:NN+1
        if x_vecf(i) + y_vecf(j) > 1
            Zmat(i,j) = NaN;
        end
    end
end
%
figure
contourf(xx,yy,Zmat);
hold on
for ii = 2:20
    for ij = [1:1:N+1]
        if isnan(U(ii,ij)) || isnan(V(ii,ij))
            continue
        else
            ah = annotation('arrow',...
                'Color', [169/255,169/255,169/255],...
                'headStyle','cback1','HeadLength',3.5,'HeadWidth',3.5);
            set(ah,'parent',gca);
            set(ah,'position',[x(ii,ij) y(ii,ij) 0.04*U(ii,ij) 0.04*V(ii,ij)]);
        end
    end
end
hold on


for ii = 1:9
    xc = 0.1*ii;
    yc =  1 - xc;
    arrow_index1 = 6;
    Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);
end
for jj = 1:5
    xc = 0.01*jj;
    yc =  0.01;
    arrow_index1 = 22 -jj;
    Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);
end
xc = 0.009; yc = 0.02;
arrow_index1 = 21;
Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);

xc = 0.015; yc = 0.85;
arrow_index1 = 5;
Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);

xc = 0.05; yc = 0.95;
arrow_index1 = 4;
Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);

xc = 0.004; yc = 0.02;
arrow_index1 = 22 -jj;
arrow_index2 = 24;
Add_trajectory2(xc,yc,arrow_index1,arrow_index2,a,fk,fs,epsilon,gamma,rks);

% Plot the trajectories passing through the saddle point
hyperbolic_traj(pt3,pt4, a,fk,fs,epsilon,gamma,rks,0.0257);
% Plot the fixed points
plot(0,0,'o','MarkerEdgeColor','k','MarkerFaceColor','b','markersize',6)
plot(pt1,0,'o','MarkerEdgeColor','k','MarkerFaceColor','b','markersize',6)
plot(0,pt2,'o','MarkerEdgeColor','k','MarkerFaceColor','b','markersize',6)
plot(pt3,pt4,'o','MarkerEdgeColor','k','MarkerFaceColor','r','markersize',6)

hold off
axis([0 1 0 1])
axis equal
shading flat

text(-0.1,0.5,'$n_S$','FontSize',15,'interpreter','latex');
text(0.45,-0.05,'$n_K$','FontSize',15,'Interpreter','latex');
set(gca,'visible','off');

set(findall(gca,'type','text'),'visible','on');
colormap("copper")


%%  parameters for delta == 0.35
epsilon = 0.2; % The cost of constitutive toxin production
gamma = 1.0; % The rescaled cost of toxin-production rate
rks = 0.85; % The rescaled ratio of the growth rate of the killer to that of the sensitive
d1 = 0.35; % rescaled basal death rate of the killer
d2 = 0.35; % rescaled basal death rate of the sensitive

% Drift part of the population growth dynamics (Eq.[S2.1] in the SI)
fk = @(nk,ns,a) rks.*(1-epsilon.*a).*nk.*(1-nk-ns) -d1*nk;
fs = @(nk,ns,a) ns.*(1-nk-ns)-a.*gamma.*nk.*ns - d2*ns;
% Fixed points of the system
pt1 = 1 - d1/(rks*(1-epsilon));
pt2 = 1 - d2;
pt3 = d1/(gamma*rks*(1-epsilon)) - d2/gamma;
pt4 = (d2 + gamma)/gamma - (d1*(1+gamma))/(gamma*rks*(1-epsilon));

% Generate the direction field on the (K,S)-triangle
a=1; % Fix the toxin production rate to 1
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
dNs=fs(x,y,a); dNk=fk(x,y,a);
% Normalize the drift
dNs(abs(dNs-0) < eps) = 0;
dNk(abs(dNk-0) < eps) = 0;
rr=sqrt(dNk.^2+dNs.^2);
u=dNk./rr; v=dNs./rr;
% Generate the direction field
U = ones(N+1)*NaN;
V = ones(N+1)*NaN;
% Background matrix for visualization
Zmat = ones(NN+1);
Zmat(1,1) = 0; % For visualization purpose, we set the first element to 0

for i = 1:N+1
    for j = 1:N+1
        if x_plot(i) + y_plot(j) <= 1
            U(i,j) = u(i,j);
            V(i,j) = v(i,j);
        else
            U(i,j) = NaN;
            V(i,j) = NaN;

        end
    end
end

for i = 1:NN+1
    for j = 1:NN+1
        if x_vecf(i) + y_vecf(j) > 1
            Zmat(i,j) = NaN;
        end
    end
end
%
figure
contourf(xx,yy,Zmat);
hold on
for ii = 2:20
    for ij = [1:1:N+1]
        if isnan(U(ii,ij)) || isnan(V(ii,ij))
            continue
        else
            ah = annotation('arrow',...
                'Color', [169/255,169/255,169/255],...
                'headStyle','cback1','HeadLength',3.5,'HeadWidth',3.5);
            set(ah,'parent',gca);
            set(ah,'position',[x(ii,ij) y(ii,ij) 0.04*U(ii,ij) 0.04*V(ii,ij)]);
        end
    end
end
hold on
for ii = 1:9
    xc = 0.105*ii;
    yc =  1 - xc;
    arrow_index1 = 8;
    Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);
end
for jj = 3:4
    xc = 0.011*jj;
    yc =  0.01;
    arrow_index1 = 22 -jj;
    Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);
end
xc = 0.15; yc = 0.85;
arrow_index1 = 10;  
Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);

xc = 0.05; yc = 0.95;
arrow_index1 = 5;
Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);

xc =0.026; yc = 0.02;
arrow_index1 = 18;
Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);

xc = 0.009; yc = 0.02;
arrow_index1 = 17;
Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);

xc = 0.07; yc = 0.007;
arrow_index1 = 14;
Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks);

% Plot the trajectories passing through the saddle point
hyperbolic_traj(pt3,pt4, a,fk,fs,epsilon,gamma,rks,0.35);
% Plot the fixed points
plot(0,0,'o','MarkerEdgeColor','k','MarkerFaceColor','b','markersize',6)
plot(pt1,0,'o','MarkerEdgeColor','k','MarkerFaceColor','b','markersize',6)
plot(0,pt2,'o','MarkerEdgeColor','k','MarkerFaceColor','b','markersize',6)
plot(pt3,pt4,'o','MarkerEdgeColor','k','MarkerFaceColor','r','markersize',6)

hold off
axis([0 1 0 1])
axis equal
shading flat

text(-0.1,0.5,'$n_S$','FontSize',15,'interpreter','latex');
text(0.45,-0.05,'$n_K$','FontSize',15,'Interpreter','latex');
set(gca,'visible','off');

set(findall(gca,'type','text'),'visible','on');
colormap("copper")
%% subfunctions
function full_traj = traj(x0,y0,a,fk,fs,epsilon,gamma,rks)
    % This function computes the trajectory of the system starting from (x0, y0)
    % using the ode45 solver.
    % The trajectory is computed for a fixed toxin production rate 'a'.
    % Input:
    %   x0：initial population of the killer
    %   y0：initial population of the sensitive
    %   a: toxin production rate
    %   fk, fs: rate (drift) functions defining the system dynamics
    %   epsilon, gamma, rks: parameters of the system
    % Output:
    %   full_traj: the trajectory of the system

    tspan = [0 500];
    [T,Z] = ode45((@(t,y)[fk(y(1),y(2),a);fs(y(1),y(2),a)]),...
    tspan,[x0 y0]);
    xval=Z(:,1);
    yval=Z(:,2);
    full_traj  = [xval';yval'];
end


function hyperbolic_traj(pt3,pt4,a,fk,fs,epsilon,gamma,rks, parameter_choice)
    % This function computes the stable and unstable manifolds of the saddle point
    % Input:
    %   pt3, pt4: the coordinates of the saddle point
    %   a: toxin production rate
    %   fk, fs: rate (drift) functions defining the system dynamics
    %   epsilon, gamma, rks: parameters of the system
    %   parameter_choice: the choice of the parameter value to plot the manifold

    tspan = [0 -200];
    tspan_foward = [0 2000];
    % Options for ode45 to include the event function
    options = odeset('Events', @myEventFunction);

    % Stable manifold
    [t_stable_forward, y_stable_forward] = ode45((@(t,y)[fk(y(1),y(2),a);fs(y(1),y(2),a)]), tspan, [pt3 + 1e-3,pt4 + 1e-3],options);
    [t_stable_backward, y_stable_backward] = ode45((@(t,y)[fk(y(1),y(2),a);fs(y(1),y(2),a)]), tspan, [pt3 - 1e-3,pt4 - 1e-3]);

    % Unstable manifold
    [t_unstable_forward, y_unstable_forward] = ode45((@(t,y)[fk(y(1),y(2),a);fs(y(1),y(2),a)]), tspan_foward, [pt3 - 1e-3,pt4 + 1e-3]);
    [t_unstable_backward, y_unstable_backward] = ode45((@(t,y)[fk(y(1),y(2),a);fs(y(1),y(2),a)]), tspan_foward, [pt3 + 1e-3,pt4 - 1e-3]);

    hold on
    plot(y_stable_forward(:, 1), y_stable_forward(:, 2), 'm-.', 'Linewidth',1.4);
    plot(y_stable_backward(:, 1), y_stable_backward(:, 2), 'm-.','Linewidth',1.4);

    % Unstable manifold
    plot(y_unstable_forward(:, 1), y_unstable_forward(:, 2), 'c-.', 'Linewidth',1.4);
    plot(y_unstable_backward(:, 1), y_unstable_backward(:, 2), 'c-.', 'Linewidth',1.4)

    % Calculate the direction of the arrows
    if parameter_choice == 0.2
        arrow_index1 = 12; %0.2
        arrow_index2 = 14;
        arrow_index3 = 98;
        arrow_index4 = 80;
    elseif parameter_choice == 0.0257
        arrow_index1 = 8; %0.0257
        arrow_index2 = 14;
        arrow_index3 = 372;
        arrow_index4 = 380;
    elseif parameter_choice == 0.35
        arrow_index1 = 17; %0.35
        arrow_index2 = 16;
        arrow_index3 = 48;
        arrow_index4 = 45;
    end

    if parameter_choice ~= 0.0257
        % Stable manifold (forward in time)
        dx1 = y_stable_forward(arrow_index1,1) - y_stable_forward(arrow_index1 +1,1);
        dy1 = y_stable_forward(arrow_index1,2) - y_stable_forward(arrow_index1 +1,2);

        ah = annotation('arrow',...
            'Color', 'm',...
            'headStyle','vback2','HeadLength',7,'HeadWidth',10,'LineStyle','none');
        set(ah,'parent',gca);
        
        set(ah,'position',[y_stable_forward(arrow_index1,1) y_stable_forward(arrow_index1,2) 0.04*dx1 0.04*dy1]);
    end
    % Stable manifold (backward in time)
    dx2 = y_stable_backward(arrow_index2,1) - y_stable_backward(arrow_index2 +1,1);
    dy2 = y_stable_backward(arrow_index2,2) - y_stable_backward(arrow_index2 +1,2);
    % Unstable manifold (forward in time)
    dx3 = y_unstable_forward(arrow_index3 + 1,1) - y_unstable_forward(arrow_index3,1);
    dy3 = y_unstable_forward(arrow_index3 + 1,2) - y_unstable_forward(arrow_index3,2);
    % Unstable manifold (backward in time)
    dx4 = y_unstable_backward(arrow_index4 + 1,1) - y_unstable_backward(arrow_index4,1);
    dy4 = y_unstable_backward(arrow_index4 + 1,2) - y_unstable_backward(arrow_index4,2);

   
    ah = annotation('arrow',...
        'Color', 'm',...
        'headStyle','vback2','HeadLength',7,'HeadWidth',10,'LineStyle','none');
    set(ah,'parent',gca);

    set(ah,'position',[y_stable_backward(arrow_index2,1) y_stable_backward(arrow_index2,2) 0.04*dx2 0.04*dy2]);
    
    ah = annotation('arrow',...
        'Color', 'c',...
        'headStyle','vback2','HeadLength',7,'HeadWidth',10,'LineStyle','none');
    set(ah,'parent',gca);

    set(ah,'position',[y_unstable_forward(arrow_index3,1) y_unstable_forward(arrow_index3,2) 0.04*dx3 0.04*dy3]);


    ah = annotation('arrow',...
        'Color', 'c',...
        'headStyle','vback2','HeadLength',7,'HeadWidth',10,'LineStyle','none');
    set(ah,'parent',gca);

    set(ah,'position',[y_unstable_backward(arrow_index4,1) y_unstable_backward(arrow_index4,2) 0.04*dx4 0.04*dy4]);


    hold on
end


function addingArrow1(arrow_index1,x,y)
    % This function adds an arrow to the plot at a specified index
    % Input:
    %   arrow_index1: index of the point where the arrow should be added
    %   x, y: coordinates of the points in the trajectory

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

function addingArrow2(arrow_index1,arrow_index2,x,y)
    % This function adds two arrows to the plot at specified indices
    % Input:
    %   arrow_index1: index of the first point where the arrow should be added
    %   arrow_index2: index of the second point where the arrow should be added
    %   x, y: coordinates of the points in the trajectory

    % Calculate the direction of the arrows
    dx1 = x(arrow_index1 + 1) - x(arrow_index1);
    dy1 = y(arrow_index1 + 1) - y(arrow_index1);

    dx2 = x(arrow_index2 + 1) - x(arrow_index2);
    dy2 = y(arrow_index2 + 1) - y(arrow_index2);


    % Plot the arrows using quiver
    ah = annotation('arrow',...
        'Color', 'k',...
        'headStyle','vback2','HeadLength',7,'HeadWidth',11,'LineStyle','none');
    set(ah,'parent',gca);
    set(ah,'position',[x(arrow_index1) y(arrow_index1) 0.04*dx1 0.04*dy1]);
    ah = annotation('arrow',...
        'Color', 'k',...
        'headStyle','vback2','HeadLength',7,'HeadWidth',11,'LineStyle','none');
    set(ah,'parent',gca);
    set(ah,'position',[x(arrow_index2) y(arrow_index2) 0.04*dx2 0.04*dy2]);
end



function [value, isterminal, direction] = myEventFunction(t, y)
    % This function defines an event for the ODE solver to stop integration
    % when the condition x + y = 1 is met.
    % Input:
    %   t: time
    %   y: solution vector
    % Output:
    %   value: value of the event function
    %   isterminal: flag indicating whether the integration should stop
    %   direction: direction of the event function

    % Event function to stop integration when x + y = 1
    value = y(1) + y(2) - 1; % Detect when x + y = 1
    isterminal = 1; % Stop the integration
    direction = 0; % Detect all zero crossings

end

function Add_trajectory(xc,yc,arrow_index1,a,fk,fs,epsilon,gamma,rks)
    % Add a trajectory to the plot with one arrow
    % Input:
    %   xc, yc: initial conditions for the trajectory
    %   arrow_index1: index of the point where the arrow should be added
    %   a: toxin production rate
    %   fk, fs: rate (drift) functions defining the system dynamics
    %   epsilon, gamma, rks: parameters
   
    % Add a trajectory to the plot with one arrow
    full_traj = traj(xc,yc,a,fk,fs,epsilon,gamma,rks);
    x_coord = full_traj(1,:);
    y_coord = full_traj(2,:);
    plot(x_coord,y_coord,'k-','LineWidth',1.5);
    addingArrow1(arrow_index1,x_coord,y_coord);
end

function Add_trajectory2(xc,yc,arrow_index1,arrow_index2,a,fk,fs,epsilon,gamma,rks)
    % This function adds a trajectory to the plot with two arrows
    % Input:
    %   xc, yc: initial conditions for the trajectory
    %   arrow_index1: index of the first point where the arrow should be added
    %   arrow_index2: index of the second point where the arrow should be added
    %   a: toxin production rate
    %   fk, fs: rate (drift) functions defining the system dynamics
    %   epsilon, gamma, rks: parameters

    % Add a trajectory to the plot with two arrows
    full_traj = traj(xc,yc,a,fk,fs,epsilon,gamma,rks);
    x_coord = full_traj(1,:);
    y_coord = full_traj(2,:);
    plot(x_coord,y_coord,'k-','LineWidth',1.5);
    addingArrow2(arrow_index1,arrow_index2,x_11,y_11);
end

end