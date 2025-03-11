%% Read in data from files
clear all; clc;
%close all;
precision1 = 'double';
precision2 = 'uint8';
filename1 = 'Max_freq_finite_horizon_valuefn.dat';
filename2 = 'Max_freq_finite_horizon_policy.dat';
filename3 = 'Max_freq_finite_horizon_inf_limit_f.dat';
filename4 = 'Max_freq_finite_horizon_inf_limit_N.dat';
% Read in parameter data
filename_parameters = 'Max_freq_finite_horizon_DomainParameters.dat';


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
% Optimal policy file
aFile = fopen(filename2);
policy = fread(aFile, precision2);
policy = reshape(policy,[Ny,Nx, fM]);
fclose(aFile);
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
ww = zeros(Ny,Nx);

slice_num = 161;

w_slice1 = w(:,:,slice_num);
policy_slice = policy(:,:,slice_num);
ww(w_slice1 < X) = 1;

f_slice1 = f(:,:,slice_num);
N_slice1 = N(:,:,slice_num);
%%
figure;
% subplot(1,2,1)
contourf(X,Y,w_slice1 - w_slice,'Edgecolor','none');
% hold on
% contour(X,Y,uu, [1], 'r--','LineWidth',2);
% hold off
xlabel('Initial fraction of the killer (x = f_0)','FontSize',14);
ylabel('Initial total population (y = N_0)','FontSize',14);
axis equal
% title('Value function','fontsize',16,'Interpreter','latex');
hcb = colorbar();
% clim([0 1]);
% cc = linspace(0,1,10);
% set(hcb,'Ticks',cc);
% myCticks = arrayfun(@(x) sprintf('%.1f',x),cc,'un',0);
% set(hcb,'Ticklabels',myCticks);
%%
% figure
% %subplot(1,2,2)
% contourf(X,Y,policy_slice,[0 1],'EdgeColor','none');
% xlabel('Fraction of the killer (f)','FontSize',14);
% ylabel('Total population (N)','FontSize',14);
% axis equal
% title('Optimal Toxin-On/Off Region','fontsize',16,'Interpreter','latex');
% colorbar();
f_slice2 = f_slice1;
f_slice2(f_slice2 > 0.999999999) = 1;
figure;
% subplot(1,2,1)
contourf(X,Y,f_slice2,'Edgecolor','none');
% hold on
% contour(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),ff(2:end-1,2:end-1), [1], 'r--','LineWidth',2);
% hold off
xlabel('Initial fraction of the killer (f_0)','FontSize',14);
ylabel('Initial total population (N_0)','FontSize',14);
axis equal
% title('Relative abundance of killer in the limit','fontsize',16,'Interpreter','latex');
% hcb = colorbar();
% clim([0 1]);
% cc = linspace(0,1,10);
% set(hcb,'Ticks',cc);
% myCticks = arrayfun(@(x) sprintf('%.1f',x),cc,'un',0);
% set(hcb,'Ticklabels',myCticks);
%%
figure
% subplot(1,2,2)
% contourf(X,Y,(f_slice1).*N_slice1,12,'Edgecolor','none');
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
% velocity function for different r_k and r_s
fx = @(x,y,a)  x.*(1-x).*((1-y).*(rks.*(1-epsilon.*a) - 1) + a.*gamma.*x.*y);

fy = @(x,y,a)  y.*(1-y).*(1 + (rks.*(1-epsilon.*a) - 1).*x) - a.*gamma.*y.^2.*x.*(1-x);

A_det = policy_slice;
% A_det(:,vic_indx+1:end) = NaN;
% A_det(:,1:defeat_indx) = NaN;
N=20;
NN = 1600;
hh = 1/NN;
q=linspace(0,1,N+1);
p=linspace(0,1,N+1);
qq=linspace(0,1,NN+1);
pp=linspace(0,1,NN+1);
[x,y]=meshgrid(q',p');

A_det_down = A_det(1:NN/N:end,1:NN/N:end);

dp=fy(x,y,0); dq=fx(x,y,0);

dp(abs(dp-0) < eps) = 0;
dq(abs(dq-0) < eps) = 0;
rr=sqrt(dp.^2+dq.^2);

u=dq./rr; v=dp./rr;


dpp=fy(x,y,1); dqq=fx(x,y,1);

dpp(abs(dp-0) < eps) = 0;
dqq(abs(dq-0) < eps) = 0;
rrr=sqrt(dpp.^2+dqq.^2);

uu=dqq./rrr; vv=dpp./rrr;
A_det_down = A_det(1:NN/N:end,1:NN/N:end);
figure
contourf(X,Y,A_det, [0, 1],'EdgeColor','none');
hold on
for ii = 2:N
    for ij = [2:1:N+1]
        if isnan(u(ii,ij)) || isnan(v(ii,ij))
            continue
        else
            ah = annotation('arrow',...
                'Color', [169/255,169/255,169/255],...
                'headStyle','cback1','HeadLength',3.5,'HeadWidth',3.5);
            set(ah,'parent',gca);
            if A_det_down(ii,ij) == 0
                set(ah,'position',[x(ii,ij) y(ii,ij) 0.03*u(ii,ij) 0.03*v(ii,ij)]);
            else
                set(ah,'position',[x(ii,ij) y(ii,ij) 0.03*uu(ii,ij) 0.03*vv(ii,ij)]);
            end
        end
    end
end
% xline(vic_thres,'m:','linewidth',2);
% xline(defeat_thres,'m:','linewidth',2);
% scatter(0.558,0.289,40,'r.');
hold off
xlabel('Initial fraction of the killer (x = f_0)','FontSize',14);
ylabel('Initial total population (y = N_0)','FontSize',14);
axis equal
colormap("copper")

%% jump discontinuity
q=linspace(0,1,NN+1);
p=linspace(0,1,NN+1);
[x,y]=meshgrid(q',p');

dp=fy(x,y,0); dq=fx(x,y,0);

dp(abs(dp-0) < eps) = 0;
dq(abs(dq-0) < eps) = 0;
rr=sqrt(dp.^2+dq.^2);

u=dq./rr; v=dp./rr;


dpp=fy(x,y,1); dqq=fx(x,y,1);

dpp(abs(dp-0) < eps) = 0;
dqq(abs(dq-0) < eps) = 0;
rrr=sqrt(dpp.^2+dqq.^2);

uu=dqq./rrr; vv=dpp./rrr;

figure
subplot(1,2,1)
contourf(x,y,dq);
hold on
contour(x(:,2:1:end-1),y(:,2:1:end-1),policy_slice(:,2:1:end-1),[1],'w--','LineWidth',2);
hold off
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
colorbar();
% clim([-0.08,0.13]);
axis square
title('a = 0','Fontsize',18);

subplot(1,2,2)
contourf(x,y,dqq);
hold on
contour(x(:,2:1:end-1),y(:,2:1:end-1),policy_slice(:,2:1:end-1),[1],'w--','LineWidth',2);
hold off
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
colorbar();
% clim([-0.08,0.13]);
axis square
title('a = 1','Fontsize',18);

sgtitle('Comparison of df','Fontsize',20);

figure
subplot(1,2,1)
contourf(x,y,dp);
hold on
contour(x(:,2:1:end-1),y(:,2:1:end-1),policy_slice(:,2:1:end-1),[1],'w--','LineWidth',2);

hold off
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
axis square
colorbar();
% clim([-0.25,0.25]);
title('a = 0','Fontsize',18);

subplot(1,2,2)
contourf(x,y,dpp);
hold on
contour(x(:,2:1:end-1),y(:,2:1:end-1),policy_slice(:,2:1:end-1),[1],'w--','LineWidth',2);

hold off
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
axis square
colorbar();
% clim([-0.25,0.25]);
title('a = 1','Fontsize',18);

sgtitle('Comparison of dN','Fontsize',20);

df_opt = dp;
dN_opt = dq;

for i = 1:Nx
    for j = 1:Ny
        if policy_slice(j,i) ~= 0
            df_opt(j,i) = dpp(j,i);
            dN_opt(j,i) = dqq(j,i);
        end
    end
end

figure
subplot(1,3,1)
contourf(x,y,df_opt);
hold on
contour(x(:,2:1:end-1),y(:,2:1:end-1),policy_slice(:,2:1:end-1),[1],'w--','LineWidth',2);

hold off
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
colorbar();
clim([-0.08,0.13]);
axis square
title('df','Fontsize',18);

subplot(1,3,2)
contourf(x,y,dN_opt);
hold on
contour(x(:,2:1:end-1),y(:,2:1:end-1),policy_slice(:,2:1:end-1),[1],'w--','LineWidth',2);

hold off
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
colorbar();
% clim([-0.25,0.25]);
axis square
title('dN','Fontsize',18);

axx = subplot(1,3,3);
contourf(x,y,policy_slice);
% hold on
% xline(defeat_thres,'m:','linewidth',2);
% xline(vic_thres,'m:','linewidth',2);
% % contour(x(:,1:1:1584),y(:,1:1:1584),A_ult(:,1:1:1584),[1],'r--','LineWidth',2);
% hold off
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
colorbar();
axis square
colormap(axx,"copper");
title('Optimal policy','Fontsize',18);

sgtitle('Contour plots of df and dN under the optimal policy','Fontsize',20);

%% jump discontinuity cont.
norm_new = zeros(1,Nx);
norm_og = zeros(1,Nx);
for i = 2:Nx-2
    for j = 1:Ny-1
        if policy_slice(j,i) == 0 && policy_slice(j+1,i) == 1
            df_og = dq(j,i);
            dN_og = dp(j,i);
            df_new = dqq(j+1,i);
            dN_new = dpp(j+1,i);
            break
        end

    end
    norm_og(i) = norm([df_og,dN_og]);
    norm_new(i) = norm([df_new,dN_new]);
    df_og_norm(i) = df_og./norm_og(i);
    dN_og_norm(i)= dN_og./norm_og(i);
    df_new_norm(i) = df_new./norm_new(i);
    dN_new_norm(i) = dN_new./norm_new(i);
end



cos_theta = zeros(1,Nx);

for i = 2:Nx-2
    cos_theta(i) = dot([df_og_norm(i),dN_og_norm(i)],[df_new_norm(i),dN_new_norm(i)]);
end
theta = acos(cos_theta);
theta_trunc = theta(2:end-2);
%per_jump = (theta./pi)*100;
max(theta_trunc)

%%
% Extract the contour lines
C = contourc(X(1,:), Y(:,1), ww, [0.5 0.5]); % Use 0.5 as the level to extract the boundary

% Extract the x and y coordinates of the boundary
x_boundary = C(1, 2:end-1);
y_boundary = C(2, 2:end-1);

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
u11 = w_slice1(ky-1,kx-1);
u12 = w_slice1(ky-1,kx);
u21 = w_slice1(ky,kx-1);
u22 = w_slice1(ky,kx);
%bilinear interpolation
w_interp = BilinearInterp(u11,u12,u21,u22,xleft,xright,...
    ytop,ybot,x0,y0,dx,dy);


%%
TT = 1;
AA = zeros(Ny,Nx,161);
for i = 1:161
    AA(:,:,i) = policy(:,:,161-i+1);
end
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