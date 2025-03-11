%% Read in data from files
clear all; clc;
%close all;
precision1 = 'double';
precision2 = 'uint8';
filename1 = 'Max_freq_valuefn.dat';
filename2 = 'Max_freq_policy.dat';

% Read in parameter data
filename_parameters = 'Max_freq_DomainParameters.dat';


% Paramter Values
uFile_params = fopen(filename_parameters);
u_params = fread(uFile_params, 9, precision1);

Nx = u_params(1) + 1
Ny = u_params(2) + 1
dx = u_params(3)
dy = u_params(4)
epsilon = u_params(5)
rks = u_params(6)
gamma = u_params(7)
lamb = u_params(8)
fclose(uFile_params);


% Value function file
wFile = fopen(filename1);
w = fread(wFile, precision1);
w =  reshape(w,[Ny,Nx]);
fclose(wFile);
% Optimal policy file
aFile = fopen(filename2);
policy = fread(aFile, precision2);
policy = reshape(policy,[Ny,Nx]);
fclose(aFile);

%% Plotting
yy = linspace(0,1,Ny);
xx = linspace(0,1,Nx);

[X,Y] = meshgrid(xx,yy);

uu = zeros(Ny, Nx);
uu(w < X) = 1;


figure;
subplot(1,2,1)
contourf(X,Y,w,'Edgecolor','none');
hold on
contour(X,Y,uu, [1], 'm-.','LineWidth',2);
hold off
% hold on
% xline(defeat_thres,'m:','linewidth',2);
% xline(vic_thres,'m:','linewidth',2);
% hold off
xlabel('Fraction of the Killer (f)','FontSize',14);
ylabel('Total population (N)','FontSize',14);
axis equal
title('Value function','fontsize',16,'Interpreter','latex');
colorbar();

subplot(1,2,2)
contourf(X,Y,policy,[0 1]);
% hold on
% contour(X,Y,policy_slice,[1],'c-.','LineWidth',2);
% hold off
% hold on
% xline(defeat_thres,'m:','linewidth',2);
% xline(vic_thres,'m:','linewidth',2);
% hold off
xlabel('Fraction of the Killer (f)','FontSize',14);
ylabel('Total population (N)','FontSize',14);
axis equal
title('Optimal Toxin-On/Off Region','fontsize',16,'Interpreter','latex');
colorbar();

%%
% w_diff = w - w_slice;
% figure;
% contourf(X,Y,w_diff,'Edgecolor','none');
% % hold on
% % contour(X,Y,uu, [1], 'm-.','LineWidth',2);
% % hold off
% % hold on
% % xline(defeat_thres,'m:','linewidth',2);
% % xline(vic_thres,'m:','linewidth',2);
% % hold off
% xlabel('Fraction of the Killer (f)','FontSize',14);
% ylabel('Total population (N)','FontSize',14);
% axis equal
% title('Difference map','fontsize',16,'Interpreter','latex');
% colorbar();

%%
% velocity function for different r_k and r_s
fx = @(x,y,a)  x.*(1-x).*((1-y).*(rks.*(1-epsilon.*a) - 1) + a.*gamma.*x.*y);

fy = @(x,y,a)  y.*(1-y).*(1 + (rks.*(1-epsilon.*a) - 1).*x) - a.*gamma.*y.^2.*x.*(1-x);

A_myo = policy;
N=20;
NN = 1600;
hh = 1/NN;
q=linspace(0,1,N+1);
p=linspace(0,1,N+1);
qq=linspace(0,1,NN+1);
pp=linspace(0,1,NN+1);
[x,y]=meshgrid(q',p');

A_det_down = A_myo(1:NN/N:end,1:NN/N:end);

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
A_det_down = A_myo(1:NN/N:end,1:NN/N:end);
figure
contourf(X,Y,A_myo, [0, 1],'EdgeColor','none');
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
hold off
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
axis equal
colormap("copper")