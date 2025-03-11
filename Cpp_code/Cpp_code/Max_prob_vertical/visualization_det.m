%% Read in data from files
clear all; clc;
%close all;
precision1 = 'double';
precision2 = 'uint8';
filename1 = 'output/Max_prob_vertical_valuefn_rho650_lamb1000.dat';
filename2 = 'output/Max_prob_vertical_policy_rho650_lamb1000.dat';

% Read in parameter data
filename_parameters = 'Max_prob_vertical_DomainParameters.dat';


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

lamb = 1
rho = 0.65
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

filename2 = ['/space/mw929/Documents/CPP/KS-competition/Performance_given_policy/output/Perf_const_valuefn_rho',num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat'];

% Value function file
wFile = fopen(filename2);
w_const = fread(wFile, precision1);
w_const=  reshape(w_const,[1601,1601]);
fclose(wFile);

filename2 = ['/space/mw929/Documents/CPP/KS-competition/Performance_given_policy/output/Perf_myo_valuefn_rho',num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat'];

% Value function file
wFile = fopen(filename2);
w_myo = fread(wFile, precision1);
w_myo=  reshape(w_myo,[1601,1601]);
fclose(wFile);

load("a_myo.mat")
%% Plotting
yy = linspace(0,1,Ny);
xx = linspace(0,1,Nx);

[X,Y] = meshgrid(xx,yy);


figure;
% subplot(1,2,1)
contourf(X,Y,w,12,'Edgecolor','none');
hold on
xline(defeat_thres,'m:','linewidth',2);
xline(vic_thres,'m:','linewidth',2);
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


% subplot(1,2,2)
% contourf(X,Y,policy,[0 1]);
% hold on
% xline(defeat_thres,'m:','linewidth',2);
% xline(vic_thres,'m:','linewidth',2);
% hold off
% xlabel('Frequency (f)','FontSize',14);
% ylabel('Total population (N)','FontSize',14);
% axis equal
% title('Optimal Toxin-On/Off Region','fontsize',16,'Interpreter','latex');
% colorbar();


%%
% velocity function for different r_k and r_s
fx = @(x,y,a)  x.*(1-x).*((1-y).*(rks.*(1-epsilon.*a) - 1) + a.*gamma.*x.*y);

fy = @(x,y,a)  y.*(1-y).*(1 + (rks.*(1-epsilon.*a) - 1).*x) - a.*gamma.*y.^2.*x.*(1-x);

A_ult = policy;
A_ult(:,vic_indx+1:end) = 0;
A_ult(:,1:defeat_indx) = 0;
N=20;
NN = 1600;
hh = 1/NN;
q=linspace(0,1,N+1);
p=linspace(0,1,N+1);
qq=linspace(0,1,NN+1);
pp=linspace(0,1,NN+1);
[x,y]=meshgrid(q',p');

A_det_down = A_ult(1:NN/N:end,1:NN/N:end);

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
A_det_down = A_ult(1:NN/N:end,1:NN/N:end);
figure
contourf(X,Y,A_ult, [0, 1],'EdgeColor','none');
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
contour(X(:,1:1:1584),Y(:,1:1:1584),a_myo(:,1:1:1584),[1],'w--','LineWidth',2.5);
xline(vic_thres,'m:','linewidth',2);
xline(defeat_thres,'m:','linewidth',2);
hold off
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
axis equal
colormap("copper")

%%
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
xline(defeat_thres,'m:','linewidth',2);
xline(vic_thres,'m:','linewidth',2);
contour(x(:,1:1:1584),y(:,1:1:1584),A_ult(:,1:1:1584),[1],'w--','LineWidth',2);
hold off
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
colorbar();
clim([-0.08,0.13]);
axis square
title('a = 0','Fontsize',18);

subplot(1,2,2)
contourf(x,y,dqq);
hold on
xline(defeat_thres,'m:','linewidth',2);
xline(vic_thres,'m:','linewidth',2);
contour(x(:,1:1:1584),y(:,1:1:1584),A_ult(:,1:1:1584),[1],'w--','LineWidth',2);
hold off
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
colorbar();
clim([-0.08,0.13]);
axis square
title('a = 1','Fontsize',18);

sgtitle('Comparison of df','Fontsize',20);

figure
subplot(1,2,1)
contourf(x,y,dp);
hold on
xline(defeat_thres,'m:','linewidth',2);
xline(vic_thres,'m:','linewidth',2);
contour(x(:,1:1:1584),y(:,1:1:1584),A_ult(:,1:1:1584),[1],'w--','LineWidth',2);
hold off
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
axis square
colorbar();
clim([-0.25,0.25]);
title('a = 0','Fontsize',18);

subplot(1,2,2)
contourf(x,y,dpp);
hold on
xline(defeat_thres,'m:','linewidth',2);
xline(vic_thres,'m:','linewidth',2);
contour(x(:,1:1:1584),y(:,1:1:1584),A_ult(:,1:1:1584),[1],'w--','LineWidth',2);
hold off
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
axis square
colorbar();
clim([-0.25,0.25]);
title('a = 1','Fontsize',18);

sgtitle('Comparison of dN','Fontsize',20);

df_opt = dq;
dN_opt = dp;

for i = 1:Nx
    for j = 1:Ny
        if A_ult(j,i) ~= 0
            df_opt(j,i) = dqq(j,i);
            dN_opt(j,i) = dpp(j,i);
        end
    end
end

figure
subplot(1,3,1)
contourf(x,y,df_opt);
hold on
xline(defeat_thres,'m:','linewidth',2);
xline(vic_thres,'m:','linewidth',2);
contour(x(:,1:1:1584),y(:,1:1:1584),A_ult(:,1:1:1584),[1],'w--','LineWidth',2);
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
xline(defeat_thres,'m:','linewidth',2);
xline(vic_thres,'m:','linewidth',2);
contour(x(:,1:1:1584),y(:,1:1:1584),A_ult(:,1:1:1584),[1],'w--','LineWidth',2);
hold off
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
colorbar();
clim([-0.25,0.25]);
axis square
title('dN','Fontsize',18);

axx = subplot(1,3,3);
contourf(x,y,A_ult);
hold on
xline(defeat_thres,'m:','linewidth',2);
xline(vic_thres,'m:','linewidth',2);
% contour(x(:,1:1:1584),y(:,1:1:1584),A_ult(:,1:1:1584),[1],'r--','LineWidth',2);
hold off
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
colorbar();
axis square
colormap(axx,"copper");
title('Optimal policy','Fontsize',18);

sgtitle('Contour plots of df and dN under the optimal policy','Fontsize',20);
%%
% velocity function for different r_k and r_s
fx = @(x,y,a)  x.*(1-x).*((1-y).*(rks.*(1-epsilon.*a) - 1) + a.*gamma.*x.*y);

fy = @(x,y,a)  y.*(1-y).*(1 + (rks.*(1-epsilon.*a) - 1).*x) - a.*gamma.*y.^2.*x.*(1-x);

A_myo = a_myo;

N=20;
NN = 1600;
hh = 1/NN;
q=linspace(0,1,N+1);
p=linspace(0,1,N+1);
qq=linspace(0,1,NN+1);
pp=linspace(0,1,NN+1);
[x,y]=meshgrid(q',p');

A_myo_down = A_myo(1:NN/N:end,1:NN/N:end);

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
A_myo_down = A_myo(1:NN/N:end,1:NN/N:end);
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
            if A_myo_down(ii,ij) == 0
                set(ah,'position',[x(ii,ij) y(ii,ij) 0.03*u(ii,ij) 0.03*v(ii,ij)]);
            else
                set(ah,'position',[x(ii,ij) y(ii,ij) 0.03*uu(ii,ij) 0.03*vv(ii,ij)]);
            end
        end
    end
end
contour(X(:,1:1:1584),Y(:,1:1:1584),A_ult(:,1:1:1584),[1],'r--','LineWidth',2.5);
hold off
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
axis equal
colormap("copper")

%%
rho_list = [0.50,0.55,0.60,0.65,0.7];
lamb_list = [0.75, 0.8,0.85, 0.9, 0.95, 1.0];

ult_table = zeros(length(lamb_list),length(rho_list));
myo_table = zeros(length(lamb_list),length(rho_list));
const_table = zeros(length(lamb_list),length(rho_list));

precision1 = 'double';
x0 = 0.5;
y0 = 0.1;
xx = linspace(0,1,1601);
xindx = find(xx == 0.5,1);
yindx = find(xx == 0.1,1);


for ii = 1:length(lamb_list)

    lamb = lamb_list(ii)



    for jj = 1:length(rho_list)

        rho = rho_list(jj)

        filename2 = ['output/Max_prob_vertical_valuefn_rho',num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat'];

        % Value function file
        wFile = fopen(filename2);
        w_ult = fread(wFile, precision1);
        w_ult=  reshape(w_ult,[1601,1601]);
        fclose(wFile);

        filename2 = ['/space/mw929/Documents/CPP/KS-competition/Performance_given_policy/output/Perf_const_valuefn_rho',num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat'];

        % Value function file
        wFile = fopen(filename2);
        w_const = fread(wFile, precision1);
        w_const=  reshape(w_const,[1601,1601]);
        fclose(wFile);

        filename2 = ['/space/mw929/Documents/CPP/KS-competition/Performance_given_policy/output/Perf_myo_valuefn_rho',num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat'];

        % Value function file
        wFile = fopen(filename2);
        w_myo = fread(wFile, precision1);
        w_myo=  reshape(w_myo,[1601,1601]);
        fclose(wFile);

        ult_table(ii, jj) = w_ult(yindx, xindx)

        myo_table(ii, jj) = w_myo(yindx, xindx)

        const_table(ii, jj) = w_const(yindx, xindx)
    end
end

%% ult

figure
imagesc(rho_list,lamb_list,ult_table);
view([0 0 1]);
xlabel('Survival fraction ($\rho$)','FontSize',20,'Interpreter','latex');
ylabel('Dilution frequency ($\lambda$)','FontSize',20,'Interpreter','latex');
colorbar();
yticks(lamb_list);
clim([0 1]);
colormap("hot")
axis square
set(gca,'YDir','normal');

%% myo

[Xmap,Ymap] = meshgrid(rho_list,lamb_list);
figure
imagesc(rho_list,lamb_list,myo_table);
view([0 0 1]);
xlabel('Survival fraction ($\rho$)','FontSize',20,'Interpreter','latex');
ylabel('Dilution frequency ($\lambda$)','FontSize',20,'Interpreter','latex');
colorbar();
yticks(lamb_list);
clim([0 1]);
colormap("hot")
axis square
set(gca,'YDir','normal');

%% const

figure
imagesc(rho_list,lamb_list,const_table);
view([0 0 1]);
xlabel('Survival fraction ($\rho$)','FontSize',20,'Interpreter','latex');
ylabel('Dilution frequency ($\lambda$)','FontSize',20,'Interpreter','latex');
colorbar();
yticks(lamb_list);
clim([0 1]);
colormap("hot")
axis square
set(gca,'YDir','normal');

%%
[Xmap,Ymap] = meshgrid(rho_list,lamb_list);
figure
imagesc(rho_list,lamb_list,ult_table - myo_table);
view([0 0 1]);
xlabel('Survival fraction ($\rho$)','FontSize',20,'Interpreter','latex');
ylabel('Dilution frequency ($\lambda$)','FontSize',20,'Interpreter','latex');
colorbar();
yticks(lamb_list);
% clim([0 1]);
colormap("cool")
axis square
set(gca,'YDir','normal');

%%
W1 = (w - w_const);
W3 = (w - w_myo);

figure
[AA,BB] = contourf(X,Y,W1,7,"ShowText",true,"LabelFormat","%.2f");
hold on
xline(vic_thres,'m:','linewidth',2);
xline(defeat_thres,'m:','linewidth',2);
hold off
% set(BB,'LineColor','none');
xlabel('Initial fraction of the killer (x = f_0)','FontSize',14);
ylabel('Initial total population (y = N_0)','FontSize',14);
axis equal
colorbar();

figure
[AA,BB] = contourf(X,Y,W3,5,"ShowText",true,"LabelFormat","%.3f");
hold on
xline(vic_thres,'m:','linewidth',2);
xline(defeat_thres,'m:','linewidth',2);
hold off
% set(BB,'LineColor','none');
xlabel('Initial fraction of the killer (x = f_0)','FontSize',14);
ylabel('Initial total population (y = N_0)','FontSize',14);
axis equal
colorbar();

