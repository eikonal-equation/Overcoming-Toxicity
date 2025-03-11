%% Read in data from files
clear all; clc;
%close all;
precision1 = 'double';
precision2 = 'uint8';

% Read in parameter data
filename_parameters = 'Max_prob_hyperbolic_DomainParameters.dat';


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


lamb = 1.0;
rho = 0.65;
filename1 = ['output/Max_prob_hyper_valuefn_rho',num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat'];
filename2 = ['output/Max_prob_hyper_policy_rho',num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat'];
filename3 = ['/space/mw929/Documents/CPP/KS-competition/Max_prob_vertical/output/Max_prob_vertical_policy_rho',num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat'];
filename4 = ['/space/mw929/Documents/CPP/KS-competition/Max_prob_vertical/output/Max_prob_vertical_valuefn_rho',num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat'];

% Value function file
wFile = fopen(filename1);
w_ult = fread(wFile, precision1);
w_ult=  reshape(w_ult,[1601,1601]);
fclose(wFile);
% Value function file
wFile = fopen(filename4);
w_ver = fread(wFile, precision1);
w_ver=  reshape(w_ver,[1601,1601]);
fclose(wFile);
% Value function file
wFile = fopen(filename2);
a_ult = fread(wFile, precision2);
a_ult=  reshape(a_ult,[1601,1601]);
fclose(wFile);
% Value function file
wFile = fopen(filename3);
a_ver = fread(wFile, precision2);
a_ver=  reshape(a_ver,[1601,1601]);
fclose(wFile);
%% Plotting
yy = linspace(0,1,Ny);
xx = linspace(0,1,Nx);

[X,Y] = meshgrid(xx,yy);

ff = @(x) 1./x;

figure;
% subplot(1,2,1)
contourf(X,Y,w_ult,12,'Edgecolor','none');
hold on
plot(xx, vic_thres.*ff(xx),'m:','linewidth',2)
plot(xx, defeat_thres.*ff(xx),'m:','linewidth',2)
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
xlim([0,1]);
ylim([0,1]);

%%
W_diff = abs(w_ult - w_ver);
mean_row = mean(W_diff');

% for i = 1:1601
%     for j = 1:1601
%         xi = (j-1)*dx;
%         yj = (i-1)*dx;
%         if xi*yj < defeat_thres || xi*yj > vic_thres
%             W_diff(i,j) = NaN;
%         end
%     end
% end
% figure;
% % subplot(1,2,1)
% [AA,BB] = contourf(X,Y,W_diff,5,"ShowText",true,"LabelFormat","%.2f");
% xlabel('Initial fraction of the killer (x)','FontSize',14);
% ylabel('Initial total population (y)','FontSize',14);
% axis equal
% % title('Value function','fontsize',16,'Interpreter','latex');
% hcb = colorbar();

figure
plot(xx(1:end),mean_row(1:end),'.-','Linewidth',2);
xlabel('Initial total population (y = N_0)','FontSize',14);
ylabel('Mean absolute difference','FontSize',14);
axis normal
grid minor
xlim([-0.05,1.05]);
ylim([0,0.26]);
%%
% velocity function for different r_k and r_s
fx = @(x,y,a)  x.*(1-x).*((1-y).*(rks.*(1-epsilon.*a) - 1) + a.*gamma.*x.*y);

fy = @(x,y,a)  y.*(1-y).*(1 + (rks.*(1-epsilon.*a) - 1).*x) - a.*gamma.*y.^2.*x.*(1-x);

A_ult = a_ult;
for i = 1:1601
    for j = 1:1601
        xi = (j-1)*dx;
        yj = (i-1)*dx;
        if xi*yj < defeat_thres || xi*yj > vic_thres
            A_ult(i,j) = 0;
        end
    end
end
% A_ult(:,vic_indx+1:end) = NaN;
% A_ult(:,1:defeat_indx) = NaN;
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
ij1 = [2:1:N+1];
ij2 = [3:1:N+1];
ij3 = [2:1:N];
for ii = 2:N
    if ii == 2
        ijorder = ij2;
    elseif ii == N
        ijorder = ij3;
    else
        ijorder = ij1;
    end
   for ij = ijorder
        

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
% contourf(X,Y,a_ult,[0 1],'Edgecolor','none');
hold on
contour(X(:,1:1:1584),Y(:,1:1:1584),a_ver(:,1:1:1584),[1],'w--','LineWidth',2.5);
plot(xx, vic_thres.*ff(xx),'m:','linewidth',2)
plot(xx, defeat_thres.*ff(xx),'m:','linewidth',2)
hold off
xlabel('Fraction of the killer (f)','FontSize',14);
ylabel('Total population (N)','FontSize',14);
axis equal
xlim([0,1]);
ylim([0,1]);
colormap("copper");


%%
%close all;
precision1 = 'double';
precision2 = 'uint8';
filename1 = 'Max_prob_hyperbolic_valuefn_100.dat';
filename2 = 'Max_prob_hyperbolic_policy_100.dat';

% Read in parameter data
filename_parameters = 'Max_prob_hyperbolic_DomainParameters_100.dat';


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
lamb = u_params(8)
rho = u_params(9)
vic_thres = u_params(10)
defeat_thres = u_params(11)
fclose(uFile_params);

defeat_indx = floor(defeat_thres/dx)+1;
vic_indx = floor(vic_thres/dx)+1;

% Value function file
wFile = fopen(filename1);
w_100 = fread(wFile, precision1);
w_100 =  reshape(w_100,[Ny,Nx]);
fclose(wFile);

%%
filename1 = 'Max_prob_hyperbolic_valuefn_200.dat';
filename2 = 'Max_prob_hyperbolic_policy_200.dat';

% Read in parameter data
filename_parameters = 'Max_prob_hyperbolic_DomainParameters_200.dat';


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
lamb = u_params(8)
rho = u_params(9)
vic_thres = u_params(10)
defeat_thres = u_params(11)
fclose(uFile_params);

defeat_indx = floor(defeat_thres/dx)+1;
vic_indx = floor(vic_thres/dx)+1;

% Value function file
wFile = fopen(filename1);
w_200 = fread(wFile, precision1);
w_200 =  reshape(w_200,[Ny,Nx]);
fclose(wFile);

%%
filename1 = 'Max_prob_hyperbolic_valuefn_400.dat';
filename2 = 'Max_prob_hyperbolic_policy_400.dat';

% Read in parameter data
filename_parameters = 'Max_prob_hyperbolic_DomainParameters_400.dat';


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
lamb = u_params(8)
rho = u_params(9)
vic_thres = u_params(10)
defeat_thres = u_params(11)
fclose(uFile_params);

defeat_indx = floor(defeat_thres/dx)+1;
vic_indx = floor(vic_thres/dx)+1;

% Value function file
wFile = fopen(filename1);
w_400 = fread(wFile, precision1);
w_400 =  reshape(w_400,[Ny,Nx]);
fclose(wFile);


%%
filename1 = 'Max_prob_hyperbolic_valuefn_800.dat';
filename2 = 'Max_prob_hyperbolic_policy_800.dat';

% Read in parameter data
filename_parameters = 'Max_prob_hyperbolic_DomainParameters_800.dat';


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
lamb = u_params(8)
rho = u_params(9)
vic_thres = u_params(10)
defeat_thres = u_params(11)
fclose(uFile_params);

defeat_indx = floor(defeat_thres/dx)+1;
vic_indx = floor(vic_thres/dx)+1;

% Value function file
wFile = fopen(filename1);
w_800 = fread(wFile, precision1);
w_800 =  reshape(w_800 ,[Ny,Nx]);
fclose(wFile);


%%
precision1 = 'double'
filename1 = 'Max_prob_hyperbolic_valuefn_1600.dat';
filename2 = 'Max_prob_hyperbolic_policy_1600.dat';

% Read in parameter data
filename_parameters = 'Max_prob_hyperbolic_DomainParameters_1600.dat';


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
lamb = u_params(8)
rho = u_params(9)
vic_thres = u_params(10)
defeat_thres = u_params(11)
fclose(uFile_params);

defeat_indx = floor(defeat_thres/dx)+1;
vic_indx = floor(vic_thres/dx)+1;

% Value function file
wFile = fopen(filename1);
w_1600 = fread(wFile, precision1);
w_1600 =  reshape(w_1600,[Ny,Nx]);
fclose(wFile);


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
lamb = u_params(8)
rho = u_params(9)
vic_thres = u_params(10)
defeat_thres = u_params(11)
fclose(uFile_params);

defeat_indx = floor(defeat_thres/dx)+1;
vic_indx = floor(vic_thres/dx)+1;

%% Linf err
l_inf = zeros(4,1);
l_inf(1) = max(max(abs(w_1600(1:16:end,1:16:end)-w_100)));
l_inf(2) = max(max(abs(w_1600(1:16:end,1:16:end)-w_200(1:2:end,1:2:end))));
l_inf(3) = max(max(abs(w_1600(1:16:end,1:16:end)-w_400(1:4:end,1:4:end))));
l_inf(4) = max(max(abs(w_1600(1:16:end,1:16:end)-w_800(1:8:end,1:8:end))));
Delta_x = [1/100,1/200,1/400,1/800];
figure
loglog(Delta_x,l_inf,'ro-','markersize',5.5,'MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',1.5);
hold on
loglog(Delta_x,Delta_x*50,'k--','LineWidth',1.2);
hold off
xlabel('$\Delta x$','FontSize',14,'Interpreter','latex');
ylabel('$l_{\infty}$-Error','FontSize',14,'Interpreter','latex');
legend('$l_{\infty}$ err','linear rate','fontsize',10,'Interpreter','latex',...
    'Fontsize',16,'location','northwest');
grid minor
axis square

ax = gca;
% Define and set xticks
% xticks = [10, 100, 500, 1000]; % Example xtick locations
xticks = [1/800, 1/400, 1/200, 1/100];
ax.XTick = xticks;

% Format x-axis tick labels
newXTickLabels = arrayfun(@(v) sprintf('%.3g', v), xticks, 'UniformOutput', false);
ax.XTickLabel = newXTickLabels;


%% L1 err
l1_err = zeros(4,1);
l1_err(1) = sum(sum(sum(abs(w_1600(1:16:end,1:16:end)-w_100))))/(101^2);
l1_err(2) = sum(sum(sum(abs(w_1600(1:16:end,1:16:end)-w_200(1:2:end,1:2:end)))))/(101^2);
l1_err(3) = sum(sum(sum(abs(w_1600(1:16:end,1:16:end)-w_400(1:4:end,1:4:end)))))/(101^2);
l1_err(4) = sum(sum(sum(abs(w_1600(1:16:end,1:16:end)-w_800(1:8:end,1:8:end)))))/(101^2);

figure
loglog(Delta_x,l1_err,'ro-','markersize',5.5,'MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',1.5);
hold on
loglog(Delta_x,Delta_x,'k--','LineWidth',1.2);
hold off
xlabel('$\Delta x$','FontSize',14,'Interpreter','latex');
ylabel('$l_1$-Error','FontSize',14,'Interpreter','latex');
legend('$l_1$ err','linear rate','fontsize',10,'Interpreter','latex',...
    'Fontsize',16,'location','northwest');
grid minor
axis square

ax = gca;


% Define and set xticks
% xticks = [10, 100, 500, 1000]; % Example xtick locations
xticks = [1/800, 1/400, 1/200, 1/100];
ax.XTick = xticks;

% Format x-axis tick labels
newXTickLabels = arrayfun(@(v) sprintf('%.3g', v), xticks, 'UniformOutput', false);
ax.XTickLabel = newXTickLabels;

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

        filename2 = ['output/Max_prob_hyper_valuefn_rho',num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat'];

        % Value function file
        wFile = fopen(filename2);
        w_ult = fread(wFile, precision1);
        w_ult=  reshape(w_ult,[1601,1601]);
        fclose(wFile);

        filename2 = ['/space/mw929/Documents/CPP/KS-competition/Performance_given_policy_hyper/output/Perf_hyper_const_valuefn_rho',num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat'];

        % Value function file
        wFile = fopen(filename2);
        w_const = fread(wFile, precision1);
        w_const=  reshape(w_const,[1601,1601]);
        fclose(wFile);

        filename2 = ['/space/mw929/Documents/CPP/KS-competition/Performance_given_policy_hyper/output/Perf_hyper_myo_valuefn_rho',num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat'];

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
% view([0 0 1]);
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


%% expected value table

rho_list = [0.50,0.55,0.60,0.65,0.7];
lamb_list = [0.75, 0.8,0.85, 0.9, 0.95, 1.0];
% rho_list = [0.5];
% lamb_list = [0.85];
mean_table_ult = zeros(length(lamb_list),length(rho_list));
mean_table_const = zeros(length(lamb_list),length(rho_list));
mean_table_myo = zeros(length(lamb_list),length(rho_list));

precision1 = 'double';
for ii = 1:length(lamb_list)

    lamb = lamb_list(ii)



    for jj = 1:length(rho_list)

        rho = rho_list(jj)
        filename2 = ['output/Max_prob_hyper_valuefn_rho',num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat'];

        % Value function file
        wFile = fopen(filename2);
        w_ult = fread(wFile, precision1);
        w_ult=  reshape(w_ult,[1601,1601]);
        fclose(wFile);

        filename2 = ['/space/mw929/Documents/CPP/KS-competition/Performance_given_policy_hyper/output/Perf_hyper_const_valuefn_rho',num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat'];

        % Value function file
        wFile = fopen(filename2);
        w_const = fread(wFile, precision1);
        w_const=  reshape(w_const,[1601,1601]);
        fclose(wFile);

        filename2 = ['/space/mw929/Documents/CPP/KS-competition/Performance_given_policy_hyper/output/Perf_hyper_myo_valuefn_rho',num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat'];

        % Value function file
        wFile = fopen(filename2);
        w_myo = fread(wFile, precision1);
        w_myo=  reshape(w_myo,[1601,1601]);
        fclose(wFile);


        filename2 = ['/space/mw929/Documents/CPP/KS-competition/Max_prob_vertical/output/Max_prob_vertical_valuefn_rho',num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat'];

        % Value function file
        wFile = fopen(filename2);
        w_ult_ver = fread(wFile, precision1);
        w_ult_ver=  reshape(w_ult_ver,[1601,1601]);
        fclose(wFile);

        filename2 = ['/space/mw929/Documents/CPP/KS-competition/Performance_given_policy/output/Perf_const_valuefn_rho',num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat'];

        % Value function file
        wFile = fopen(filename2);
        w_const_ver = fread(wFile, precision1);
        w_const_ver=  reshape(w_const_ver,[1601,1601]);
        fclose(wFile);

        filename2 = ['/space/mw929/Documents/CPP/KS-competition/Performance_given_policy/output/Perf_myo_valuefn_rho',num2str(1000*rho),'_lamb',num2str(1000*lamb),'.dat'];

        % Value function file
        wFile = fopen(filename2);
        w_myo_ver = fread(wFile, precision1);
        w_myo_ver=  reshape(w_myo_ver,[1601,1601]);
        fclose(wFile);


      
        
        %valslice_ult = abs(w_ult(kk,:) - w_ult_ver(kk,:));
        %mean_val_ult = trapz(xx,valslice_ult);
        mean_val_ult = abs(w_ult(yindx,xindx) - w_ult_ver(yindx,xindx));
        mean_table_ult(ii,jj) = mean_val_ult;

        %valslice_const = abs(w_const(kk,:) - w_const_ver(kk,:));
        %mean_val_const = trapz(xx,valslice_const);
        mean_val_const = abs(w_const(yindx,xindx) - w_const_ver(yindx,xindx));
        mean_table_const(ii,jj) = mean_val_const;

        %valslice_myo = abs(w_myo(kk,:) - w_myo_ver(kk,:));
        %mean_val_myo = trapz(xx,valslice_myo);
        mean_val_myo = abs(w_myo(yindx,xindx) - w_myo_ver(yindx,xindx));
        mean_table_myo(ii,jj) = mean_val_myo;
        

    end
end


%% ult 
[Xmap,Ymap] = meshgrid(rho_list,lamb_list);
figure
imagesc(rho_list,lamb_list,mean_table_ult);
view([0 0 1]);
xlabel('Survival fraction ($\rho$)','FontSize',20,'Interpreter','latex');
ylabel('Dilution frequency ($\lambda$)','FontSize',20,'Interpreter','latex');
colorbar();
yticks(lamb_list);
clim([0 0.018]);
axis square
set(gca,'YDir','normal');

%% myo 
[Xmap,Ymap] = meshgrid(rho_list,lamb_list);
figure
imagesc(rho_list,lamb_list,mean_table_myo);
view([0 0 1]);
xlabel('Survival fraction ($\rho$)','FontSize',20,'Interpreter','latex');
ylabel('Dilution frequency ($\lambda$)','FontSize',20,'Interpreter','latex');
colorbar();
yticks(lamb_list);
clim([0 0.018]);
axis square
set(gca,'YDir','normal');

%% myo 
[Xmap,Ymap] = meshgrid(rho_list,lamb_list);
figure
imagesc(rho_list,lamb_list,mean_table_const);
view([0 0 1]);
xlabel('Survival fraction ($\rho$)','FontSize',20,'Interpreter','latex');
ylabel('Dilution frequency ($\lambda$)','FontSize',20,'Interpreter','latex');
colorbar();
yticks(lamb_list);
clim([0 0.018]);
axis square
set(gca,'YDir','normal');