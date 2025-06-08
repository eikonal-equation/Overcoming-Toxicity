%% Killer vs Sensitive; P(Ts < Td) with many dilutions
% Phase diagram generator before first dilution
clear;
clc;
%% Initialization
n = 800; % number of mesh points along on side
h = 1/n;


epsilon = 0.2;
rks = 0.85;
lamb = 1.0;
gamma = 1.0;
rho = 0.65;

thres_fail = 0.01;
thres_success = 1 - thres_fail;
tol = 1e-6;
delta = 0.7;


W = zeros(n+1);
A = zeros(n+1);


% velocity function for different r_k and r_s
fx = @(x,y,a)  x.*(1-x).*((1-y).*(r_ks.*(1-epsilon.*a) - 1) + a.*gamma.*x.*y);

fy = @(x,y,a)  y.*(1-y).*(1 + (r_ks.*(1-epsilon.*a) - 1).*x) - a.*gamma.*y.^2.*x.*(1-x);

% cdf of exponential distribution
ft = @(t,dt,lamb) exp(-lamb.*t) - exp(-lamb.*(t+dt));
% pdf of exponential distribution
ftt = @(t,lamb) lamb.*exp(-lamb.*t);

xc_list = [0.5,0.9,0.1,0.1];
yc_list = [0.5,0.1,0.9,0.1];
xx = linspace(0,1,n+1);
[X,Y] = meshgrid(xx,xx);
%%
rho_list = linspace(0.5,0.7,5);
lamb_list = linspace(0.8,1,5);
xloc = 0.5;
yloc = 0.1;
xindx = find(xx == xloc);
yindx = find(xx == yloc);
tic
for ii = 1:length(lamb_list)
    for jj = 1:length(rho_list)
        rho = rho_list(jj)
        lamb = lamb_list(ii)
        % computing the optimal value function / policy
        [W,A] = MaxProb_with_dilution(n,gamma,epsilon,rks,lamb,tol,delta,rho,thres_fail,thres_success);
        str = sprintf('Prob ValFunc and Policy with rate=%.2f rks=%.2f g=%.1f rho=%.2f.mat',lamb,rks,gamma,rho);
        save(str,'W','A');
    end
end
%%
[Xmap,Ymap] = meshgrid(rho_list,lamb_list);
figure
surf(Xmap,Ymap,ult_map,'EdgeColor','none');
% shading faceted;
view([0 0 1]);
xlabel('$\rho$','FontSize',20,'Interpreter','latex');
ylabel('$\lambda$','FontSize',20,'Interpreter','latex');
colorbar();
clim([0 1]);
yticks(lamb_list);
colormap("hot")
axis square
%%
str = sprintf('Performance map for xc=%g yc=%g.mat',xloc, yloc);
save(str,'const_map','det_map','rho_list','lamb_list');

%% Max Prob
ff = @(x) 1./x;
mygrey = [192 192 192]/255;

% str = sprintf('Prob ValFunc and Policy with rate=%.1f rks=%.1f g=%.1f rho=%.1f.mat',lamb,r_ks,gamma,rho);
% load(str,'W','A');

[X,Y]=meshgrid(0:h:1,0:h:1);
hh = figure('units','pixels','outerposition',[0 0 1600 800]);
subplot(1,2,1)
contourf(X,Y,A, [0, 1]);
hold on
xline(thres_success,'m:','linewidth',2);
xline(thres_fail,'m:','linewidth',2);
hold off
xlabel('Fraction of killer cells (f)');
ylabel('Total population (N)');
axis equal
title1 = sprintf('Optimal Toxin-On/Off Region');
title(title1);
colorbar;

subplot(1,2,2)
[AA,BB] = contourf(X,Y,W,15);
% [M,ctr] = contour(X,Y,W,'LineWidth',2);
hold on
xline(thres_success,'m:','linewidth',2);
xline(thres_fail,'m:','linewidth',2);
xline(0.5,'w:','LineWidth',2);
% for i = 1:30
%     plot(xx,0.02*i*ff(xx),'-.','linewidth',1.2,'Color',mygrey);
%     hold on
% end
% hold off
set(BB,'LineColor','none');
xlabel('Fraction of killer cells (f)');
ylabel('Total population (N)');
axis equal
title2 = sprintf('Value Function (probabilistic)');
% title2 = sprintf('Value Function (max E[f(T)])');
% title2 = sprintf('Minimal Time (after taking log) to reach %g%% N',thres_large*100);
title22 = sprintf('Contour plots for \x03bb=%.1f, \x03b3=%.1f, r_k/r_s=%.2f, and \x03c1=%.2f',lamb, gamma,r_ks,rho);
title(title2);
colorbar;
sgtitle(title22);
xlim([0,1]);
ylim([0,1]);
%%
num_pts = M(2,1);
ctr_x = M(1,2:num_pts + 1 - 600);
ctr_y = M(2,2:num_pts + 1 - 600);
figure
plot(ctr_y,ctr_x,'.-','linewidth',2);
% xlabel('Fraction of killer cells (f)');
% ylabel('Total population (N)');
xlabel('Total population (N)');
ylabel('Fraction of killer cells (f)');
title('On the contour $w = 0.1$','Interpreter','latex','FontSize',20);
grid on
figure
plot(ctr_y,ctr_y.*ctr_x,'.-','linewidth',2);
% xlabel('Fraction of killer cells (f)');
% ylabel('Killer population (fN)');
xlabel('Total population (N)');
ylabel('Killer population (fN)');
title('On the contour $w = 0.1$','Interpreter','latex','FontSize',20);
grid on
%% Monte Carlo Simulations
% r_ks = 0.9;
% gamma = 1.0;
% lamb_list =  0.1:0.1:0.7
num_samples = 1e5;
dt = 1e-3
xc_list = [0.5,0.5,0.5];
yc_list = [0.2,0.5,0.8];
xc = 0.5
yc = 0.5
choice = 'conservative'
A_const = ones(n+1);
mychoice = 2
tic
% for jj = 1:4
%     mychoice = jj;
%     for ii = 1:3
%         xc = xc_list(ii);
%         yc = yc_list(ii);
if mychoice == 1
    str = sprintf('Prob ValFunc and Policy with rate=%.2f rks=%.2f g=%.1f rho=%.2f.mat',lamb,r_ks,gamma,rho);
    load(str);
    [prob_success,DilutionNum_list,T_succ_list,T_fail_list,Dilution_succ_list,Dilution_fail_list] ...
        = MC_prob_dilution(num_samples,h,A,lamb,r_ks,gamma,epsilon,dt,xc,yc,choice,rho,thres_fail,thres_success);
    str = sprintf('MC prob at (x,y)=(%.1f,%.1f) with rate=%.2f rks=%.2f g=%.1f rho=%.2f.mat',xc,yc,lamb,r_ks,gamma,rho);

    save(str,'prob_success','DilutionNum_list','T_succ_list','T_fail_list','Dilution_succ_list','Dilution_fail_list');
elseif mychoice == 2
    str = sprintf('ValFunc Prob_const with rate=%.2f rks=%.2f g=%.1f rho=%.2f.mat',lamb,r_ks,gamma,rho);
    load(str);
    [prob_success,DilutionNum_list,T_succ_list,T_fail_list,Dilution_succ_list,Dilution_fail_list] ...
        = MC_prob_dilution(num_samples,h,A_const,lamb,r_ks,gamma,epsilon,dt,xc,yc,choice,rho,thres_fail,thres_success);
    str = sprintf('MC const at (x,y)=(%.1f,%.1f) with rate=%.2f rks=%.2f g=%.1f rho=%.2f.mat',xc,yc,lamb,r_ks,gamma,rho);

    save(str,'prob_success','DilutionNum_list','T_succ_list','T_fail_list','Dilution_succ_list','Dilution_fail_list');
elseif mychoice == 3
    str = sprintf('ValFunc Prob_1st with rate=%.2f rks=%.2f g=%.1f rho=%.2f.mat',lamb,r_ks,gamma,rho);
    load(str);
    [prob_success,DilutionNum_list,T_succ_list,T_fail_list,Dilution_succ_list,Dilution_fail_list] ...
        = MC_prob_dilution(num_samples,h,A_1st,lamb,r_ks,gamma,epsilon,dt,xc,yc,choice,rho,thres_fail,thres_success);
    str = sprintf('MC 1st at (x,y)=(%.1f,%.1f) with rate=%.2f rks=%.2f g=%.1f rho=%.2f.mat',xc,yc,lamb,r_ks,gamma,rho);

    save(str,'prob_success','DilutionNum_list','T_succ_list','T_fail_list','Dilution_succ_list','Dilution_fail_list');
elseif mychoice == 4
    str = sprintf('ValFunc Prob_det with rate=%.2f rks=%.2f g=%.1f rho=%.2f.mat',lamb,r_ks,gamma,rho);
    load(str);
    [prob_success,DilutionNum_list,T_succ_list,T_fail_list,Dilution_succ_list,Dilution_fail_list] ...
        = MC_prob_dilution(num_samples,h,A_det,lamb,r_ks,gamma,epsilon,dt,xc,yc,choice,rho,thres_fail,thres_success);
    str = sprintf('MC det at (x,y)=(%.1f,%.1f) with rate=%.2f rks=%.2f g=%.1f rho=%.2f.mat',xc,yc,lamb,r_ks,gamma,rho);

    save(str,'prob_success','DilutionNum_list','T_succ_list','T_fail_list','Dilution_succ_list','Dilution_fail_list');
end
%     end
% end
tt = toc

%%
xc = 0.5
yc = 0.5
mychoice = 2;
figure
DilutionNum_list_2 = DilutionNum_list(~(DilutionNum_list==-1));
% DilutionNum_list_2(DilutionNum_list == -1) = NaN;
histogram(DilutionNum_list_2,'Normalization','pdf','BinWidth',1)
xlabel('Number of Dilutions');
str1 = sprintf('Histogram of number of dilutions from (x,y) = (%g,%g) ',xc,yc);
str2 = sprintf('where \x03bb=%.1f, \x03b3=%.1f, r_k/r_s=%.1f, and \x03c1=%.1f',lamb,gamma,r_ks,rho);
title({str1,str2});
grid on
if mychoice == 1
    filename = sprintf('images_2_29/paths/num_dilution_opt_0%g0%g.png',xc*10,yc*10);
    print(filename,'-dpng');
elseif mychoice == 2
    filename = sprintf('images_2_29/paths/num_dilution_const_0%g0%g.png',xc*10,yc*10);
    print(filename,'-dpng');
elseif mychoice == 3
    filename = sprintf('images_2_29/paths/num_dilution_1st_0%g0%g.png',xc*10,yc*10);
    print(filename,'-dpng');
elseif mychoice == 4
    filename = sprintf('images_2_29/paths/num_dilution_det_0%g0%g.png',xc*10,yc*10);
    print(filename,'-dpng');
end


figure
Dilution_succ_list_2 = Dilution_succ_list(~(Dilution_succ_list==-1));
% Dilution_succ_list_2(Dilution_succ_list == -1) = NaN;
histogram(Dilution_succ_list_2,'Normalization','pdf','BinWidth',1)
xlabel('Number of Dilutions');
str1 = sprintf('Histogram of number of dilutions (victory) from (x,y) = (%g,%g) ',xc,yc);
str2 = sprintf('where \x03bb=%.1f, \x03b3=%.1f, r_k/r_s=%.1f, and \x03c1=%.1f',lamb,gamma,r_ks,rho);
title({str1,str2});
grid on
if mychoice == 1
    filename = sprintf('images_2_29/paths/victory_dilution_opt_0%g0%g.png',xc*10,yc*10);
    print(filename,'-dpng');
elseif mychoice == 2
    filename = sprintf('images_2_29/paths/victory_dilution_const_0%g0%g.png',xc*10,yc*10);
    print(filename,'-dpng');
elseif mychoice == 3
    filename = sprintf('images_2_29/paths/victory_dilution_1st_0%g0%g.png',xc*10,yc*10);
    print(filename,'-dpng');
elseif mychoice == 4
    filename = sprintf('images_2_29/paths/victory_dilution_det_0%g0%g.png',xc*10,yc*10);
    print(filename,'-dpng');
end

figure
Dilution_fail_list_2 = Dilution_fail_list(~(Dilution_fail_list==-1));
% Dilution_fail_list_2(Dilution_fail_list == -1) = NaN;
histogram(Dilution_fail_list_2,'Normalization','pdf','BinWidth',1)
xlabel('Number of Dilutions');
str1 = sprintf('Histogram of number of dilutions (defeat) from (x,y) = (%g,%g) ',xc,yc);
str2 = sprintf('where \x03bb=%.1f, \x03b3=%.1f, r_k/r_s=%.1f, and \x03c1=%.1f',lamb,gamma,r_ks,rho);
title({str1,str2});
grid on
if mychoice == 1
    filename = sprintf('images_2_29/paths/defeat_dilution_opt_0%g0%g.png',xc*10,yc*10);
    print(filename,'-dpng');
elseif mychoice == 2
    filename = sprintf('images_2_29/paths/defeat_dilution_const_0%g0%g.png',xc*10,yc*10);
    print(filename,'-dpng');
elseif mychoice == 3
    filename = sprintf('images_2_29/paths/defeat_dilution_1st_0%g0%g.png',xc*10,yc*10);
    print(filename,'-dpng');
elseif mychoice == 4
    filename = sprintf('images_2_29/paths/defeat_dilution_det_0%g0%g.png',xc*10,yc*10);
    print(filename,'-dpng');
end

figure
T_succ_list_2 = T_succ_list(~(T_succ_list == 1e6));
% T_succ_list_2(T_succ_list == 1e6) = NaN;
pts = linspace(min(T_succ_list_2),max(T_succ_list_2),201);
[fn,xi] = ksdensity(T_succ_list_2,pts,'Kernel','normal');
hold on
histogram(T_succ_list_2,'Normalization','pdf')
plot(xi,fn,'linewidth',2);
hold off
str1 = sprintf('EPDF of T_s starting from (x,y) = (%g,%g)',xc,yc);
str2 = sprintf('where \x03bb=%.1f, \x03b3=%.1f, r_k/r_s=%.1f, and \x03c1=%.1f',lamb,gamma,r_ks,rho);
title({str1,str2});
xlabel('Time (t)');
legend('histogram','ksdensity','Fontsize',15)
grid on
if mychoice == 1
    filename = sprintf('images_2_29/paths/victory_pdf_opt_0%g0%g.png',xc*10,yc*10);
    print(filename,'-dpng');
elseif mychoice == 2
    filename = sprintf('images_2_29/paths/victory_pdf_const_0%g0%g.png',xc*10,yc*10);
    print(filename,'-dpng');
elseif mychoice == 3
    filename = sprintf('images_2_29/paths/victory_pdf_1st_0%g0%g.png',xc*10,yc*10);
    print(filename,'-dpng');
elseif mychoice == 4
    filename = sprintf('images_2_29/paths/victory_pdf_det_0%g0%g.png',xc*10,yc*10);
    print(filename,'-dpng');
end


T_fail_list_2 = T_fail_list(~(T_fail_list == 1e6));
figure
pts = linspace(min(T_fail_list_2),max(T_fail_list_2),201);
[fn,xi] = ksdensity(T_fail_list_2,pts,'Kernel','normal');
hold on
histogram(T_fail_list_2,'Normalization','pdf')
plot(xi,fn,'linewidth',2);
hold off
str1 = sprintf('EPDF of T_d starting from (x,y) = (%g,%g)',xc,yc);
str2 = sprintf('where \x03bb=%.1f, \x03b3=%.1f, r_k/r_s=%.1f, and \x03c1=%.1f',lamb,gamma,r_ks,rho);
title({str1,str2});
xlabel('Time (t)');
legend('histogram','ksdensity','Fontsize',15)
grid on
if mychoice == 1
    filename = sprintf('images_2_29/paths/defeat_pdf_opt_0%g0%g.png',xc*10,yc*10);
    print(filename,'-dpng');
elseif mychoice == 2
    filename = sprintf('images_2_29/paths/defeat_pdf_const_0%g0%g.png',xc*10,yc*10);
    print(filename,'-dpng');
elseif mychoice == 3
    filename = sprintf('images_2_29/paths/defeat_pdf_1st_0%g0%g.png',xc*10,yc*10);
    print(filename,'-dpng');
elseif mychoice == 4
    filename = sprintf('images_2_29/paths/defeat_pdf_det_0%g0%g.png',xc*10,yc*10);
    print(filename,'-dpng');
end
%% sample path visualization
xc = 0.5;
yc = 0.5;
choice = 'conservative';
dt = 1e-3;
mychoice = 2;
[X,Y] = meshgrid(0:h:1,0:h:1);
if mychoice == 1
    [xlist_full,ylist_full,policy_list,T_success,T_fail,total_dilution] ...
        = Optimal_path_ProbDilution(h,A,lamb,r_ks,gamma,epsilon,dt,xc,yc,choice,rho,thres_fail,thres_success);
    figure
    contourf(X,Y,A,[0,1]);
    hold on
    plot(xlist_full(policy_list==1),ylist_full(policy_list==1),'r.');
    plot(xlist_full(policy_list==0),ylist_full(policy_list==0),'g.');
    plot(xlist_full(1),ylist_full(1),'p','markersize',7,'linewidth',1.5,'markerfacecolor','w','MarkerEdgeColor','w');
    xline(thres_success,'m:','linewidth',2);
    xline(thres_fail,'m:','linewidth',2);
    xlabel('Fraction of killer cells (f)');
    ylabel('Total population (N)');
    axis equal
    xlim([0,1]);
    ylim([0,1]);
elseif mychoice == 2
    [xlist_full,ylist_full,policy_list,T_success,T_fail,total_dilution] ...
        = Optimal_path_ProbDilution(h,A_const,lamb,r_ks,gamma,epsilon,dt,xc,yc,choice,rho,thres_fail,thres_success);
    figure
    contourf(X,Y,A_const,[0,1]);
    hold on
    plot(xlist_full(policy_list==1),ylist_full(policy_list==1),'r.');
    plot(xlist_full(policy_list==0),ylist_full(policy_list==0),'g.');
    plot(xlist_full(1),ylist_full(1),'p','markersize',7,'linewidth',1.5,'markerfacecolor','w','MarkerEdgeColor','w');
    xline(thres_success,'m:','linewidth',2);
    xline(thres_fail,'m:','linewidth',2);
    xlabel('Fraction of killer cells (f)');
    ylabel('Total population (N)');
    axis equal
    xlim([0,1]);
    ylim([0,1]);
elseif mychoice == 3
    [xlist_full,ylist_full,policy_list,T_success,T_fail,total_dilution] ...
        = Optimal_path_ProbDilution(h,A_1st,lamb,r_ks,gamma,epsilon,dt,xc,yc,choice,rho,thres_fail,thres_success);

    figure
    contourf(X,Y,A_1st,[0,1]);
    hold on
    plot(xlist_full(policy_list==1),ylist_full(policy_list==1),'r.');
    plot(xlist_full(policy_list==0),ylist_full(policy_list==0),'g.');
    plot(xlist_full(1),ylist_full(1),'p','markersize',7,'linewidth',1.5,'markerfacecolor','w','MarkerEdgeColor','w');
    xline(thres_success,'m:','linewidth',2);
    xline(thres_fail,'m:','linewidth',2);
    xlabel('Fraction of killer cells (f)');
    ylabel('Total population (N)');
    axis equal
    xlim([0,1]);
    ylim([0,1]);
elseif mychoice == 4
    [xlist_full,ylist_full,policy_list,T_success,T_fail,total_dilution] ...
        = Optimal_path_ProbDilution(h,A_det,lamb,r_ks,gamma,epsilon,dt,xc,yc,choice,rho,thres_fail,thres_success);

    figure
    contourf(X,Y,A_det,[0,1]);
    hold on
    plot(xlist_full(policy_list==1),ylist_full(policy_list==1),'r.');
    plot(xlist_full(policy_list==0),ylist_full(policy_list==0),'g.');
    plot(xlist_full(1),ylist_full(1),'p','markersize',7,'linewidth',1.5,'markerfacecolor','w','MarkerEdgeColor','w');
    xline(thres_success,'m:','linewidth',2);
    xline(thres_fail,'m:','linewidth',2);
    xlabel('Fraction of killer cells (f)');
    ylabel('Total population (N)');
    axis equal
    xlim([0,1]);
    ylim([0,1]);
end


%% Max E[f(T)]
str = sprintf('ValFunc and Policy with rate=%.1f rks=%.1f g=%.1f.mat',lamb,r_ks,gamma);
load(str,'U','A');
%%
[U,A] = MaxFreq_expT(n,gamma,epsilon,r_ks,lamb,tol,delta);
str = sprintf('ValFunc and Policy with rate=%.1f rks=%.2f g=%.1f.mat',lamb,r_ks,gamma);
save(str,'U','A');


[X,Y]=meshgrid(0:h:1,0:h:1);
hh = figure('units','pixels','outerposition',[0 0 1600 800]);
subplot(1,2,1)
contourf(X,Y,A, [0, 1]);
hold on
% xline(thres_success,'m:','linewidth',2);
% xline(thres_fail,'m:','linewidth',2);
hold off
xlabel('Fraction of killer cells (f)');
ylabel('Total population (N)');
axis equal
title1 = sprintf('Optimal Toxin-On/Off Region');
title(title1);
colorbar;

subplot(1,2,2)
[AA,BB] = contourf(X,Y,U);
hold on
% xline(thres_success,'m:','linewidth',2);
% xline(thres_fail,'m:','linewidth',2);
xline(0.5,'w:','LineWidth',2);
hold off
set(BB,'LineColor','none');
xlabel('Fraction of killer cells (f)');
ylabel('Total population (N)');
axis equal
title2 = sprintf('Value Function (max E[f])');
% title2 = sprintf('Minimal Time (after taking log) to reach %g%% N',thres_large*100);
title22 = sprintf('Contour plots for \x03bb=%.1f, \x03b3=%.1f, and r_k/r_s=%.2f',lamb, gamma,r_ks);
title(title2);
colorbar;
sgtitle(title22);
colorbar;

%% Prob and constitutive killer
% rho_list = [linspace(0.5,0.59,10),linspace(0.6,0.7,5)];
% rho_list = linspace(0.51,0.59,9);
rho_list = [0.52,0.525,0.53,0.535,0.54,0.545,0.55,0.555,0.56,0.565,0.57,0.575,0.58,0.585,0.59,0.595,0.6,0.65,0.7];
% lamb_list = linspace(0.8,1,5);
% lamb_list = [0.7 0.825 0.875 0.925 0.975];
lamb_list = 0.75;
% const_map = zeros(length(lamb_list),length(rho_list));
xloc = 0.5;
yloc = 0.1;
xindx = find(xx == xloc);
yindx = find(xx == yloc);
A_const = constitutive_policy(n,thres_fail,thres_success);
tic
for ii = 1:length(lamb_list)
    for jj = 1:length(rho_list)
        rho = rho_list(jj)
        lamb = lamb_list(ii)

        W_const = Prob_given_policy(A_const,n,gamma,epsilon,rks,lamb,rho,thres_fail,thres_success);

        str = sprintf('ValFunc Prob_const with rate=%.3f rks=%.2f g=%.1f rho=%.3f.mat',lamb,rks,gamma,rho);
        save(str,'W_const','A_const');

        %     temp = W_const(yindx,xindx);
        %     const_map(ii,jj) = temp;

    end
end
tt = toc

%%
save('heatmap_const.mat','const_map');

%%
[Xmap,Ymap] = meshgrid(rho_list,lamb_list);
figure
imagesc(rho_list,lamb_list,const_map);
view([0 0 1]);
xlabel('$\rho$','FontSize',20,'Interpreter','latex');
ylabel('$\lambda$','FontSize',20,'Interpreter','latex');
colorbar();
yticks(lamb_list);
clim([0 1]);
colormap("hot")
axis square
set(gca,'YDir','normal');
%%
[X,Y]=meshgrid(0:h:1,0:h:1);
hh = figure('units','pixels','outerposition',[0 0 1600 800]);
subplot(1,2,1)
contourf(X,Y,A_const, [0, 1]);
hold on
xline(thres_success,'m:','linewidth',2);
xline(thres_fail,'m:','linewidth',2);
hold off
xlabel('Fraction of killer cells (f)');
ylabel('Total population (N)');
axis equal
title1 = sprintf('Toxin-On/Off Region');
title(title1);
colorbar;

subplot(1,2,2)
[AA,BB] = contourf(X,Y,W_const);
hold on
xline(thres_success,'m:','linewidth',2);
xline(thres_fail,'m:','linewidth',2);
xline(0.5,'w:','LineWidth',2);
hold off
set(BB,'LineColor','none');
xlabel('Fraction of killer cells (f)');
ylabel('Total population (N)');
axis equal
title2 = sprintf('Performance of the policy (Probabilistic)');
% title2 = sprintf('Minimal Time (after taking log) to reach %g%% N',thres_large*100);
title22 = sprintf('Contour plots for Constitutive Killer \x03bb=%.1f, \x03b3=%.1f, r_k/r_s=%.2f, and \x03c1=%.2f',lamb, gamma,r_ks,rho);
title(title2);
colorbar;
sgtitle(title22);
colorbar;

%%
[AA,BB] = contourf(X,Y,W_const,12);
hold on
xline(thres_success,'m:','linewidth',2);
xline(thres_fail,'m:','linewidth',2);
% xline(0.5,'w:','LineWidth',2);
hold off
set(BB,'LineColor','none');
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
axis equal
colorbar();
clim([0 1]);
cc = linspace(0,1,11);
set(hcb,'Ticks',cc);
myCticks = arrayfun(@(x) sprintf('%.1f',x),cc,'un',0);
set(hcb,'Ticklabels',myCticks);

%% Prob and deterministic policy
rho_list = linspace(0.5,0.7,5);
lamb_list = linspace(0.8,1,5);
const_map = zeros(length(lamb_list),length(rho_list));
xloc = 0.5;
yloc = 0.1;
xindx = find(xx == xloc);
yindx = find(xx == yloc);
det_map = zeros(length(lamb_list),length(rho_list));
[U,A] = min_time_VPI(n,gamma,epsilon,r_ks,tol,delta,thres_fail,thres_success);
A_det = A;
for ii = 1:length(lamb_list)
    for jj = 1:length(rho_list)
        rho = rho_list(jj)
        lamb = lamb_list(ii)


        W_det = Prob_given_policy(A_det,n,gamma,epsilon,r_ks,lamb,rho,thres_fail,thres_success);

        str = sprintf('ValFunc Prob_det with rate=%.1f rks=%.2f g=%.1f rho=%.2f.mat',lamb,r_ks,gamma,rho);
        save(str,'W_det','A_det');
        det_map(ii,jj) = W_det(yindx,xindx);

    end
end
%%
[Xmap,Ymap] = meshgrid(rho_list,lamb_list);
figure
surf(Xmap,Ymap,det_map,'EdgeColor','none');
% shading faceted;
view([0 0 1]);
xlabel('$\rho$','FontSize',20,'Interpreter','latex');
ylabel('$\lambda$','FontSize',20,'Interpreter','latex');
colorbar();
clim([0 1]);
yticks(lamb_list);
colormap("hot")
axis square

%%
[X,Y]=meshgrid(0:h:1,0:h:1);
hh = figure('units','pixels','outerposition',[0 0 1600 800]);
subplot(1,2,1)
contourf(X,Y,A_det, [0, 1]);
hold on
xline(thres_success,'m:','linewidth',2);
xline(thres_fail,'m:','linewidth',2);
hold off
xlabel('Fraction of killer cells (f)');
ylabel('Total population (N)');
axis equal
title1 = sprintf('Toxin-On/Off Region');
title(title1);
colorbar;

subplot(1,2,2)
[AA,BB] = contourf(X,Y,W_det);
hold on
xline(thres_success,'m:','linewidth',2);
xline(thres_fail,'m:','linewidth',2);
xline(0.5,'w:','LineWidth',2);
hold off
set(BB,'LineColor','none');
xlabel('Fraction of killer cells (f)');
ylabel('Total population (N)');
axis equal
title2 = sprintf('Performance of the policy (Probabilistic)');
% title2 = sprintf('Minimal Time (after taking log) to reach %g%% N',thres_large*100);
title22 = sprintf('Contour plots with deterministic policy \x03bb=%.1f, \x03b3=%.1f, r_k/r_s=%.2f, and \x03c1=%.2f',lamb, gamma,r_ks,rho);
title(title2);
colorbar;
sgtitle(title22);
colorbar;

%%
N=20;
NN = 800;
hh = 1/NN;
q=linspace(0,1,N+1);
p=linspace(0,1,N+1);
qq=linspace(0,1,NN+1);
pp=linspace(0,1,NN+1);
[x,y]=meshgrid(q',p');

A_const_down = A_const(1:NN/N:end,1:NN/N:end);

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
contourf(X,Y,A_det, [0, 1]);
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
xline(thres_success,'m:','linewidth',2);
xline(thres_fail,'m:','linewidth',2);
hold off
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
axis equal


%%
[AA,BB] = contourf(X,Y,W_det,12);
hold on
xline(thres_success,'m:','linewidth',2);
xline(thres_fail,'m:','linewidth',2);
% xline(0.5,'w:','LineWidth',2);
hold off
set(BB,'LineColor','none');
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
axis equal
colorbar();
clim([0 1]);
cc = linspace(0,1,11);
set(hcb,'Ticks',cc);
myCticks = arrayfun(@(x) sprintf('%.1f',x),cc,'un',0);
set(hcb,'Ticklabels',myCticks);

%% Prob and E[f(T)]
str = sprintf('ValFunc and Policy with rate=%.1f rks=%.2f g=%.1f.mat',lamb,r_ks,gamma);
load(str,'U','A');
A_1st = A;
W_1st = Prob_given_policy(A_1st,n,gamma,epsilon,r_ks,lamb,rho,thres_fail,thres_success);


str = sprintf('ValFunc Prob_1st with rate=%.1f rks=%.2f g=%.1f rho=%.2f.mat',lamb,r_ks,gamma,rho);
save(str,'W_1st','A_1st');
%%
% rho_list = linspace(0.5,0.7,5);
%rho_list = [0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.65, 0.7];
rho_list = [0.525, 0.535, 0.545, 0.555, 0.565, 0.575, 0.585, 0.595];
lamb_list = [0.75, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0];

xloc = 0.5;
yloc = 0.1;
xindx = find(xx == xloc);
yindx = find(xx == yloc);

% myo_map = zeros(length(lamb_list),length(rho_list));
tic
for ii = 1:length(lamb_list)

    lamb = lamb_list(ii)

    [U_1st,A_1st] = MaxFreq_expT(n,gamma,epsilon,rks,lamb,tol,delta);
    str = sprintf('ValFunc and Policy with rate=%.3f rks=%.2f g=%.1f.mat',lamb,rks,gamma);
    save(str,'U_1st','A_1st');

    for jj = 1:length(rho_list)

        rho = rho_list(jj)



        W_1st = Prob_given_policy(A_1st,n,gamma,epsilon,rks,lamb,rho,thres_fail,thres_success);

        str = sprintf('ValFunc Prob_1st with rate=%.3f rks=%.2f g=%.1f rho=%.3f.mat',lamb,rks,gamma,rho);
        save(str,'W_1st','A_1st');

        %         myo_map(ii,jj) = W_1st(yindx,xindx);

        ii
        jj

    end
end
tt = toc
%%
[Xmap,Ymap] = meshgrid(rho_list,lamb_list);
figure
surf(Xmap,Ymap,myo_map,'EdgeColor','none');
% shading faceted;
view([0 0 1]);
xlabel('$\rho$','FontSize',20,'Interpreter','latex');
ylabel('$\lambda$','FontSize',20,'Interpreter','latex');
colorbar();
clim([0 1]);
yticks(lamb_list);
colormap("hot")
axis square
%%
save('heatmap_myo.mat','myo_map');

%%
[X,Y]=meshgrid(0:h:1,0:h:1);
hh = figure('units','pixels','outerposition',[0 0 1600 800]);
subplot(1,2,1)
contourf(X,Y,A_1st, [0, 1]);
hold on
xline(thres_success,'m:','linewidth',2);
xline(thres_fail,'m:','linewidth',2);
hold off
xlabel('Fraction of killer cells (f)');
ylabel('Total population (N)');
axis equal
title1 = sprintf('Toxin-On/Off Region');
title(title1);
colorbar;

subplot(1,2,2)
[AA,BB] = contourf(X,Y,W_1st);
hold on
xline(thres_success,'m:','linewidth',2);
xline(thres_fail,'m:','linewidth',2);
xline(0.5,'w:','LineWidth',2);
hold off
set(BB,'LineColor','none');
xlabel('Fraction of killer cells (f)');
ylabel('Total population (N)');
axis equal
title2 = sprintf('Performance of the policy (Probabilistic)');
% title2 = sprintf('Minimal Time (after taking log) to reach %g%% N',thres_large*100);
title22 = sprintf('Contour plots with 1st dilution policy \x03bb=%.1f, \x03b3=%.1f, r_k/r_s=%.2f, and \x03c1=%.2f',lamb, gamma,r_ks,rho);
title(title2);
colorbar;
sgtitle(title22);
colorbar;
%%
% W(W < 1e-4) = NaN;
% W_const(W_const < 1e-4) = NaN;
% W_1st(W_1st < 1e-4) = NaN;
% W_det(W_det < 1e-4) = NaN;

% W1 = log(W - W_const)./W);
% W2 = log((W - W_det)./W);
% W3 = log((W - W_1st)./W);
% W1 = (W - W_const)./W;
% W2 = (W - W_det)./W;
% W3 = (W - W_1st)./W;




W1 = (W - W_const);
% W2 = (W - W_det);
W3 = (W - W_1st);

figure
[AA,BB] = contourf(X,Y,W1,5,"ShowText",true,"LabelFormat","%.2f");
hold on
xline(thres_success,'m:','linewidth',2);
xline(thres_fail,'m:','linewidth',2);
% xline(0.5,'w:','LineWidth',2);
hold off
% set(BB,'LineColor','none');
xlabel('Fraction of killer cells (f)','FontSize',14);
ylabel('Total population (N)','FontSize',14);
axis equal
colorbar();
% clim([0,0.7]);

% figure
% [AA,BB] = contourf(X,Y,W2,12);
% hold on
% xline(thres_success,'m:','linewidth',2);
% xline(thres_fail,'m:','linewidth',2);
% xline(0.5,'w:','LineWidth',2);
% hold off
% set(BB,'LineColor','none');
% xlabel('Fraction of killer cells (f)','FontSize',14);
% ylabel('Total population (N)','FontSize',14);
% axis equal
% colorbar();
% % clim([0,0.007]);

figure
[AA,BB] = contourf(X,Y,W3,5,"ShowText",true,"LabelFormat","%.3f");
hold on
xline(thres_success,'m:','linewidth',2);
xline(thres_fail,'m:','linewidth',2);
% xline(0.5,'w:','LineWidth',2);
hold off
% set(BB,'LineColor','none');
xlabel('Fraction of killer cells (f)','FontSize',14);
ylabel('Total population (N)','FontSize',14);
axis equal
colorbar();



%%
% velocity function for different r_k and r_s
fx = @(x,y,a)  x.*(1-x).*((1-y).*(rks.*(1-epsilon.*a) - 1) + a.*gamma.*x.*y);

fy = @(x,y,a)  y.*(1-y).*(1 + (rks.*(1-epsilon.*a) - 1).*x) - a.*gamma.*y.^2.*x.*(1-x);


N=20;
NN = 800;
hh = 1/NN;
q=linspace(0,1,N+1);
p=linspace(0,1,N+1);
qq=linspace(0,1,NN+1);
pp=linspace(0,1,NN+1);
[x,y]=meshgrid(q',p');

A = ones(NN+1);
A_down = A(1:NN/N:end,1:NN/N:end);

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
contourf(X,Y,A, [0, 1],'EdgeColor','none');
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
            if A_down(ii,ij) == 0
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
%%
A = ones(n+1);
A(1,1) = 0;
% lamb = 1;
rho = 0.65;
%rho = 0.001015;
num_dilution = 1;
dt = h;
xc = 0.1;
yc = 1;
choice = 'conservative'
[xlist_full,ylist_full,policy_list,xstart, ystart, xend, yend] ...
    = Optimal_path_DilutionIndexed(num_dilution,h,A,lamb,rks,gamma,epsilon,dt,xc,yc,choice,rho);
xx = 0:h:1;
[X,Y] = meshgrid(xx,xx);

%%
figure
contourf(X,Y,A,[0,1],'EdgeColor','none');
hold on

plot(xlist_full,ylist_full,'.','LineWidth',3)
hold off
axis equal
%%
f_lim = @(r,rho,T) 1./(rho+(1-rho)./(1-exp(-(r.*T+log(rho)))));
mypink = [148 0 211]/255;
figure
contourf(X,Y,A,[0,1],'EdgeColor','none');
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
            if A_down(ii,ij) == 0
                set(ah,'position',[x(ii,ij) y(ii,ij) 0.03*u(ii,ij) 0.03*v(ii,ij)]);
            else
                set(ah,'position',[x(ii,ij) y(ii,ij) 0.03*uu(ii,ij) 0.03*vv(ii,ij)]);
            end
        end
    end
end

plot(xlist_full(policy_list==1),ylist_full(policy_list==1),'.','markersize',3,'Color',mypink);
plot(xlist_full(policy_list==0),ylist_full(policy_list==0),'g.','markersize',3);
% plot(xstart(1),ystart(1),'o','markersize',4.5,'linewidth',1.5,'markerfacecolor','w','MarkerEdgeColor','w');
% plot(xend(1),yend(1),'s','markersize',4,'linewidth',1.5,'markerfacecolor','k','MarkerEdgeColor','k');
% %
% %
% plot(xstart(2),ystart(2),'o','markersize',4,'linewidth',1.5,'markerfacecolor','c','MarkerEdgeColor','c');
% plot(xend(2),yend(2),'s','markersize',4,'linewidth',1.5,'markerfacecolor','k','MarkerEdgeColor','k');
% %
% %
% plot(xstart(3),ystart(3),'o','markersize',4,'linewidth',1.5,'markerfacecolor','c','MarkerEdgeColor','c');
% plot(xend(3),yend(3),'s','markersize',4,'linewidth',1.5,'markerfacecolor','k','MarkerEdgeColor','k');
% plot(xstart(4),ystart(4),'o','markersize',4,'linewidth',1.5,'markerfacecolor','c','MarkerEdgeColor','c');
% plot(xend(4),yend(4),'s','markersize',4,'linewidth',1.5,'markerfacecolor','k','MarkerEdgeColor','k');
% plot(xstart(5),ystart(5),'o','markersize',4,'linewidth',1.5,'markerfacecolor','c','MarkerEdgeColor','c');
% plot(xend(5),yend(5),'s','markersize',4,'linewidth',1.5,'markerfacecolor','k','MarkerEdgeColor','k');
% plot(xstart(6),ystart(6),'o','markersize',4,'linewidth',1.5,'markerfacecolor','c','MarkerEdgeColor','c');
% plot(xend(6),yend(6),'s','markersize',4,'linewidth',1.5,'markerfacecolor','k','MarkerEdgeColor','k');

% % Regular kilers win (0.5,0.7)
% ylast_end = f_lim(rks*(1-epsilon),rho,1);
% ylast_init = ylast_end*rho;
% ylast = ylast_init:0.001:ylast_end;
% xlast = 0.998*ones(1,length(ylast));
%
% plot(xlast,ylast,'.-','markersize',4,'linewidth',3.5,'Color',mypink);
% plot(xlast(1),ylast_init,'o','markersize',4.9,'linewidth',1.5,'markerfacecolor','c','MarkerEdgeColor','c');
% plot(xlast(1),ylast_end,'s','markersize',4.9,'linewidth',1.5,'markerfacecolor','k','MarkerEdgeColor','k');
%
% annotation('arrow',[0.563,0.563],[0.705,0.515],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',7);
% annotation('arrow',[0.588,0.588],[0.595,0.445],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',7);
% annotation('arrow',[0.606,0.606],[0.555,0.415],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',7);
% annotation('arrow',[0.821,0.821],[0.44,0.45],'color','k','linewidth',1.3,'linestyle','--','headwidth',13,'headlength',6);
% annotation('arrow',[0.821,0.821],[0.39,0.38],'color','k','linewidth',1.3,'linestyle','--','headwidth',13,'headlength',6);
% text(0.685,0.38,'The process','fontsize',12);
% text(0.685,0.33,'continues...','fontsize',12);
% text(xstart(1)-0.03,ystart(1)-0.03,'Init','fontsize',12);
% text(xstart(2)-0.02,ystart(2)-0.04,'1','fontsize',12);
% text(xstart(3)-0.02,ystart(3)-0.04,'2','fontsize',12);
% text(xstart(4)-0.02,ystart(4)-0.04,'3','fontsize',12);


% % Regular kilers lose (0.5,0.1)
% ylast_end = f_lim(1,rho,1);
% ylast_init = ylast_end*rho;
% ylast = ylast_init:0.001:ylast_end;
% xlast = 0.001*ones(1,length(ylast));
%
% plot(xlast,ylast,'.-','markersize',4,'linewidth',3.5,'Color',mypink);
% plot(xlast(1),ylast_init,'o','markersize',4.9,'linewidth',1.5,'markerfacecolor','c','MarkerEdgeColor','c');
% plot(xlast(1),ylast_end,'s','markersize',4.9,'linewidth',1.5,'markerfacecolor','k','MarkerEdgeColor','k');
%
% annotation('arrow',[xend(1)+0.039,xend(1)+0.039],[yend(1)+0.08,ystart(2)+0.095],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',7);
% annotation('arrow',[xend(2)+0.055,xend(2)+0.055],[yend(2)+0.07,ystart(3)+0.09],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',7);
% annotation('arrow',[xend(3)+0.071,xend(3)+0.071],[yend(3)+0.05,ystart(4)+0.085],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',7);
% annotation('arrow',[0.213,0.213],[0.63,0.64],'color','k','linewidth',1.3,'linestyle','--','headwidth',13,'headlength',6);
% annotation('arrow',[0.213,0.213],[0.515,0.505],'color','k','linewidth',1.3,'linestyle','--','headwidth',13,'headlength',6);
% text(0.06,0.4,'The process','fontsize',12);
% text(0.06,0.35,'continues...','fontsize',12);
% text(xstart(1)-0.03,ystart(1)-0.03,'Init','fontsize',12);
% text(xstart(2)-0.02,ystart(2)-0.04,'1','fontsize',12);
% text(xstart(3)-0.02,ystart(3)-0.04,'2','fontsize',12);
% text(xstart(4)-0.02,ystart(4)-0.04,'3','fontsize',12);


% % Random unlucky (0.5,0.1)
% annotation('arrow',[xend(1)+0.04,xend(1)+0.04],[yend(1)+0.07,ystart(2)+0.096],'linewidth',1.3,'linestyle','--','headwidth',7,'headlength',6);
% annotation('arrow',[xend(2)+0.05,xend(2)+0.05],[yend(2)+0.069,ystart(3)+0.092],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',6);
% annotation('arrow',[xend(3)+0.069,xend(3)+0.069],[yend(3)+0.05,ystart(4)+0.085],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',7);
% annotation('arrow',[xend(4)+0.084,xend(4)+0.084],[yend(4)+0.03,ystart(5)+0.07],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',7);
% text(0.02,0.35,'The process','fontsize',12);
% text(0.02,0.28,'continues...','fontsize',12);
% text(xstart(1)-0.02,ystart(1)-0.03,'Init','fontsize',12);
% text(xstart(2)-0.01,ystart(2)-0.035,'1','fontsize',12);
% text(xstart(3)-0.02,ystart(3)-0.04,'2','fontsize',12);
% text(xstart(4)-0.02,ystart(4)-0.04,'3','fontsize',12);
% text(xstart(5)-0.02,ystart(5)-0.04,'4','fontsize',12);


% % Random common (0.5,0.1)
% annotation('arrow',[xend(1)+0.044,xend(1)+0.044],[yend(1)+0.07,ystart(2)+0.093],'linewidth',1.3,'linestyle','--','headwidth',7,'headlength',6);
% annotation('arrow',[xend(2)+0.0575,xend(2)+0.0575],[yend(2)+0.055,ystart(3)+0.090],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',6);
% annotation('arrow',[xend(3)+0.021,xend(3)+0.021],[yend(3)-0.025,ystart(4)+0.03],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',7);
% annotation('arrow',[xend(4)-0.0025,xend(4)-0.0025],[yend(4)-0.01,ystart(5)+0.045],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',7);
% text(0.6,0.35,'The process','fontsize',12);
% text(0.6,0.28,'continues...','fontsize',12);
% text(xstart(1)-0.02,ystart(1)-0.03,'Init','fontsize',12);
% text(xstart(2)-0.01,ystart(2)-0.035,'1','fontsize',12);
% text(xstart(3)-0.02,ystart(3)-0.04,'2','fontsize',12);
% text(xstart(4)-0.02,ystart(4)-0.04,'3','fontsize',12);
% text(xstart(5)-0.02,ystart(5)-0.04,'4','fontsize',12);


% % Random lucky (0.5,0.1)
% annotation('arrow',[xend(1)+0.0018,xend(1)+0.0018],[yend(1)-0.031,ystart(2)+0.04],'linewidth',1.3,'linestyle','--','headwidth',7,'headlength',6);
% annotation('arrow',[xend(2)-0.063,xend(2)-0.063],[yend(2)-0.032,ystart(3)+0.032],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',6);
% annotation('arrow',[xend(3)-0.138,xend(3)-0.138],[yend(3)-0.03,ystart(4)+0.03],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',7);
% annotation('arrow',[xend(4)-0.161,xend(4)-0.161],[yend(4)-0.03,ystart(5)+0.03],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',7);
% text(0.7,0.35,'The process','fontsize',12);
% text(0.7,0.28,'continues...','fontsize',12);
% text(xstart(1)-0.02,ystart(1)-0.03,'Init','fontsize',12);
% text(xstart(2)-0.01,ystart(2)-0.035,'1','fontsize',12);
% text(xstart(3)-0.02,ystart(3)-0.04,'2','fontsize',12);
% text(xstart(4)-0.02,ystart(4)-0.04,'3','fontsize',12);
% text(xstart(5)-0.015,ystart(5)-0.04,'4','fontsize',12);


xlabel('Fraction of the killer (f)','fontsize',14);
% ylabel('Total population (N)','fontsize',14);
axis equal

xlim([0,1]);
ylim([0,1]);
colormap("copper")

%%
figure
nk_win = xlist_full.*ylist_full;
ns_win = (1-xlist_full).*ylist_full;
nn = length(xlist_full);
t_list = dt.*(1:nn);
plot(t_list,nk_win,'b-','linewidth',1.5);
hold on
plot(t_list,ns_win,'r-','linewidth',1.5);
hold off
xlabel("Time (t)",'fontsize',18);
ylabel("Normalized population",'fontsize',18);

yyaxis right
plot(t_list,xlist_full,'-.','linewidth',2.5);
ylabel("Frequency of the killer",'fontsize',18);


legend('Constitutive killer','Sensitive','freq of the killer (f)','fontsize',15,'location','best');

title('Killers win','fontsize',20)
grid on
grid minor

f_win = xlist_full;
N_win = ylist_full;

save('killer win.mat','nk_win','ns_win','f_win','N_win');

%%
figure
nk_lose = xlist_full.*ylist_full;
ns_lose = (1-xlist_full).*ylist_full;
nn = length(xlist_full);
t_list = dt.*(1:nn);
plot(t_list,nk_lose,'b-','linewidth',1.5);
hold on
plot(t_list,ns_lose,'r-','linewidth',1.5);
hold off
xlabel("Time (t)",'fontsize',18);
ylabel("Normalized population",'fontsize',18);


yyaxis right
plot(t_list,xlist_full,'-.','linewidth',2.5);
ylabel("Frequency of the killer",'fontsize',18);

legend('Constitutive killer','Sensitive','freq of the killer (f)','fontsize',15,'location','best');
title('Sensitives win','fontsize',15)
grid on
grid minor


f_lose = xlist_full;
N_lose = ylist_full;

save('killer lose.mat','nk_lose','ns_lose','f_lose','N_lose');

%%
mygreen = [50,205,50]/255;
nn = length(nk_win);
t_list = dt.*(1:nn);
% figure
tiledlayout(2,1)
nexttile
hold on
plot(t_list,nk_win,'r-','linewidth',1.5);
plot(t_list,ns_win,'b-','linewidth',1.5);
hold off
% xlabel("Time (t)",'fontsize',18);
% ylabel("Normalized population",'fontsize',18);

yyaxis right
plot(t_list,f_win,'-.','linewidth',2.5,'color',mygreen);

legend('Constitutive killer ($n_K$)','Sensitive ($n_S$)','Fraction of the killer ($f$)',...
    'fontsize',12,'location','best','interpreter','latex');
title('Killers win','fontsize',15)
grid on
% grid minor
xlim([0,82]);
ylim([0.5,1.01]);
ax = gca;
ax.YAxis(2).Color = mygreen;


nexttile
hold on
plot(t_list,nk_lose,'r-','linewidth',1.5);
plot(t_list,ns_lose,'b-','linewidth',1.5);
hold off
xlabel("Rescaled time (t)",'fontsize',18);
% ylabel("Normalized population",'fontsize',18);

yyaxis right
plot(t_list,f_lose,'-.','linewidth',2.5,'color',mygreen);

% legend('Constitutive killer','Sensitive','Freq of the killer (f)','fontsize',15,'location','best');
title('Sensitives win','fontsize',15)
grid on
% grid minor
xlim([0,82]);
ylim([0,0.44]);
ax = gca;
ax.YAxis(2).Color = mygreen;

ttt = text(-9,0.1,'Normalized population','fontsize',18);
ttt2 = text(89,0.17,'Fraction of the killer','fontsize',18,'Color',mygreen);
ttt.Rotation = 90;
ttt2.Rotation = 90;


%%
figure
nk = xlist_full.*ylist_full;
ns = (1-xlist_full).*ylist_full;
nn = length(xlist_full);
t_list = dt.*(1:nn);
plot(t_list,nk,'r-','linewidth',2.5);
hold on
plot(t_list,ns,'b-','linewidth',2.5);
hold off
xlabel("Rescaled time (t)",'fontsize',15);
ylabel("Normalized population",'fontsize',15);


legend('Toxin-producer','Toxin-sensitive','fontsize',15);
% title('Killers win','fontsize',20)
grid on
grid minor
axis fill
pbaspect([2 1 1])
xlim([0 31])


f_traj = xlist_full;
N_traj = ylist_full;

save('abs_pop_traj.mat','nk','ns','f_traj','N_traj');

%% two paths together
N=20;
A = ones(n+1);
A(1,1) = 0;
% lamb = 1;
rho = 0.65;
%rho = 0.001015;
num_dilution = 3
dt = h;
choice = 'conservative'
figure
contourf(X,Y,A,[0,1],'EdgeColor','none');
xx = 0:h:1;
[X,Y] = meshgrid(xx,xx);

f_lim = @(r,rho,T) 1./(rho+(1-rho)./(1-exp(-(r.*T+log(rho)))));
mypink = [148 0 211]/255;
mygreen  =[46 139 87]/255;
mygreen2  =[10 110 90]/255*1.2;
for ii = 2:N
        for ij = [2:1:N+1]
            if isnan(u(ii,ij)) || isnan(v(ii,ij))
                continue
            else
                ah = annotation('arrow',...
                    'Color', [169/255,169/255,169/255],...
                    'headStyle','cback1','HeadLength',3.5,'HeadWidth',3.5);
                set(ah,'parent',gca);
                if A_down(ii,ij) == 0
                    set(ah,'position',[x(ii,ij) y(ii,ij) 0.03*u(ii,ij) 0.03*v(ii,ij)]);
                else
                    set(ah,'position',[x(ii,ij) y(ii,ij) 0.03*uu(ii,ij) 0.03*vv(ii,ij)]);
                end
            end
        end
    end
hold on
for jj = 1:2
    if jj == 1
        xc = 0.5;
        yc = 0.7;
    else
        xc = 0.5;
        yc = 0.1;
    end
    [xlist_full,ylist_full,policy_list,xstart, ystart, xend, yend] ...
        = Optimal_path_DilutionIndexed(num_dilution,h,A,lamb,rks,gamma,epsilon,dt,xc,yc,choice,rho);

    
    if jj == 1
        mycolor = mypink;
        load('kwin.mat');
    else
        mycolor = mygreen;
        load('klose.mat');
    end
    plot(xlist_full(policy_list==1),ylist_full(policy_list==1),'.','markersize',3,'Color',mycolor);
    plot(xlist_full(policy_list==0),ylist_full(policy_list==0),'g.','markersize',3);
    plot(xstart(1),ystart(1),'o','markersize',4.5,'linewidth',1.5,'markerfacecolor','w','MarkerEdgeColor','w');
    plot(xend(1),yend(1),'s','markersize',4,'linewidth',1.5,'markerfacecolor','k','MarkerEdgeColor','k');
    %
    %
    plot(xstart(2),ystart(2),'o','markersize',4,'linewidth',1.5,'markerfacecolor','c','MarkerEdgeColor','c');
    plot(xend(2),yend(2),'s','markersize',4,'linewidth',1.5,'markerfacecolor','k','MarkerEdgeColor','k');
    %
    %
    plot(xstart(3),ystart(3),'o','markersize',4,'linewidth',1.5,'markerfacecolor','c','MarkerEdgeColor','c');
    plot(xend(3),yend(3),'s','markersize',4,'linewidth',1.5,'markerfacecolor','k','MarkerEdgeColor','k');
    plot(xstart(4),ystart(4),'o','markersize',4,'linewidth',1.5,'markerfacecolor','c','MarkerEdgeColor','c');
    if jj == 1
        % Regular kilers win (0.5,0.7)
        ylast_end = f_lim(rks*(1-epsilon),rho,1);
        ylast_init = ylast_end*rho;
        ylast = ylast_init:0.001:ylast_end;
        xlast = 0.998*ones(1,length(ylast));

        plot(xlast,ylast,'.-','markersize',4,'linewidth',3.5,'Color',mypink);
        plot(xlast(1),ylast_init,'o','markersize',4.9,'linewidth',1.5,'markerfacecolor','c','MarkerEdgeColor','c');
        plot(xlast(1),ylast_end,'s','markersize',4.9,'linewidth',1.5,'markerfacecolor','k','MarkerEdgeColor','k');

        annotation('arrow',[0.563,0.563],[0.705,0.515],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',7);
        annotation('arrow',[0.588,0.588],[0.595,0.445],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',7);
        annotation('arrow',[0.606,0.606],[0.555,0.415],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',7);
        annotation('arrow',[0.821,0.821],[0.44,0.45],'color','k','linewidth',1.3,'linestyle','--','headwidth',13,'headlength',6);
        annotation('arrow',[0.821,0.821],[0.39,0.38],'color','k','linewidth',1.3,'linestyle','--','headwidth',13,'headlength',6);
        text(0.685,0.38,'The process','fontsize',12,'color',mypink);
        text(0.685,0.33,'continues...','fontsize',12,'color',mypink);
        text(xstart(1)-0.03,ystart(1)-0.03,'Init','fontsize',12);
        text(xstart(2)-0.02,ystart(2)-0.04,'1','fontsize',12);
        text(xstart(3)-0.02,ystart(3)-0.04,'2','fontsize',12);
        text(xstart(4)-0.02,ystart(4)-0.04,'3','fontsize',12);
    else

        % Regular kilers lose (0.5,0.1)
        ylast_end = f_lim(1,rho,1);
        ylast_init = ylast_end*rho;
        ylast = ylast_init:0.001:ylast_end;
        xlast = 0.001*ones(1,length(ylast));

        plot(xlast,ylast,'.-','markersize',4,'linewidth',3.5,'Color',mygreen);
        plot(xlast(1),ylast_init,'o','markersize',4.9,'linewidth',1.5,'markerfacecolor','c','MarkerEdgeColor','c');
        plot(xlast(1),ylast_end,'s','markersize',4.9,'linewidth',1.5,'markerfacecolor','k','MarkerEdgeColor','k');

        annotation('arrow',[xend(1)+0.039,xend(1)+0.039],[yend(1)+0.08,ystart(2)+0.095],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',7);
        annotation('arrow',[xend(2)+0.055,xend(2)+0.055],[yend(2)+0.07,ystart(3)+0.09],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',7);
        annotation('arrow',[xend(3)+0.071,xend(3)+0.071],[yend(3)+0.05,ystart(4)+0.085],'linewidth',1.3,'linestyle','--','headwidth',8,'headlength',7);
        annotation('arrow',[0.213,0.213],[0.63,0.64],'color','k','linewidth',1.3,'linestyle','--','headwidth',13,'headlength',6);
        annotation('arrow',[0.213,0.213],[0.515,0.505],'color','k','linewidth',1.3,'linestyle','--','headwidth',13,'headlength',6);
        text(0.06,0.4,'The process','fontsize',12,'color',mygreen2);
        text(0.06,0.35,'continues...','fontsize',12,'color',mygreen2);
        text(xstart(1)-0.03,ystart(1)-0.03,'Init','fontsize',12);
        text(xstart(2)-0.02,ystart(2)-0.04,'1','fontsize',12);
        text(xstart(3)-0.02,ystart(3)-0.04,'2','fontsize',12);
        text(xstart(4)-0.02,ystart(4)-0.04,'3','fontsize',12);
    end




end
xlabel('Fraction of the killer (f)','fontsize',14);
ylabel('Total population (N)','fontsize',14);
axis equal

xlim([0,1]);
ylim([0,1]);
colormap("copper")