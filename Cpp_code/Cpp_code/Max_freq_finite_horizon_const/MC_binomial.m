%% MC simulations -- Binomial sampling
clear;
clc;
%% Initialization
n = 800; % number of mesh points along on side
h = 1/n;
dt = h;

% epsilon = 0.2;
% rks = 0.85;

epsilon = 0.1;
rks = 1.0;

lamb = 1.0;
gamma = 1.0;
rho = 0.65;

thres_fail = 0.01;
thres_success = 1 - thres_fail;
tol = 1e-6;

Kc = 1e5;
Num_Samples = 1e5;
num_dilution = 200;
% num_dilution = 50;
Tp = 1;


% velocity function for different r_k and r_s
fx = @(x,y,a)  x.*(1-x).*((1-y).*(rks.*(1-epsilon.*a) - 1) + a.*gamma.*x.*y);

fy = @(x,y,a)  y.*(1-y).*(1 + (rks.*(1-epsilon.*a) - 1).*x) - a.*gamma.*y.^2.*x.*(1-x);

% cdf of exponential distribution
ft = @(t,dt,lamb) exp(-lamb.*t) - exp(-lamb.*(t+dt));
% pdf of exponential distribution
ftt = @(t,lamb) lamb.*exp(-lamb.*t);

xc_list = linspace(0.1,0.9,9);
yc_list = linspace(0.1,0.9,9);
xx = linspace(0,1,n+1);
[X,Y] = meshgrid(xx,xx);

%% MC
xc_list = linspace(0.05,0.15,5);
diff_x = (xc_list(2) - xc_list(1))/2
diff_y = (yc_list(2) - yc_list(1))/2
%%
% tic
% pdf_f_list = cell(length(xc_list),length(yc_list));
% pdf_N_list = cell(length(xc_list),length(yc_list));
% 
% for i = 1:length(xc_list)
%     xc = xc_list(i)
%     for j = 1:length(yc_list)
%         yc = yc_list(j);
%         pdf_f = zeros(1,Num_Samples);
%         pdf_N = zeros(1,Num_Samples);
%         parfor ii = 1:Num_Samples
%             
%             [xlist_full,ylist_full,xstart, ystart, xend, yend] ...
%                 = Const_path_DilutionBinomial(Tp,num_dilution,h,rks,gamma,epsilon,dt,xc,yc,rho,Kc);
% 
%             pdf_f(ii) = xlist_full(end);
%             pdf_N(ii) = ylist_full(end);
% 
%             ii
%         end
%         pdf_f_list{i,j} = pdf_f;
%         pdf_N_list{i,j} = pdf_N;
%     end
%     
% end
% t_total = toc


tic
pdf_f_list_unif = cell(length(xc_list),length(yc_list));
pdf_N_list_unif = cell(length(xc_list),length(yc_list));

for i = 1:length(xc_list)
    xc = xc_list(i);
    for j = 1:length(yc_list)
        yc = yc_list(j);
        pdf_f = zeros(1,Num_Samples);
        pdf_N = zeros(1,Num_Samples);
        parfor ii = 1:Num_Samples

%             xc_low = xc - 0.05;
%             xc_up = xc + 0.05;
%             yc_low = yc - 0.05;
%             yc_up = yc + 0.05;
            xc_low = xc - diff_x;
            xc_up = xc + diff_x;
            yc_low = yc - diff_y;
            yc_up = yc + diff_y;

            xc_sample = xc_low + (xc_up - xc_low)*rand;
            yc_sample = yc_low + (yc_up - yc_low)*rand;
            
            [xlist_full,ylist_full,xstart, ystart, xend, yend] ...
                = Const_path_DilutionBinomial(Tp,num_dilution,h,rks,gamma,epsilon,dt,xc_sample,yc_sample,rho,Kc);

            pdf_f(ii) = xlist_full(end);
            pdf_N(ii) = ylist_full(end);

            ii
        end
        pdf_f_list_unif{i,j} = pdf_f;
        pdf_N_list_unif{i,j} = pdf_N;
    end
    
end
t_total = toc

%% plotting
mean_table = zeros(length(yc_list),length(xc_list));
med_table = zeros(length(yc_list),length(xc_list));
[X,Y] = meshgrid(xc_list,yc_list);
nbins = 100;
binwidth = 0.05;
for i = 1:length(xc_list)
    xc = xc_list(i);
    for j = 1:length(yc_list)
        yc = yc_list(j);
% for i = 1
%     xc = xc_list(i);
%     for j = 9
%         yc = yc_list(j);

        pdf_f =  pdf_f_list{i,j} ; 
        pdf_N = pdf_N_list{i,j} ;
        mean_table(j,i) = mean(pdf_f);
        med_table(j,i) = median(pdf_f);
%         [f1,z1] = plotting_expT_v2(pdf_f,pdf_N, nbins,binwidth, xc, yc);


%         pdf_f =  pdf_f_list_unif{i,j} ; 
%         pdf_N = pdf_N_list_unif{i,j} ;
%         mean_table(j,i) = mean(pdf_f);
%         med_table(j,i) = median(pdf_f);
%         [f1,z1] = plotting_expT_v2(pdf_f,pdf_N, nbins,binwidth, xc, yc);



    end
end

%%

figure
imagesc(xc_list,yc_list,mean_table);

% hold on
% contour(X,Y,f_slice2,'k--','linewidth',2);
% hold off

xlabel('Initial fraction of the killer (x = f_0)','FontSize',14);
ylabel('Initial total population (y = N_0)','FontSize',14);
% xlabel('Mean initial fraction of the killer (x = f_0)','FontSize',14);
% ylabel('Mean initial total population (y = N_0)','FontSize',14);
xticks(xc_list);
yticks(yc_list);
colorbar();
colormap("jet");
clim([0 1]);
% axis square
pbaspect([1 2 1])
set(gca,'YDir','normal');


% figure
% imagesc(xc_list,yc_list,med_table);
% % hold on
% % contour(X,Y,f_slice2,'k--','linewidth',2);
% % 
% % hold off
% % xlabel('Initial fraction of the killer (x)','FontSize',14);
% % ylabel('Initial total population (y)','FontSize',14);
% xlabel('Mean initial fraction of the killer (x)','FontSize',14);
% ylabel('Mean initial total population (y)','FontSize',14);
% xticks(xc_list);
% yticks(yc_list);
% colorbar();
% colormap("jet");
% clim([0 1]);
% % axis square
% pbaspect([1 2 1])
% set(gca,'YDir','normal');

%%
xc = 0.1;
yc = 0.2;

% xc_low = xc - 0.05;
% xc_up = xc + 0.05;
% yc_low = yc - 0.05;
% yc_up = yc + 0.05;
% 
% xc_sample = xc_low + (xc_up - xc_low)*rand;
% yc_sample = yc_low + (yc_up - yc_low)*rand;
xc_sample = xc;
yc_sample = yc;

% for ii = 1:100
[xlist_full,ylist_full,xstart, ystart, xend, yend] ...
                = Const_path_DilutionBinomial(Tp,num_dilution,h,rks,gamma,epsilon,dt,xc_sample,yc_sample,rho,Kc);
% if xlist_full(end) > 0.9
%     xc_sample
%     yc_sample
%     break
% end
% ii
% end

%%
figure
nk_win = xlist_full.*ylist_full;
ns_win = (1-xlist_full).*ylist_full;
nn = length(xlist_full);
t_list = dt.*(1:nn);
plot(t_list,nk_win,'r-','linewidth',1.5);
hold on 
plot(t_list,ns_win,'b-','linewidth',1.5);
hold off
xlabel("Rescaled time (t)",'fontsize',18);
ylabel("Normalized population",'fontsize',18);

yyaxis right
plot(t_list,xlist_full,'-.','linewidth',2.5);
ylabel("Fraction of the killer",'fontsize',18);


legend('Constitutive killer','Sensitive','Frac of the killer (f)','fontsize',15,'location','best');

title('Killers win','fontsize',20)
grid on
grid minor

f_win = xlist_full;
N_win = ylist_full;

% save('killer win.mat','nk_win','ns_win','f_win','N_win');
% save('killer win unif.mat','nk_win','ns_win','f_win','N_win');
save('killer win rks1.mat','nk_win','ns_win','f_win','N_win');
%%
figure
nk_lose = xlist_full.*ylist_full;
ns_lose = (1-xlist_full).*ylist_full;
nn = length(xlist_full);
t_list = dt.*(1:nn);
plot(t_list,nk_lose,'r-','linewidth',1.5);
hold on 
plot(t_list,ns_lose,'b-','linewidth',1.5);
hold off
xlabel("Rescaled time (t)",'fontsize',18);
ylabel("Normalized population",'fontsize',18);


yyaxis right
plot(t_list,xlist_full,'-.','linewidth',2.5);
ylabel("Fraction of the killer",'fontsize',18);

legend('Constitutive killer','Sensitive','Frac of the killer (f)','fontsize',15,'location','best');
title('Sensitives win','fontsize',15)
grid on
grid minor


f_lose = xlist_full;
N_lose = ylist_full;

% save('killer lose.mat','nk_lose','ns_lose','f_lose','N_lose');
% save('killer lose unif.mat','nk_lose','ns_lose','f_lose','N_lose');
save('killer lose rks1.mat','nk_lose','ns_lose','f_lose','N_lose');
%%
mygreen = [50,205,50]/255;
nn = length(xlist_full);
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
    'fontsize',12,'location','southoutside','interpreter','latex');
title('Killers win','fontsize',15)
grid on
% grid minor
xlim([0,205]);
% ylim([0.49,1.01]);
ylim([0.05,1.03]);
% title('Killers lead','fontsize',15)
% grid on
% % grid minor
% xlim([0,52]);
% ylim([0.49,0.66]);
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
xlim([0,205]);
% ylim([0,0.6]);
% title('Sensitives lead','fontsize',15)
% grid on
% % grid minor
% xlim([0,52]);
% ylim([0.38,0.55]);
ax = gca;
ax.YAxis(2).Color = mygreen;

ttt = text(-20,0.03,'Normalized population','fontsize',18);
ttt2 = text(225,0.06,'Fraction of the killer','fontsize',18,'Color',mygreen);
ttt.Rotation = 90;
ttt2.Rotation = 90;
