function [f1,z1] = plotting_expT_v2(Xf_list, Nf_list, nbins, binwidth,xc, yc)
[f1,z1] = ksdensity(Xf_list,'Kernel','epanechnikov');
hh = figure;
histogram(Xf_list,'binwidth',binwidth,'Normalization','pdf');
% hold on;
% plot(z1,f1,'-.r','linewidth',1.3);
% hold off
% xlabel('Fraction of killer cells (f)');
% title11 = sprintf('Empirical distribution of f(T) where T~Exp(%g) \n starting from (x,y) = (%g,%g) with \x03b3=%g',lamb,xc,yc,gamma);
% title(title11);
xlim([-0.1,1.1]);
ylim([0,21]);
grid on;
axis square
pause(0.01);
filename = sprintf('rks1/binomial_pdf_unif_rks1_xc%gyc%g.png',xc*1000,yc*1000);
%     saveas(gcf,filename);
print(filename,'-dpng');
close(hh);
pause(0.01);
% 
% 
% [f2,z2] = ksdensity(Nf_list,'Support',[0 1],'Kernel','epanechnikov');
% hh = figure;
% histogram(Nf_list,nbins,'Normalization','pdf');
% hold on;
% plot(z2,f2,'-.r','linewidth',1.3)
% hold off;
% xlabel('Total population (N)');
% title11 = sprintf('Empirical distribution of N(T) where T~Exp(%g) \n starting from (x,y) = (%g,%g) with \x03b3=%g',lamb,xc,yc,gamma);
% title(title11);
% grid on;
% pause(0.01);
% filename = sprintf('Population rate=%g gamma=%g.png',lamb,gamma);
% %     saveas(gcf,filename);
% print(filename,'-dpng');
% close(hh);
% pause(0.01);

end