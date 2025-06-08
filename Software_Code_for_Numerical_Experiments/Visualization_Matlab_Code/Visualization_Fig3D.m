function Visualization_Fig3D()
% Visualization_Fig3D.m
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
rho = 0.65; % rescaled survival rate
%% two paths together
% Compute the vector field
N=20; % Number of grid points we intended to plot the vector field
NN = 800; % Number of grid points we used to compute the vector field
h = 1/ NN; % Step size for the discretization
dt = h; % Time step for the ODE solver
choice = 'conservative' % choice of determining wheter to produce toxin or not (in this case, the killer always produces toxin)
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
        num_dilution = 4; % Number of dilutions
        [xlist_full,ylist_full,policy_list,xstart, ystart, xend, yend] ...
            = Optimal_path_DilutionIndexed(num_dilution,h,A,lamb,rks,gamma,epsilon,dt,xc,yc,choice,rho);
        save('kwin.mat','xstart','ystart','xend','yend');
    else
        xc = 0.5;
        yc = 0.1;
        num_dilution = 4; % Number of dilutions
        [xlist_full,ylist_full,policy_list,xstart, ystart, xend, yend] ...
            = Optimal_path_DilutionIndexed(num_dilution,h,A,lamb,rks,gamma,epsilon,dt,xc,yc,choice,rho);
        save('klose.mat','xstart','ystart','xend','yend');   
    end

    num_dilution = 3; % Number of dilutions
    [xlist_full,ylist_full,policy_list,xstart, ystart, xend, yend] ...
        = Optimal_path_DilutionIndexed(num_dilution,h,A,lamb,rks,gamma,epsilon,dt,xc,yc,choice,rho);
    
    if jj == 1
        mycolor = mypink;
        load('kwin.mat','xstart','ystart','xend','yend');
    else
        mycolor = mygreen;
        load('klose.mat','xstart','ystart','xend','yend');   
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

end