function Visualization_Fig3B()
% Visualization_Fig3B.m
% This function generates the Fig. 3B in the paper, which visualizes
% the separatrix of the population dynamics of the system (Eq. [2] in the main text)  
% under different death rates.
%
% Author: MingYi Wang, Cornell University
% Last Modified: 06/2025
%
clear; clc;
%% Initialization
NN = 800; % Number of grid points we used to compute the direction field
death_arr = [0.05,0.215,0.3,0.43,0.5];
yy = linspace(0,1,NN+1);
xx = linspace(0,1,NN+1);
[X,Y] = meshgrid(xx,yy);

% Background matrix for visualization
Zmat = ones(NN+1);
Zmat(1,1) = 0; % For visualization purpose, we set the first element to 0
%% Plotting
figure
contourf(xx,yy,Zmat);
hold on
for ii = 1:length(death_arr)
    str = sprintf('separatrix_%g.mat',ii);
    load(str,'y_unstable_backward',"y_stable_backward",'y_unstable_forward',"y_stable_forward");
    nn1 = y_stable_forward(:,1) + y_stable_forward(:,2);
    ff1 = y_stable_forward(:,1)./nn1;
    ff1(ff1>1) = 1;
    nn2 = y_stable_backward(:,1) + y_stable_backward(:,2);
    ff2 = y_stable_backward(:,1)./nn2;
    ff2(ff2>1) = 1;
    indx = find(nn2 < 0,1) - 1;
    ff2(ff2<1e-3) = 0;
    nn2(nn2<1e-3) = 0;
    plot(ff1,nn1,'-.','LineWidth',2.7,'Color','m');
    plot(ff2(1:indx),nn2(1:indx),'-.','LineWidth',2.7,'Color','m');

    hold on
end
for jj = 1:length(death_arr)
    str = sprintf('separatrix_%g.mat',jj);
    load(str,'y_unstable_backward',"y_stable_backward",'y_unstable_forward',"y_stable_forward");
    nn1 = y_stable_forward(:,1) + y_stable_forward(:,2);
    ff1 = y_stable_forward(:,1)./nn1;
    ff1(ff1>1) = 1;
    strr = sprintf('$\\delta=%g$',death_arr(jj));
    if jj == 1
        text(ff1(1)-0.01,0.025,strr,'FontSize',11,'Interpreter','latex','FontWeight','bold');
    elseif jj == 2
        text(ff1(2)-0.1,nn1(1),strr,'FontSize',11,'Interpreter','latex','FontWeight','bold');
    else

        text(ff1(1)-0.06,nn1(1),strr,'FontSize',11,'Interpreter','latex','FontWeight','bold');
    end
end
hold off
axis equal
xlabel('Fraction of the killer (f)','FontSize',14);
ylabel('Total population (N)','FontSize',14);
set(findall(gca,'type','text'),'visible','on');
colormap("copper")

end
