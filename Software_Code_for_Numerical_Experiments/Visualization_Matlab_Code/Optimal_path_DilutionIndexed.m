function [xlist_full,ylist_full,policy_list,xstart, ystart, xend, yend] ...
    = Optimal_path_DilutionIndexed(num_dilution,h,policy_mat,lamb,r_ks,gamma,epsilon,dt,xc,yc,choice,rho)
% This function computes the optimal path for a given number of dilutions
% Inputs:
%   num_dilution: number of dilutions
%   h: step size for the grid
%   policy_mat: policy matrix containing the policies for each grid point
%   lamb: growth rate of the population
%   r_ks: rescaled ratio of the growth rate of the killer to that of the sensitive
%   gamma: rescaled cost of toxin-production rate
%   epsilon: cost of constitutive toxin production
%   dt: time step for the simulation
%   xc: initial fraction of the killer
%   yc: initial total population size
%   choice: policy determination choice ('conservative', 'aggressive', or 'majority')
%   rho: rescaled survival rate
% Outputs:
%   xlist_full: list of x-coordinates (fraction of the killer) over time
%   ylist_full: list of y-coordinates (total population size) over time
%   policy_list: list of policies chosen at each time step
%   xstart: starting position before each dilution (including the IC) for the x-coordinates
%   ystart: starting points before each dilution (including the IC) for the y-coordinates
%   xend: ending points after each dilution for the x-coordinates
%   yend: ending points after each dilution for the y-coordinates
%
% Author: MingYi Wang, Cornell University
% Last Modified: 06/2025
%

%% Initialization
% Drift part of the population growth dynamics (Eq.[2] in the paper with "delta == 0")
fx = @(x,y,a)  x.*(1-x).*((1-y).*(r_ks.*(1-epsilon.*a) - 1) + a.*gamma.*x.*y);
fy = @(x,y,a)  y.*(1-y).*(1 + (r_ks.*(1-epsilon.*a) - 1).*x) - a.*gamma.*y.^2.*x.*(1-x);

xx = 0:h:1; % array of sample gridpts
N = length(policy_mat) - 1; % number of gridpts
T_list = ones(500); %List of inter-dilution times; default is 1 for each dilution

T0_list = zeros(1,num_dilution); % List of initial times of each dilution
Tf_list = zeros(1,num_dilution); % List of final times of each dilution
T0_list(1) = 0; % Initial time is 0
Tf_list(1) = T_list(1); % Initialize the first final time

for n = 2:num_dilution
    T0_list(n) = T0_list(n-1) + T_list(n-1); % Initial time of the nth dilution
    Tf_list(n) = Tf_list(n-1) + T_list(n); % Final time of the nth dilution
end

% Initialize lists to store results
ylist_full=[];
xlist_full=[];
policy_list = [];
xstart = [xc];
ystart = [yc];
xend = [];
yend = [];

% loop over all dilutions
for ii = 1:num_dilution 
    T0 = T0_list(ii);
    Tf = Tf_list(ii);
    tt = T0:dt:Tf;
    tnum = length(tt);%number of time steps
    ylist = zeros(1,tnum);
    xlist = zeros(1,tnum);
    
    if ii == 1
        xlist(1) = xc;
        ylist(1) = yc;
    else
        xlist(1) = xlist_full(end);
        ylist(1) = rho*ylist_full(end); %dilution occurs
        xstart = [xstart, xlist(1)];
        ystart = [ystart, ylist(1)];
 
    end
    
    for nn = 1:tnum-1 % loop over all time steps within each dilution
        yval= ylist(nn);
        xval= xlist(nn);
        if yval > 0
            kx=find(xval<= xx',1);
            ky=find(yval<= xx',1);
            
            Edge_flag = If_At_Edge(kx, ky, N); % check if at edge
            
            if Edge_flag
                % If at the edge, find the nearest neighbor in the grid
                [nearest_x, nearest_y] = nearest_neighbor(kx,ky,h,xval,yval);
                % Get the policy at the nearest neighbor
                [policy] = Return_Policy_At_Edge(policy_mat,nearest_x, nearest_y);
            else
                % If not at the edge, use the 4-point stencil to determine the policy
                d1=policy_mat(ky-1,kx-1);d2=policy_mat(ky-1,kx);
                d3=policy_mat(ky,kx-1);d4=policy_mat(ky,kx);
                d_square=[d1,d2,d3,d4];
                d_square(isnan(d_square)) = 0;
                [policy] = Return_Policy(d_square, choice);
            end
            policy_list=[policy_list,policy];
            
            %% Compute the drift at the current point
            dxx = fx(xval,yval,policy);
            dyy = fy(xval,yval,policy);
            
            % Dealing with numerical issues
            if abs(dxx) < 1e-15
                dxx = 0;
            elseif abs(dyy) < 1e-15
                dyy = 0;
            end
            
            ylist(nn+1) = min(yval+dyy*dt,1);
            xlist(nn+1) = min(xval+dxx*dt,1);
            if nn == tnum - 1
                policy_list=[policy_list,nan];
            end
            
            if (xlist(nn+1) > 0.9) && (abs(xlist(nn+1) - xlist(nn)) < 0.5*dt^2)
                xlist(nn+1) = xlist(nn);
            end
            
            
        end
    end
    % Store the results for the current dilution
    xend = [xend, xlist(end)];
    yend = [yend, ylist(end)];
    xlist_full = [xlist_full,xlist];
    ylist_full = [ylist_full,ylist];
end


%% subfucntions
function Edge_flag = If_At_Edge(kx, ky, N)
    % This function checks if the current grid point is at the edge of the grid
    if (ky == N + 1) && (kx == N + 1)
        Edge_flag = true;
    elseif (ky == N + 1)
        Edge_flag = true;
    elseif (kx == N + 1)
        Edge_flag = true;
    else
        Edge_flag = false;
    end
end

function [policy] = Return_Policy(d_square, choice)
    % This function determines the policy based on the values in d_square
    % d_square: a 4-pt stencil containing the values of the policy at the neighboring grid points
    switch choice
        case 'conservative'
            if (sum(d_square) == 4)
                policy = 1;
            else
                policy = 0;
            end
            
        case 'aggressive'
            if (sum(d_square) > 0)
                policy = 1;
            else
                policy = 0;
            end
            
        case 'majority'
            if (sum(d_square) >= 3)
                policy = 1;
            else
                policy = 0;
            end
            
    end
end

function Indx = nearest_index(right_indx, xloc, dx)
    % This function finds the nearest index to xloc in the grid defined by dx
    % right_indx: the index of the grid point that is greater than or equal to xloc
    % xloc: the location for which we want to find the nearest index
    % dx: the step size of the grid

    if ((right_indx - 1)*dx - xloc) <= 0.5*dx
        Indx = right_indx;
    else
        Indx = right_indx - 1;
    end
end

function [nearest_x, nearest_y] = nearest_neighbor(idx_x,idx_y,h,xloc,yloc)
    % This function finds the nearest grid point to (xloc, yloc) in the grid defined by h
    % idx_x: index of the grid point in the x-direction
    % idx_y: index of the grid point in the y-direction
    % h: step size of the grid
    % xloc: x-coordinate of the location for which we want to find the nearest grid point
    % yloc: y-coordinate of the location for which we want to find the nearest grid point
    % Outputs:
    % nearest_x: index of the nearest grid point in the x-direction
    % nearest_y: index of the nearest grid point in the y-direction
    nearest_x = nearest_index(idx_x, xloc, h);
    nearest_y = nearest_index(idx_y, yloc, h);
end

function [policy] = Return_Policy_At_Edge(Policy_mat,nearest_x, nearest_y)
    % This function returns the policy at the nearest grid point
    % Policy_mat: matrix containing the policies for each grid point
    % nearest_x: index of the nearest grid point in the x-direction
    % nearest_y: index of the nearest grid point in the y-direction
    % Outputs:
    % policy: policy at the nearest grid point

    if nearest_x == 1 && nearest_y == 1
        nearest_x = 2;
        nearest_y = 2;
    elseif nearest_x == 1
        nearest_x = 2;
    elseif nearest_y == 1
        nearest_y = 2;
    end

    policy = Policy_mat(nearest_y, nearest_x);

end 

end