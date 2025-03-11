function [xlist_full,ylist_full,xstart, ystart, xend, yend] ...
    = Const_path_DilutionBinomial(Tp,num_dilution,h,r_ks,gamma,epsilon,dt,xc,yc,rho,Kc)

% velocity function for different r_k and r_s
fx = @(x,y,a)  x.*(1-x).*((1-y).*(r_ks.*(1-epsilon.*a) - 1) + a.*gamma.*x.*y);
fy = @(x,y,a)  y.*(1-y).*(1 + (r_ks.*(1-epsilon.*a) - 1).*x) - a.*gamma.*y.^2.*x.*(1-x);

xx = 0:h:1;
% N = length(A) - 1;

T_list = Tp*ones(num_dilution);
T0_list(1) = 0;
Tf_list(1) = T_list(1);

for n = 2:num_dilution
    T0_list(n) = T0_list(n-1) + T_list(n-1);
    Tf_list(n) = Tf_list(n-1) + T_list(n);
end
total_time = Tf_list(end);

% break_flag = false;
% Euler
ylist_full=[];
xlist_full=[];
% policy_list = [];
xstart = [xc];
ystart = [yc];
xend = [];
yend = [];
policy = 1;

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
        [frac,pop] = sampling_binomial(xlist_full(end),ylist_full(end),Kc,rho);
        xlist(1) = frac;
        ylist(1) = pop;
        xstart = [xstart, xlist(1)];
        ystart = [ystart, ylist(1)];
 
    end
    
    for nn = 1:tnum-1
        
        yval= ylist(nn);
        xval= xlist(nn);
        if yval > 0
            kx=find(xval<= xx',1);
            ky=find(yval<= xx',1);
            
            
            
            dxx = fx(xval,yval,policy);
            dyy = fy(xval,yval,policy);
            
            if abs(dxx) < 1e-15
                dxx = 0;
            elseif abs(dyy) < 1e-15
                dyy = 0;
            end
            
            ylist(nn+1) = min(yval+dyy*dt,1);
            xlist(nn+1) = min(xval+dxx*dt,1);
%             if nn == tnum - 1
%                 policy_list=[policy_list,nan];
%             end
            
            if (xlist(nn+1) > 0.9) && (abs(xlist(nn+1) - xlist(nn)) < 0.5*dt^2)
                xlist(nn+1) = xlist(nn);
            end
            
            
        end
    end
    xend = [xend, xlist(end)];
    yend = [yend, ylist(end)];
    xlist_full = [xlist_full,xlist];
    ylist_full = [ylist_full,ylist];
    
    
end
% xlist_full = [xlist_full,xlist(end)];
% ylist_full = [ylist_full,rho*ylist(end)];
% policy_list = [policy_list, nan];
end

%% subfucntions
function Edge_flag = If_At_Edge(kx, ky, N)
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

function [policy] = ...
    Return_Policy(d_square, choice)

switch choice
    case 'conservative'
        if (sum(d_square) == 4)
            policy = 1;
            
            %             a1_set=[a1_set,nn];
            
        else
            policy = 0;
            
            %             a0_set=[a0_set,nn];
            
        end
        
    case 'aggressive'
        if (sum(d_square) > 0)
            policy = 1;
            
            %             a1_set=[a1_set,nn];
            
        else
            policy = 0;
            
            %             a0_set=[a0_set,nn];
            
        end
        
    case 'majority'
        if (sum(d_square) >= 3)
            policy = 1;
            
            %             a1_set=[a1_set,nn];
            
        else
            policy = 0;
            
            %             a0_set=[a0_set,nn];
            
        end
        
end
end

function Indx = nearest_index(right_indx, xloc, dx)
if ((right_indx - 1)*dx - xloc) <= 0.5*dx
    Indx = right_indx;
else
    Indx = right_indx - 1;
end
end

function [nearest_x, nearest_y] = nearest_neighbor(idx_x,idx_y,h,xloc,yloc)
nearest_x = nearest_index(idx_x, xloc, h);
nearest_y = nearest_index(idx_y, yloc, h);
end

function [policy] = ...
    Return_Policy_At_Edge(Policy_mat,nearest_x, nearest_y)

if nearest_x == 1 && nearest_y == 1
    nearest_x = 2;
    nearest_y = 2;
elseif nearest_x == 1
    nearest_x = 2;
elseif nearest_y == 1
    nearest_y = 2;
end

policy = Policy_mat(nearest_y, nearest_x);


% if (policy == 1)
%     a1_set=[a1_set,nn];
% else
%     a0_set=[a0_set,nn];
%
% end

end