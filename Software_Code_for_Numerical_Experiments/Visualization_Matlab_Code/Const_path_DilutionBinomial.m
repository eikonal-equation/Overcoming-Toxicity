function [xlist_full,ylist_full,xstart, ystart, xend, yend] ...
    = Const_path_DilutionBinomial(Tp,num_dilution,h,r_ks,gamma,epsilon,dt,xc,yc,rho,Kc)
    % This function simulates the population growth dynamics under periodic dilutions, 
    % using the "Binomial sampling method (described in Section S9 of the SI Appendix)".
    % Here, we assume constitutive toxin production rate for the killer strain.
    %
    % Input:
    %   Tp: time period of each dilution
    %   num_dilution: number of dilutions
    %   h: grid size for the population fractions
    %   r_ks: rescaled growth rate of the killer strain
    %   gamma: rescaled cost of toxin production rate
    %   epsilon: cost of constitutive toxin production
    %   dt: time step size
    %   xc: initial fraction of the killer strain
    %   yc: initial normalized the total population
    %   rho: rescaled survival rate
    %   Kc: shared carrying capacity 
    % Output:
    %   xlist_full: full list of the fraction of the killer strain over time
    %   ylist_full: full list of the normalized the total population over time
    %   xstart: starting fraction of the killer strain
    %   ystart: normalized the total population
    %   xend: ending fraction of the killer strain
    %   yend: ending normalized the total population
    %
    % Author: MingYi Wang, Cornell University
    % Last Modified: 06/2025
    %

%% Initialization
% Drift part of the population growth dynamics (Eq.[2] in the paper with "delta == 0")
fx = @(x,y,a)  x.*(1-x).*((1-y).*(rks.*(1-epsilon.*a) - 1) + a.*gamma.*x.*y);
fy = @(x,y,a)  y.*(1-y).*(1 + (rks.*(1-epsilon.*a) - 1).*x) - a.*gamma.*y.^2.*x.*(1-x);

xx = 0:h:1; % sample points in x-direction

T_list = Tp*ones(num_dilution); % sequence of inter-dilution times
T0_list(1) = 0; % starting time of the first dilution
Tf_list(1) = T_list(1); % ending time of the first dilution

% Calculate the total time for each dilution
for n = 2:num_dilution
    T0_list(n) = T0_list(n-1) + T_list(n-1);
    Tf_list(n) = Tf_list(n-1) + T_list(n);
end

% Initialize the lists for storing the population dynamics
ylist_full=[];
xlist_full=[];
xstart = [xc];
ystart = [yc];
xend = [];
yend = [];
policy = 1; % constitutive toxin production policy

for ii = 1:num_dilution % loop over all dilutions
    T0 = T0_list(ii); % get starting time of the current dilution
    Tf = Tf_list(ii); % get ending time of the current dilution
    tt = T0:dt:Tf; % time vector for the current dilution
    tnum = length(tt);% number of time steps
    ylist = zeros(1,tnum); % initialize the normalized total population for the current dilution
    xlist = zeros(1,tnum); % initialize the fraction of the killer strain for the current dilution
    
    if ii == 1 % For the first dilution, use the initial conditions
        xlist(1) = xc;
        ylist(1) = yc;
    else
        % For subsequent dilutions, use "Binomial sampling method" to obtain
        % post-dilution populations to initialize the next round
        [frac,pop] = sampling_binomial(xlist_full(end),ylist_full(end),Kc,rho);
        xlist(1) = frac;
        ylist(1) = pop;
        xstart = [xstart, xlist(1)];
        ystart = [ystart, ylist(1)];
 
    end
    
    for nn = 1:tnum-1 % loop over all time steps in the current dilution
        yval= ylist(nn); % current normalized total population
        xval= xlist(nn); % current fraction of the killer strain
        if yval > 0

            % Compute the drift part of the population growth dynamics
            dxx = fx(xval,yval,policy); 
            dyy = fy(xval,yval,policy);
            
            % Dealing with numerical issues
            if abs(dxx) < 1e-15
                dxx = 0;
            elseif abs(dyy) < 1e-15
                dyy = 0;
            end
            % Update the dynamics using a first-order Euler method
            ylist(nn+1) = min(yval+dyy*dt,1);
            xlist(nn+1) = min(xval+dxx*dt,1);
  
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
end