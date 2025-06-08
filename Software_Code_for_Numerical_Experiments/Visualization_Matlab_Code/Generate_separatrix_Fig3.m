function Generate_separatrix_Fig3()
% This function generates the separatrix for the model in Fig. 3 of the paper

%% separatrix
death_arr = [0.05,0.215,0.3,0.43,0.5];

for ii = 1:length(death_arr)
    epsilon = 0.2; % The cost of constitutive toxin production
    gamma = 1.0; % The rescaled cost of toxin-production rate
    rks = 0.85; % The rescaled ratio of the growth rate of the killer to that of the sensitive  
    d1 = death_arr(ii); % rescaled basal death rate of the killer
    d2 = d1; % rescaled basal death rate of the sensitive

    % Drift part of the population growth dynamics (Eq.[S2.1] in the SI)
    fk = @(nk,ns,a) rks.*(1-epsilon.*a).*nk.*(1-nk-ns) -d1*nk;
    fs = @(nk,ns,a) ns.*(1-nk-ns)-a.*gamma.*nk.*ns - d2*ns;
    % Fixed points of the system
    pt1 = 1 - d1/(rks*(1-epsilon));
    pt2 = 1 - d2;
    pt3 = d1/(gamma*rks*(1-epsilon)) - d2/gamma;
    pt4 = (d2 + gamma)/gamma - (d1*(1+gamma))/(gamma*rks*(1-epsilon));
  
    [y_stable_forward, y_unstable_forward,y_stable_backward,y_unstable_backward] = ...
        hyperbolic_traj2(pt3,pt4, a,fk,fs,epsilon,gamma,rks);
    str = sprintf('separatrix_%g.mat',ii);
    save(str,'y_unstable_backward',"y_stable_backward",'y_unstable_forward',"y_stable_forward","pt1","pt2","pt3","pt4");

end


%% subfunctions
function [value, isterminal, direction] = myEventFunction(t, y)
    % This function defines an event for the ODE solver to stop integration
    % when the condition x + y = 1 is met.
    % Input:
    %   t: time
    %   y: solution vector
    % Output:
    %   value: value of the event function
    %   isterminal: flag indicating whether the integration should stop
    %   direction: direction of the event function

    % Event function to stop integration when x + y = 1
    value = y(1) + y(2) - 1; % Detect when x + y = 1
    isterminal = 1; % Stop the integration
    direction = 0; % Detect all zero crossings

end

function [y_stable_forward, y_unstable_forward,y_stable_backward,y_unstable_backward] = ...
    hyperbolic_traj2(pt3,pt4,a,fk,fs,epsilon,gamma,rks)
    % This function computes the stable and unstable manifolds of the system
    % given the fixed points pt3 and pt4.
    % Input:
    %   pt3: x-coordinate of the fixed point
    %   pt4: y-coordinate of the fixed point
    %   a: toxin production rate
    %   fk: function handle for the drift of the killer population
    %   fs: function handle for the drift of the sensitive population
    %   epsilon: cost of constitutive toxin production
    %   gamma: rescaled cost of toxin-production rate
    %   rks: rescaled ratio of the growth rate of the killer to that of the sensitive
    % Output:
    %   y_stable_forward: stable manifold forward trajectory
    %   y_unstable_forward: unstable manifold forward trajectory
    %   y_stable_backward: stable manifold backward trajectory
    %   y_unstable_backward: unstable manifold backward 
    
    tspan = [0 -200];
    tspan_foward = [0 2000];
    % Options for ode45 to include the event function
    options = odeset('Events', @myEventFunction);

    % Stable manifold
    [t_stable_forward, y_stable_forward] = ode45((@(t,y)[fk(y(1),y(2),a);fs(y(1),y(2),a)]), tspan, [pt3 + 1e-3,pt4 + 1e-3],options);
    [t_stable_backward, y_stable_backward] = ode45((@(t,y)[fk(y(1),y(2),a);fs(y(1),y(2),a)]), tspan, [pt3 - 1e-3,pt4 - 1e-3]);

    % Unstable manifold
    [t_unstable_forward, y_unstable_forward] = ode45((@(t,y)[fk(y(1),y(2),a);fs(y(1),y(2),a)]), tspan_foward, [pt3 - 1e-3,pt4 + 1e-3]);
    [t_unstable_backward, y_unstable_backward] = ode45((@(t,y)[fk(y(1),y(2),a);fs(y(1),y(2),a)]), tspan_foward, [pt3 + 1e-3,pt4 - 1e-3]);

end

end