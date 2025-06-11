function single_pop_dilutionIndex = singlePop_randDilution(num_dilution,r,xinit,rho,lamb)
    % This function simulates one sample of the population dynamics of a single strain
    % under random dilutions, given the initial population size (xinit),
    % the rescaled growth rate (r), the survival fraction (rho), and the
    % rescaled dilution frequency (lamb).
    % The function returns the population size list after a series of random dilutions.
    %
    %
    % Author: MingYi Wang, Cornell University
    % Last Modified: 06/2025
    %

% Theoretial result after a single dilution (Eq. [S3.1.4] in the SI Appendix)
pop_single = @(t,x,r) 1./(1 + ((1-x)./x)*exp(-r.*t));


single_pop_dilutionIndex=[xinit]; % initialize the vector
tlist = exprnd(1/lamb,num_dilution,1); % Generate random time intervals for dilutions

for ii = 1:num_dilution
    tt = tlist(ii);
    if ii == 1
        single_pop_dilutionIndex(ii + 1) = pop_single(tt,xinit,r);
    else
        single_pop_dilutionIndex(ii + 1) = pop_single(tt,rho*single_pop_dilutionIndex(ii),r);
    end
    
end
end
