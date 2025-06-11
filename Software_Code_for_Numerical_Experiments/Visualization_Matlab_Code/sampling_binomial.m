function [frac,pop] = sampling_binomial(f,N,K,rho)
% This function use a Binomial sampling method to obtain the
% post-dilution populations to initialize the next round.
% Inputs:
%   f: fraction of the killer strain in the population before dilution
%   N: total population before dilution
%   K: shared carrying capacity
%   rho: rescaled survival rate of the killer strain    

tot_pop = round(N*K); % Total population before dilution (round to nearest integer)
Killer_pop = round(f*tot_pop); % Killer population before dilution
Sensitive_pop = tot_pop - Killer_pop; % Sensitive population before dilution

Killer_pop_left = binornd(Killer_pop,rho); % Killer population after dilution using Binomial sampling
Sensitive_pop_left = binornd(Sensitive_pop,rho); % Sensitive population after dilution using Binomial sampling
pop_left = Killer_pop_left + Sensitive_pop_left; % Absolute Total population after dilution

frac = Killer_pop_left/(pop_left); % Fraction of the killer strain after dilution
pop = pop_left/K; % Normalized total population after dilution

end