function [frac,pop] = sampling_binomial(f,N,K,rho)
tot_pop = round(N*K);
pop1 = round(f*tot_pop);
pop2 = tot_pop - pop1;

pop1_left = binornd(pop1,rho);
pop2_left = binornd(pop2,rho);
pop_left = pop1_left + pop2_left;

frac = pop1_left/(pop_left);
pop = pop_left/K;

end