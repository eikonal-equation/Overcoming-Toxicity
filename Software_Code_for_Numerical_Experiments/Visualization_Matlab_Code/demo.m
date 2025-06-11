%% Numerical Experiments Demo file for 
% "Overcoming toxicity: how non-antagonistic microbes manage to thrive in boom-and-bust environments" 
%
% This script is usued to demonstrate how to parse data by running our
% C++ code; then run our Matlab code to generate the figures in our manuscript.
%
% The user needs to run the C++ code themselves to obtain the data files
% before running this "Demo" script.
%
% Author: MingYi Wang, Cornell University
% Last modified: 06/2025
%
clear all;
close all;
clc;
%% Please note that the users need to run the Code sequentially
% as there are dependencies between the codes and data files.
File_valfn_4A = "Regular_constitutive_valuefn.dat";
File_parameters_4A = "Regular_dilutions_constitutive_DomainParameters.dat";
File_parameters_5A = "Regular_dilutions_myopic_optimal_DomainParameters.dat";

% For Fig. 7
File_parameters_7 = 'Sto_tactic_opt_DomainParameters.dat';
Filenames_7 = ["Perf_vertical_const_valuefn_rho650_lamb1000.dat",...
"Perf_vertical_tactic_valuefn_rho650_lamb1000.dat",...
"Sto_tactic_policy_lamb1000.dat", "Sto_strategic_vertical_policy_rho650_lamb1000.dat"];

% For Fig. 8
xc = 0.5; % x-coordinate of the point of interest 
yc = 0.1; % y-coordinate of the point of interest 
N = 1601; % number of points in each direction
lamb_list = [0.75,0.8,0.85,0.9,0.95,1.0]; % list of arrival rates (dilution frequencies)
rho_list = [0.5,0.55,0.6,0.65,0.7]; % list of survival rates (survival fractions)

% For Fig. S8
File_parameters_S8 = 'Sto_strategic_opt_vertical_DomainParameters.dat';
Filenames_S8 = ["Perf_vertical_const_valuefn_rho650_lamb1000.dat",...
"Perf_vertical_tactic_valuefn_rho650_lamb1000.dat", "Sto_strategic_vertical_valuefn_rho650_lamb1000.dat",...
"Sto_tactic_policy_lamb1000.dat", "Sto_strategic_vertical_policy_rho650_lamb1000.dat"];

Visualization_Fig3A();
Visualization_Fig3B();
Visualization_Fig3C();
Visualization_Fig3D();
Visualization_Fig4A(File_valfn_4A,File_parameters_4A);
Visualization_Fig4B();
Visualization_Fig4CtoF();
Visualization_Fig5(File_parameters_5A);
Visualization_Fig7(Filenames_7, File_parameters_7);
Visualization_Fig8(xc, yc, N, lamb_list, rho_list);
Visualization_FigS3();
Visualization_FigS4();
Visualization_FigS6();
Visualization_FigS8(Filenames_S8,File_parameters_S8);
Visualization_FigS12();
Visualization_FigS15();
