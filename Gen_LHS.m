%% This portion of the code generates the LHS

clear all 
close all 
clc

% Chem-atr-rad-land Chem-atr-rad-col initial-starch init-energy-land initial-energy-col max-col-starch-harvest max-land-starch-harvest replenish-style col-fert-thresh land-fert-thresh  
% prop-fert-cutoff freq-starch-replenish prop-starch-replenish step-size
% land-movement-cost col-movement-cost H1P1
lb = [1 1   0   0   0  0  0  0  0  0 0   1 0 0.05  0  0  0]; % Specify lower bounds 
ub = [6 6 100 100 100 20 20  1 20 20 1 720 1 5 10 10 10]; % Specify upper bounds
ns = 10000; % number of samples
p = 17; % number of parameters
xn = lhsdesign(ns,p); % generate normalized design
x = bsxfun(@plus,lb,bsxfun(@times,xn,(ub-lb)));
writematrix(x,'LHS.csv')
