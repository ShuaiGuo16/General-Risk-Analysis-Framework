clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Generate samples for Monte Carlo simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Sept. 2019
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MC = 20000;    % Monte Carlo sample number
% 1-Flame parameter samples (gain)
load './data/realization_gain.mat'
FDF_gain_MC = mvnrnd(f,C,MC);
% 2-Flame parameter samples (phase)
load './data/realization_phase.mat'
FDF_phase_MC = mvnrnd(f,C,MC);
% 2-R_in/R_out/alpha samples
X = lhsdesign(MC,3,'criterion','correlation');
R_in_MC = X(:,1)*0.1+0.9;  % Range: [0.9,1]
R_out_MC = X(:,2)*0.3+0.7; % Range: [0.7,1]
alpha_MC = X(:,3)*50+110;  % Range: [110,160]
% save './data/MC_samples_A11.mat'  FDF_gain_MC FDF_phase_MC R_in_MC R_out_MC alpha_MC