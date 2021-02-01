clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Compare results predicted by Helmholtz solver and 
%        proposed iterative GP-based scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETAILS
%   ===> Three situations exist: linearly/nonlinearly stable,
%        linearly unstable, linearly stable, nonlinearly unstable
%   ===> Linearly/nonlinearly stable case will not be shown in the
%        histogram plot since limit cycle is not occurred. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Sept. 2019
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b
% Toolbox: UQLab (uqlab.com), optimization
% Ref: [1] S. Guo et al, A Gaussian-Process-based framework for
% high-dimensional uncertainty quantification analysis in thermoacoustic
% instability prediction, 38th international symposium on Combustion, 2020,
% Adelaide, Australia.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load results
load './data/Eigenmode_nonlinear_GP.mat'
load './data/Eigenmode_nonlinear_MC.mat'
load './data/realization_gain.mat'

% Limit cycle calculation
MC = size(freq_iter,1);
[Amp_GP,LC_freq_GP,index_GP] = LC_calculator(freq,gref,Amp_realization);
[Amp_MC,LC_freq_MC,index_MC] = LC_calculator(freq_iter,gref_iter,Amp_realization);

% Compare histogram of limit cycle amplitude
figure(1)
hold on

% GP-based iterative scheme results
% Only plot samples which lead to limit cycle oscillation
Amp_GP = Amp_GP(Amp_GP~=-0.08); unstable_number_GP = length(Amp_GP);
pts = linspace(min(Amp_GP),max(Amp_GP),200);
[f_Pre,xi_Pre] = ksdensity(Amp_GP,pts);
% Correct probability density
plot(xi_Pre,f_Pre*unstable_number_GP/MC,'k-','LineWidth',2)

% Helmholtz results
edges = -0.08:0.05:0.84;
histogram(Amp_MC,edges,'Normalization','pdf')

% Experimental results
plot(0.35,0,'rd','markerFaceColor','r','MarkerSize',10)

% Minimal amplitude line
plot([0.07 0.07],[0 0.45],'r--','LineWidth',4)

hold off
yticks([0 0.15 0.3 0.45])
axis([0 0.84 0 0.45])
h = gca;
h.FontSize = 14;

figure(2)
hold on

% GP-based iterative scheme results
% Only plot samples which lead to limit cycle oscillation
LC_freq_GP = LC_freq_GP(LC_freq_GP~=50); unstable_number_GP = length(LC_freq_GP);
pts = linspace(102,120,100);
[f_Pre,xi_Pre] = ksdensity(LC_freq_GP,pts);
plot(xi_Pre,f_Pre*unstable_number_GP/MC,'k-','LineWidth',2)

% Helmholtz results
edges = 50:0.8:128;
histogram(LC_freq_MC,edges,'Normalization','pdf')

% Experimental results
plot(115,0,'rd','markerFaceColor','r','MarkerSize',10)

hold off
axis([100 120 0 0.02])
yticks([0 0.005 0.01 0.015 0.02])
h = gca;
h.FontSize = 14;