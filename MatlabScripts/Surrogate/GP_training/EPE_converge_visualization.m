clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Generate Fig. 1 in Sec 4.2 GP model training
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Sept. 2019
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data
load './data/EPE_history.mat'
% Set convergence threshold as 5% of initial EPE value
iteration = 1:size(EPE_freq,1);
threshold_freq = EPE_freq(1)*0.05;
threshold_gref = EPE_gref(1)*0.05;


% For GP model approximating frequency
figure(1)
hold on
plot(iteration,EPE_freq,'LineWidth',1.2)
plot([0 size(EPE_freq,1)],[threshold_freq threshold_freq],'r--')
set(gca, 'YScale', 'log')
xlabel('Iterations')
ylabel('EPE')
axis([0 50 0 100])
h = gca;
h.FontSize = 14;

% For GP model approximating growth rate
figure(2)
hold on
plot(iteration,EPE_gref,'LineWidth',1.2)
plot([0 size(EPE_freq,1)],[threshold_gref threshold_gref],'r--')
set(gca, 'YScale', 'log')
xlabel('Iterations')
ylabel('EPE')
h = gca;
h.FontSize = 14;
axis([0 50 0 10000])