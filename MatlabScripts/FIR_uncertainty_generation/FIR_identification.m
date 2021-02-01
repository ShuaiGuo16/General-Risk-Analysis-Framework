clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Identify FIR model from u' & q' signals
%   ===> Generate Fig. 2 in section 4.3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Sept. 2019
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b
% Toolbox: system identification toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FIR identification
load './data/data4SysID_A.mat'
[FIR_mean,FIR_var,FIR_model] = LengthFIR_est(0.4, data_timeseries);   % Use 400ms time series length
% save ./data/'FIR_fit_A.mat' FIR_model

% Post-processing
figure(1)

hold on
N = size(FIR_mean,2);
delta = data_timeseries.Ts;
time = 0:delta:(N-1)*delta;
stem(time,FIR_mean,'k','filled','LineWidth',1,'MarkerSize',6)
plot(time,FIR_mean'+1.96*sqrt(diag(FIR_var)),'r--','LineWidth',1.5)
plot(time,FIR_mean'-1.96*sqrt(diag(FIR_var)),'r--','LineWidth',1.5)
hold off

h = gca;
h.FontSize = 14;
xlabel('Time (s)')
ylabel('Amplitude')