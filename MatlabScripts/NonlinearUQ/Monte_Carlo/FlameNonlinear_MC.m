clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Monte Carlo simulation to assess the impact of FDF dataset 
%        uncertainty on the variation of limit cycle amplitude and
%        frequency
%   ===> Apply Monte Carlo directly on Helmholtz solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITHM
%   ===> Track growth rate trajectory over velocity perturbation amplitude
%        to determin when limit cycle happens and its associated frequency and
%        amplitude
%   ===> For each amplitude level, firstly, extract the associated 26 gain
%        and phase data; secondly, employ rationalfit to fit a flame
%        model; finally, employ Helmholtz solver to calculate the modal
%        frequency and growth rate corresponding to this amplitude level
%   ===> Save frequency and growth rate trajectory for post-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Sept. 2019
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b
% Ref: [1] S. Guo et al, A Gaussian-Process-based framework for
% high-dimensional uncertainty quantification analysis in thermoacoustic
% instability prediction, 38th international symposium on Combustion, 2020,
% Adelaide, Australia.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
addpath('./SolverFunctions_rf/')     % Helmholtz solver utility functions 
Exp = load('./data/FDF_A_ori.mat');
s_init = 1i*112.28*2*pi;             % Initial values for iteration
config_index = 11;                   % A11 configuration
ArgRs=0; ArgRn=pi;                   % Phase of reflection coefficient: ArgRs (inlet), ArgRn (outlet)

%% Load uncertain samples
load './data/realization_gain.mat'
load './data/realization_phase.mat'
load './data/MC_samples_A11.mat'

MC = size(FDF_gain_MC,1); % Monte Carlo sample number
% Structure: each column corresponding to one amplitude level
freq = zeros(MC,length(Amp_realization));  gref = zeros(MC,length(Amp_realization));

for Amp_loop = 1:length(Amp_realization)
    
    % Extract 26 gain and phase dataset associated with the current
    % amplitude level
    FDF_gain_MC_Amp = ExtractAmp(Amp_realization(Amp_loop),FDF_gain_MC,Amp_realization,Freq_realization);   % Gain
    FDF_phase_MC_Amp = ExtractAmp(Amp_realization(Amp_loop),FDF_phase_MC,Amp_realization,Freq_realization);   % Gain
    
    % Fit the rational function as the flame model 
    fit = cell(MC,1);
    for fit_loop = 1:MC
        FDF_complex = FDF_gain_MC_Amp(fit_loop,:)'.*exp(-1i*FDF_phase_MC_Amp(fit_loop,:)');
        fit{fit_loop} = rationalfit(Exp.Freq(1:26),FDF_complex);
        fit_loop
    end
    
    iterator = ['Monte Carlo Simulation',' at',' ',num2str(Amp_loop),' Amplitude level'];
    % Set up parallel computing progress bar
    hbar = parfor_progressbar(MC,iterator);
    for MC_loop = 1:MC
        % Helmholtz solver calculation 
        [freq(MC_loop,Amp_loop),gref(MC_loop,Amp_loop)] = Helmholtz_rf('Secant',fit{MC_loop},...
            R_in_MC(MC_loop),ArgRs,R_out_MC(MC_loop),ArgRn,alpha_MC(MC_loop),config_index,s_init);
        hbar.iterate(1);
        
    end
    close(hbar);

end

% save './data/Eigenmode_nonlinear_MC.mat' freq gref
    
    