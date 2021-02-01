clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Monte Carlo simulation to assess the impact of FIR model 
%            coefficients uncertainty on the variation of modal frequency
%            and growth rate
%   ===> Apply Monte Carlo directly on Helmholtz solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITHM
%   ===> Generate random samples for FIR coefficients, R-in, R-out 
%            and damping coefficient alpha
%   ===> For each sample/realization, call Helmholtz solver to calculate
%            its corresponding modal frequency and growth rate
%   ===> Save samples & responses, post-processing
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
addpath('./SolverFunctions/')    % Helmholtz solver utility functions 
load './data/FIR_Fit_A.mat'      % Load uncertain FIR model
config_index = 11;                   % A11 configuration
s_init = 1i*112.28*2*pi;            % Initial values for iteration
ArgRs=0; ArgRn=pi;                % Phase of reflection coefficient: ArgRs (inlet), ArgRn (outlet)

%% Generate uncertain samples
MC = 20000;  % Monte Carlo sample number 
% 1-FIR parameter samples
FIR_mean = FIR_model.Numerator;
nb = size(FIR_mean,2);
FIR_var = getcov(FIR_model);
FIR_var = FIR_var(1:nb,1:nb);
FIR_MC = lhsnorm(FIR_mean,FIR_var,MC);
% 2-R_in/R_out/alpha samples (Latin-Hypercube)
X = lhsdesign(MC,3,'criterion','correlation');
R_in_MC = X(:,1)*0.3+0.7;  % Range: [0.7,1]
R_out_MC = X(:,2)*0.4+0.6; % Range: [0.6,1]
alpha_MC = X(:,3)*60+100;  % Range: [100,160]
% save './data/MC_samples_A11.mat' FIR_MC R_in_MC R_out_MC alpha_MC

%% Start MC routine
Eigenmode = zeros(MC,2);
for MC_loop = 1:MC
    % Helmholtz solver calculation
    [freq,greq] = Helmholtz_FIR('Secant',FIR_MC(MC_loop,:),...
        R_in_MC(MC_loop),ArgRs,R_out_MC(MC_loop),ArgRn,alpha_MC(MC_loop),config_index,s_init);
    Eigenmode(MC_loop,:) = [freq,greq];
    MC_loop
end
% save './data/Eigenmode_A11_MC.mat' Eigenmode    % Eigenvalue realizations

%% Post-processing
figure(1)
plot(Eigenmode(:,1),Eigenmode(:,2),'o')