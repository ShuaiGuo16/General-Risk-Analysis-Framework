clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Generate initial training samples for GP model training
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Sept. 2019
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b
% Toolbox: UQLab (uqlab.com)
% Ref: [1] S. Guo et al, A Gaussian-Process-based framework for
% high-dimensional uncertainty quantification analysis in thermoacoustic
% instability prediction, 38th international symposium on Combustion, 2020,
% Adelaide, Australia.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem set-up
addpath('./SolverFunctions/')    % Add solver location
config_index = 11;                    % A11 configuration
s_init = 1i*112*2*pi;                  % Initial value for iterative algorithm
ArgRs=0; ArgRn=pi;                 % Phase of reflection coefficient: ArgRs (inlet), ArgRn (outlet)

% Generate Halton samples
uqlab
for ii = 1:5
    Input.Marginals(ii).Type = 'Uniform';
    Input.Marginals(ii).Parameters = [0 1];
end
ParaInput = uq_createInput(Input);

% Start with 50 samples
sample_number = 50;
training_X = uq_getSample(sample_number,'Halton');
training_Y = zeros(sample_number,2);
% Scale to solver-acceptable values
G = 2.5*training_X(:,1)+0.5;        % Flame gain [0.5,3]
phi = pi*training_X(:,2);                % Flame phase [0 pi]
R_in = 0.3*training_X(:,3)+0.7;     % R-in [0.7,1]
R_out = 0.4*training_X(:,4)+0.6;   % R-out [0.6,1]
alpha = 60*training_X(:,5)+100;    % Range: [100,160]


% Solver calculation
for cal_index = 1:sample_number
    cal_index
    [training_Y(cal_index,1),training_Y(cal_index,2)] = Helmholtz_gain_phase('Secant',G(cal_index),phi(cal_index),...
        R_in(cal_index),ArgRs,R_out(cal_index),ArgRn,alpha(cal_index),config_index,s_init);
end

% Post-processing
figure(1)
plot(training_Y(:,1),training_Y(:,2),'ob')
% save './data/Training.mat' training_X training_Y