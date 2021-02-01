clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Iterative scheme to propagate FDF uncertainty to limit 
%        cycle frequency and amplitude predictions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITHM
%   ===> Construct GP models based on the adaptively selected training
%        samples
%   ===> For each amplitude level, firstly, extract the associated 26 gain
%        and phase data; secondly, employ rationalfit to fit a flame
%        model; finally, employ proposed iterative scheme to calculate the modal
%        frequency and growth rate corresponding to this amplitude level
%   ===> Save frequency and growth rate trajectory for post-processing
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

% Initialization
uqlab

% Build surrogate model
load './data/adp_training.mat'
required_sample_num = 50*2+50;    % Required number of samples to reach target accuracy
Metaopts_freq = CreateMetaOpts_Halton(adp_training_X(1:required_sample_num,:),...
    adp_training_Y(1:required_sample_num,1));
GP_freq = uq_createModel(Metaopts_freq);

Metaopts_gref = CreateMetaOpts_Halton(adp_training_X(1:required_sample_num,:),...
    adp_training_Y(1:required_sample_num,2));
GP_gref = uq_createModel(Metaopts_gref);

% Load samples
load './data/MC_samples_A11.mat'
load './data/realization_gain.mat'
load './data/realization_phase.mat'
Exp = load('./data/FDF_A_ori.mat');
    
    % Set Monte Carlo parameters
    MC_num = 100;
    % Structure: each column corresponding to one amplitude level
    freq_iter = zeros(MC_num,length(Amp_realization)); gref_iter = zeros(MC_num,length(Amp_realization));
    x0 = 112;
    % Disable gradient
    options = optimoptions(@fsolve,'Display','off','Algorithm','trust-region');

    for Amp_loop = 1:length(Amp_realization)

        if Amp_loop == 1
            x0 = x0*ones(MC_num,1);  x0 = [x0,zeros(MC_num,1)];
        else
            x0 = [freq_iter(:,Amp_loop-1),gref_iter(:,Amp_loop-1)];
        end
        
        % Extract 26 gain and phase dataset associated with the current
        % amplitude level
        FDF_gain_MC_Amp = ExtractAmp(Amp_realization(Amp_loop),FDF_gain_MC,Amp_realization,Freq_realization);   % Gain
        FDF_phase_MC_Amp = ExtractAmp(Amp_realization(Amp_loop),FDF_phase_MC,Amp_realization,Freq_realization); % Phase

        % Fit the rational function as the flame model 
        fit = cell(MC_num,1);
        for fit_loop = 1:MC_num
            FDF_complex = FDF_gain_MC_Amp(fit_loop,:)'.*exp(-1i*FDF_phase_MC_Amp(fit_loop,:)');
            fit{fit_loop} = rationalfit(Exp.Freq(1:26),FDF_complex);
            fit_loop
        end

        iterator = ['Monte Carlo Simulation',' at',' ',num2str(Amp_loop),' Amplitude level'];
        hbar = parfor_progressbar(MC_num,iterator);

        parfor solve_loop = 1:MC_num
            EigenFun = @(x) Eigenmode_solver_FDF(x,fit{solve_loop},R_in_MC(solve_loop),...
                R_out_MC(solve_loop),alpha_MC(solve_loop),GP_freq,GP_gref);
            Eigen = fsolve(EigenFun,x0(solve_loop,:),options);
            freq_iter(solve_loop,Amp_loop) = Eigen(1);
            gref_iter(solve_loop,Amp_loop) = Eigen(2);
            hbar.iterate(1);
        end
        close(hbar);

    end

% Output results
% save './data/Eigenmode_nonlinear_GP.mat' freq_iter gref_iter