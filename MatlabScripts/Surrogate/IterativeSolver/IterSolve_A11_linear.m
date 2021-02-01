clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Iterative scheme to propagate FIR uncertainty to modal 
%            frequency and growth rate predictions
%   ===> Generate Fig. 3 in Sec 4.3 Linear thermoacoustic UQ analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITHM
%   ===> Construct GP models based on the adaptively selected training
%            samples
%   ===> Since GP models are just an algebraic model, it is possible to 
%            analytically derive the gradients to facilitate a faster
%            calculation to Eq. (3)-(4) in Sec 2 Surrogate-based UQ
%            strategy
%   ===> Nevertheless, it is also possible to simply disable gradient
%            calculation, which will lead to a slightly slower convergence
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
uqlab

% Build surrogate model from previously determined training samples
load './data/adp_training.mat'
required_sample_num = 50*2+50;    % Required number of samples to reach target accuracy

% GP model for frequency
Metaopts_freq = CreateMetaOpts_Halton(adp_training_X(1:required_sample_num,:),...
    adp_training_Y(1:required_sample_num,1));
GP_freq = uq_createModel(Metaopts_freq);
% GP model for growth rate
Metaopts_gref = CreateMetaOpts_Halton(adp_training_X(1:required_sample_num,:),...
    adp_training_Y(1:required_sample_num,2));
GP_gref = uq_createModel(Metaopts_gref);
% Extract gradient info
freq_R = GP_freq.Internal.Kriging.GP.R;  
gref_R = GP_gref.Internal.Kriging.GP.R;  
freq_gradConst = freq_R\(GP_freq.ExpDesign.Y-GP_freq.Kriging.beta*GP_freq.Internal.Kriging.Trend.F);
gref_gradConst = gref_R\(GP_gref.ExpDesign.Y-GP_gref.Kriging.beta*GP_gref.Internal.Kriging.Trend.F);

% Load samples/reference responses
load './data/Eigenmode_A11_MC.mat'
load './data/MC_samples_A11.mat'
MC_num = size(FIR_MC,1);
freq_iter = ones(MC_num,1);  gref_iter = ones(MC_num,1);

x0 = [mean(adp_training_Y(1:required_sample_num,1)),0];  % Initial value for iterative scheme
remaining_step = MC_num;
% Enable gradient
options = optimoptions(@fsolve,'Display','off','Algorithm','trust-region',...
    'SpecifyObjectiveGradient',true);

for solve_loop = 1:MC_num
    EigenFun = @(x) Eigenmode_solver_gradient(x,FIR_MC(solve_loop,:),R_in_MC(solve_loop),...
        R_out_MC(solve_loop),alpha_MC(solve_loop),GP_freq,GP_gref,freq_gradConst,gref_gradConst);
    Eigen = fsolve(EigenFun,x0,options);
    freq_iter(solve_loop) = Eigen(1);
    gref_iter(solve_loop) = Eigen(2);
    remaining_step = remaining_step-1
end

% Post-processing
figure(1)
hold on
pts = min(freq_iter):0.01:max(freq_iter);
[f_Pre,xi_Pre] = ksdensity(freq_iter,pts);
plot(xi_Pre,f_Pre,'k-','LineWidth',2)
H = histogram(Eigenmode(:,1),'Normalization','pdf');
hold off

set(H,'FaceColor',[0.75 0.75 0.75]);
xlabel('Frequency (Hz)')
ylabel('PDF')
% yticks([0 0.02 0.04 0.06 0.08])
% axis([250 285 0 0.08])
h = gca;
h.FontSize = 12;

figure(2)
hold on
pts = min(gref_iter):0.01:max(gref_iter);
[f_Pre,xi_Pre] = ksdensity(gref_iter,pts);
plot(xi_Pre,f_Pre,'k-','LineWidth',2)
H = histogram(Eigenmode(:,2),'Normalization','pdf');
hold off

set(H,'FaceColor',[0.75 0.75 0.75]);
xlabel('Growth Rate (rad/s)')
ylabel('PDF')
% yticks([0 0.02 0.04 0.06 0.08])
% axis([250 285 0 0.08])
h = gca;
h.FontSize = 12;