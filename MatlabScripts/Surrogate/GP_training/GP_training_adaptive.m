clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Train two GP models adaptively
%   ===> Both GP models input 
%                                     flame gain G, flame phase phi, 
%                                     magnitude of reflection coefficient at inlet R_in
%                                     magnitude of reflection coefficient at outlet R_out
%                                     system acoustic damping coefficient alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITHM
%   ===> Active learning strategy based on bias-variance decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Sept. 2019
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b
% Toolbox: UQLab (uqlab.com), optimization
% Ref: [1] H. Liu et al., An adaptive sampling approach for Kriging
% metamodeling by maximizing expected prediction error, Computers &
% Chemical Engineering, 106(2), 2017, 171-182
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('./SolverFunctions/')
uqlab

load './data/Candidate.mat'            % A pool of candidate samples for sample enrichment

%% (a) Step-1: Load initial training set 
load './data/Training.mat'
adp_training_X = training_X(1:50,:);     % 50 initial training number
adp_training_Y = training_Y(1:50,:);

%% (b) Step-2: Build two initial GP models
Metaopts_freq = CreateMetaOpts_Halton(adp_training_X, adp_training_Y(:,1));
GP_freq = uq_createModel(Metaopts_freq);

Metaopts_gref = CreateMetaOpts_Halton(adp_training_X, adp_training_Y(:,2));
GP_gref = uq_createModel(Metaopts_gref);

%% (c) Step-3: Start adaptive loop
max_iteration = 70;
iteration = 1;
% EPE values recorder
EPE_freq = zeros(max_iteration,1);
EPE_gref = zeros(max_iteration,1);

while iteration <= max_iteration
    
    % Show current index
    iteration
    
     % Employ adaptive scheme
      [freq_sample,freq_freq,freq_gref,EPE_freq(iteration)] = SampleUpdate(GP_freq,candidate);
      [gref_sample,gref_freq,gref_gref,EPE_gref(iteration)] = SampleUpdate(GP_gref,candidate);
      % Enrich training sample
      adp_training_X = [adp_training_X;freq_sample;gref_sample];
      adp_training_Y = [adp_training_Y;freq_freq,freq_gref;gref_freq,gref_gref];

     % Train new Kriging model
     clear GP_freq GP_gref   % delete old GP model
     
    Metaopts_freq = CreateMetaOpts_Halton(adp_training_X, adp_training_Y(:,1));
    GP_freq = uq_createModel(Metaopts_freq);

    Metaopts_gref = CreateMetaOpts_Halton(adp_training_X, adp_training_Y(:,2));
    GP_gref = uq_createModel(Metaopts_gref);
     
    % Iteration goes further
    iteration = iteration + 1;
    
end

% save './data/adp_training.mat' adp_training_X adp_training_Y
% save './data/EPE_history.mat' EPE_freq EPE_gref
