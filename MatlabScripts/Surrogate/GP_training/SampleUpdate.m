function [candidate,freq,gref,max_EPE] = SampleUpdate(Krig_model,MonteCarlo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Active learning scheme based on maximizing expected 
%            prediction error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%   ===> Krig_model: current GP model
%   ===> MonteCarlo: the pool of samples to select from
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
%   ===> candidate: the sample with the maximum EPE value
%                            select this sample to enrich the current
%                            training sample
%   ===> freq:  modal frequency of the selected sample
%   ===> gref:  modal growth rate of the selected sample
%   ===> max_EPE:   EPE value of the selected sample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Oct. 2019
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b
% Toolbox: UQLab (uqlab.com), optimization
% Ref: [1] H. Liu et al., An adaptive sampling approach for Kriging
% metamodeling by maximizing expected prediction error, Computers &
% Chemical Engineering, 106(2), 2017, 171-182
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate cross-validation error
cv_error_bank = CrossValidation(Krig_model);

% Search through all the candidate samples
[~,variance] = uq_evalModel(Krig_model,MonteCarlo);
cv_error = zeros(size(MonteCarlo,1),1);
for search_loop = 1:size(MonteCarlo,1)
    % Find the closest training sample
    index = Find_p_points(MonteCarlo(search_loop,:),Krig_model.ExpDesign.X); 
    cv_error(search_loop) = cv_error_bank(index);
end
% Calculate EPE values for all potential samples
EPE = cv_error+variance;

% To prevent selecting the same samples from previous adaptive steps
replicate_index = ismember(MonteCarlo,Krig_model.ExpDesign.X,'rows');
EPE = EPE(~replicate_index);
MonteCarlo = MonteCarlo(~replicate_index,:);
[max_EPE,max_index] = max(EPE);
candidate = MonteCarlo(max_index,:);

% Rescale the samples to facilitate Helmholtz solver calculation
candidate(1) = candidate(1)*2.5+0.5;
candidate(2) = candidate(2)*pi;
candidate(3) = candidate(3)*0.3+0.7;
candidate(4) = candidate(4)*0.4+0.6;
candidate(5) = candidate(5)*60+100;

% Helmholtz solver set up
config_index = 11;         % A11 configuration
s_init = 1i*112*2*pi;       % Initial value for iterative algorithm
ArgRs=0; ArgRn=pi;      % Phase of reflection coefficient: ArgRs (inlet), ArgRn (outlet)

[freq,gref] = Helmholtz_gain_phase('Secant',candidate(1),candidate(2),...
        candidate(3),ArgRs,candidate(4),ArgRn,candidate(5),config_index,s_init);

% Scale candidate sample back to [0 1]
candidate(1) = (candidate(1)-0.5)/2.5;
candidate(2) = candidate(2)/pi;
candidate(3) = (candidate(3)-0.7)/0.3;
candidate(4) = (candidate(4)-0.6)/0.4;
candidate(5) = (candidate(5)-100)/60;

end

