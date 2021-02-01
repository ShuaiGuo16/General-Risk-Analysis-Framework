function [Fun] = Eigenmode_solver_FDF(x,fit,R_in,R_out,alpha,GP_freq,GP_gref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Specify options for Gaussian Process model training
%            in UQLab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%   ===>   x:  x(1)=frequency,  x(2)=growth rate 
%   ===>   fit:  rationfit object
%   ===>   R_in: reflection coefficient at inlet
%   ===>   R_out: reflection coefficient at outlet
%   ===>   alpha: damping coefficient
%   ===>   GP_freq:  GP surrogate model for frequency
%   ===>   GP_gref:  GP surrogate model for growth rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
%   ===> Fun:   Fun(1)=frequency match, Fun(2)=growth rate match
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Sept. 2019
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b
% Ref: [1] S. Guo et al, A Gaussian-Process-based framework for
% high-dimensional uncertainty quantification analysis in thermoacoustic
% instability prediction, 38th international symposium on Combustion, 2020,
% Adelaide, Australia.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Nominal equations
s = x(2)+1i*x(1)*2*pi;
FTF = sum(fit.C./(s-fit.A));

% Rescale the samples
input = zeros(1,5);
input(1) = (abs(FTF)-0.5)/2.5;     % Gain
input(2) = angle(FTF)/pi;    % Phase
input(3) = (R_in-0.7)/0.3;           % R_in
input(4) = (R_out-0.6)/0.4;          % R_out
input(5) = (alpha-100)/60;           % alpha
% Nominal prediction
predict_freq = uq_evalModel(GP_freq,input);
predict_gref = uq_evalModel(GP_gref,input);

% Construct equations
Fun(1) = x(1)-predict_freq;
Fun(2) = x(2)-predict_gref;


end

