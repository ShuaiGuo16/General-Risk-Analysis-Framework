function [GP_Model] = GP_noise(training_X,training_Y,noise_var)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Stochastic Gaussian Process regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%   ===> training_X: matrix, each row represents a training sample
%   ===> training_Y: matrix, each row represents the response of the 
%                corresponding training sample
%   ===> noise_var: vector, variances of the sample response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
%   ===> GP_Model: trained GP model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Oct. 2018
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b, Optimization, Kriging scripts provided in the companion code of [1]
% Ref: [1] A. Forrester, Engineering Design via Surrogate Modelling: A Practical Guide
%             2008, Wiley.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1-Training data preparation
GP_Model.X = training_X;
GP_Model.y = training_Y;
GP_Model.noise = diag(noise_var);    % An add-on for covariance matrix

%% 2-Optimize under noise condition
k=3;        % 3 parameters to optimize
ub = [2 2 1];  lb = [-3 -3 0];   % Set up bounds

[x, val] = ...
    ga(@(x)likelihood_noise(x,GP_Model),k,[],[],[],[],lb,ub);

    GP_Model.Theta = x(1:2); 
    GP_Model.SigmaSqr = x(3); 
    [~,GP_Model.Psi,GP_Model.U] = likelihood_noise(x,GP_Model);
    
end

