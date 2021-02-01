function [cv_error] = CrossValidation(Kriging_model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Calculate cross-validation error at each training sample location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%   ===> Kriging_model: current GP model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
%   ===> cv_error: vector, cross-validation error at each training sample location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Oct. 2018
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b
% Toolbox: UQlab (uqlab.com)
% Ref: [1] H. Liu et al., An adaptive sampling approach for Kriging
% metamodeling by maximizing expected prediction error, Computers &
% Chemical Engineering, 106(2), 2017, 171-182
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sample_number = Kriging_model.ExpDesign.NSamples;
cv_error = zeros(sample_number,1);

for i = 1:sample_number
    % Obtain R & R_inv
    R = Kriging_model.Internal.Kriging.GP.R;
    R_inverse = inv(R);
    
    % Construct d & H
    F = ones(sample_number,1);
    H = F*inv(F'*F)*F';
    
    beta = Kriging_model.Kriging.beta;
    d = Kriging_model.ExpDesign.Y - F*beta;
    
    % Calculate cross validation error
    cv_error(i) = R_inverse(i,:)*(d+H(:,i)*d(i)/(1-H(i,i)))/R_inverse(i,i);
end

cv_error = cv_error.^2;
