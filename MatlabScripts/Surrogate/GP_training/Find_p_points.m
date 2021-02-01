function [index] = Find_p_points(target_sample,sample_data_bank)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Find the nearest sample among the data bank for target sample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%   ===> target_sample: sample point needs to be processed
%   ===> sample_data_bank: the provided data bank of samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
%   ===> index: the index of the nearest sample in the data bank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Oct. 2018
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pool_num = size(sample_data_bank,1);
dis = sqrt(sum((repmat(target_sample,pool_num,1)-sample_data_bank).^2,2));

% sort distance
[Y,index] = min(dis);

end

