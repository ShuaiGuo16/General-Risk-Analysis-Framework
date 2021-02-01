function [combined_array] = mesh2array(dim1,dim2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Construct a 2D matrix based on two vectors 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%   ===> dim1: vector, the first dimension content
%   ===> dim2: vector, the second dimension content
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
%   ===> combined_array: the resulting 2D matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Oct. 2018
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b
% Ref: [1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
combined_array = zeros(length(dim1)*length(dim2),2);
for i = 1:size(dim1,1)
    combined_array((i-1)*length(dim2)+1:i*length(dim2),:) = [dim1(i)*ones(length(dim2),1),dim2];
end

end

