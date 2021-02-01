function [FDF_realizations_Amp] = ExtractAmp(target_Amp,FDF_realizations,Amp_list,Freq_list)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Extract FDF data associated with the specific amplitude level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%   ===> target_Amp: current amplitude level
%   ===> FDF_realizations: full FDF dataset (260)
%   ===> Amp_list: amplitude levels
%   ===> Freq_list: frequency levels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
%   ===> FDF_realizations_Amp: FDF data associated with the target amplitude level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Oct. 2018
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b
% Ref: [1] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_index = find(Amp_list==target_Amp);
FDF_realizations_Amp = FDF_realizations(:,start_index:length(Amp_list):start_index+(length(Freq_list)-1)*length(Amp_list));

end

