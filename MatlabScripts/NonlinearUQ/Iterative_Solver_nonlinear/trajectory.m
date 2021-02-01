function [Amp,LC_freq] = trajectory(freq,gref,Amp_list,changes)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Calculate limit cycle properties based on growth rate trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%   ===> freq: matrix, each row represent frequency trajectory of
%                    individual samples
%   ===> gref: matrix, each row represent growht rate trajectory of
%                    individual samples
%   ===> Amp_list: amplitude levels
%   ===> changes: a vector recording sign change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
%   ===> Amp: limit cycle amplitude
%   ===> LC_freq: limit cycle frequency 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Sept. 2019
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b
% Ref: [1] N. Noiray et al., A unified framework for nonlinear combustion
%        instability analysis based on the flame describing function, Journal of
%        Fluid Mechanics 615 (2008) pp. 139-167.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate derivative
gref_diff = diff(gref)./diff(Amp_list);

for i = 1:length(gref)-1
    if changes(i)==1 & gref(i)>0 & gref_diff(i)<0        
       LC_freq = freq(i)-gref(i)*...
            (freq(i)-freq(i+1))/(gref(i)-gref(i+1));
       Amp = Amp_list(i)-gref(i)*...
            (Amp_list(i)-Amp_list(i+1))/(gref(i)-gref(i+1));        
    end
end
        
end

