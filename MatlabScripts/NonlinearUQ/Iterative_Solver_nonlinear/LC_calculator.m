function [limit_cycle_amp,limit_cycle_freq,index] = LC_calculator(freq,gref,Amp_list)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
%   ===> limit_cycle_amp: vector, limit cycle amplitude for each sample
%   ===> limit_cycle_freq: vector, limit cycle frequency for each sample
%   ===> index: logic matrix, recording the growth rate trajectory pattern (type 1, 2a, 2b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Sept. 2019
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b
% Ref: [1] N. Noiray et al., A unified framework for nonlinear combustion
%        instability analysis based on the flame describing function, Journal of
%        Fluid Mechanics 615 (2008) pp. 139-167.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MC = size(freq,1);
sign_change = zeros(MC,1);
changes = zeros(MC,length(Amp_list)-1);
% For each sample, record how many sign changes in growht rate trajectory
for i = 1:MC
    pos = gref(i,:)>0;
    changes(i,:) = xor(pos(1:end-1),pos(2:end));
    sign_change(i) = sum(changes(i,:));
end

% 4 situations
index0 = sign_change==0 & gref(:,1)<0;    % Always stable (Type 2b)
index1 = sign_change==1;                        % One cross (Type 1)
index2 = sign_change==2;                        % Two cross (Type 2a)
index3 = sign_change==0 & gref(:,1)>0;    % Always unstable (Type 1, extrapolation required)
index = [index0,index1,index2,index3];

limit_cycle_amp = ones(MC,1);
limit_cycle_freq = ones(MC,1);

for i = 1:MC
    
    if index0(i)
        limit_cycle_amp(i) = -0.08;   % Assign some arbitrary negative value for amplitude
        limit_cycle_freq(i) = 50;       % Assign some arbitrary small value for frequency
   
    elseif index3(i)   % Always unstable, extrapolation required
        limit_cycle_freq(i) = freq(i,end-1)-gref(i,end-1)*...
            (freq(i,end-1)-freq(i,end))/(gref(i,end-1)-gref(i,end));
        limit_cycle_amp(i) = Amp_list(end-1)-gref(i,end-1)*...
            (Amp_list(end-1)-Amp_list(end))/(gref(i,end-1)-gref(i,end));
        
    else
        [limit_cycle_amp(i),limit_cycle_freq(i)] = ...
            trajectory(freq(i,:),gref(i,:),Amp_list,changes(i,:));
    end
        
end

end

