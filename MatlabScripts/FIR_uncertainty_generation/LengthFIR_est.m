function [ FIR_vector, FIR_var, FIR_model ] = LengthFIR_est(time_length, dataseries)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Objective:
%       =====> Generate FIR based on the desired length                   
%  Input:
%       =====> time_length: Desired time length
%                    dataseries: Full length time series (type: iddata)
%  Output:
%       =====> FIR_vector: FIR coefficients (row vector)
%       =====> FIR_var: Variance matrix of FIR coefficients 
%                                  (column vector)
%       =====> FIR_model: FIR model in idtf format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Sept. 2019
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b
% Toolbox: system identification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start from 1s
begin_time = 1;    end_time = begin_time + time_length;

% Cut data according to the time length
selected_data = dataseries(ceil(begin_time/dataseries.Ts):ceil(end_time/dataseries.Ts));

% FIR identification
nb = 65; % number of non-zero impulse response coefficients
FIR_model = impulseest(selected_data,nb);  % idtf format

FIR_vector = FIR_model.Numerator;
FIR_var = getcov(FIR_model);
FIR_var = FIR_var(1:nb,1:nb);

end

