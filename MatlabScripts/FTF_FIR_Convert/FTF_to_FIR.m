clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Generate reference FIR model from discrete frequency 
%            response data
% ALGORITHM
%   ===> Inverse Discrete Fourier Transform (iDFT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Oct. 2018
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load FTF data
Exp = load('./data/FTF_A.mat');
% Set impulse response model coefficient number
N = 250;    
Exp.phase = -Exp.phase;

%% Append the FTF data for conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency response data will be in 10Hz interval, with 250 FIR
% coefficients, FIR sampling interval will be 4e-4s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1-interpolate gain/phase (200-400 Hz) to ensure evenly distribution
Freq_interp = 200:10:400;
gain_app = interp1(Exp.Freq(21:end),Exp.gain(21:end),Freq_interp,'spline');
phase_app = interp1(Exp.Freq(21:end),Exp.phase(21:end),Freq_interp,'spline');
% 2-Append Freq/gain/phase (0-400 Hz)
Freq_list = [Exp.Freq(1:20);Freq_interp'];
gain_list = [Exp.gain(1:20);gain_app'];
phase_list = [Exp.phase(1:20);phase_app'];
% 3-Append Freq/gain/phase (400-1250 Hz)
gain_half = [gain_list;zeros(N/2-size(gain_list,1)+1,1)];  % For gain, add zero
phase_coeff = polyfit(Freq_interp,phase_app,1);           % For phase, linear extrapolation
phase_extrap = polyval(phase_coeff,410:10:1250);
phase_half = [phase_list;phase_extrap'];
% Generate conjugate
FTF_half = gain_half.*exp(1i*phase_half);
FTF_conj = conj(FTF_half);
FTF_full = [FTF_half;FTF_conj(end-1:-1:2)];

%% Convert back to FIR model
complex_vector = zeros(N,N);
for index_i = 1:N
    for index_j = 1:N
        complex_vector(index_i,index_j) = exp(1i*2*pi/N*(index_i-1)*(index_j-1));
    end
end
FIR_coeff = real(1/N*FTF_full.'*complex_vector);

%% Post-processing
% Impulse response
time = 0:4e-4:(N-1)*4e-4;    
figure(1)
stem(time,FIR_coeff,'k','filled','LineWidth',1.2,'MarkerSize',6)
xlabel('Time (s)')
ylabel('Magnitude')
title('Impuse Response')
h = gca;
h.FontSize = 14;

% Frequency response
Freq_evaluate = 0:10:400;
[gain, phase] = FTF_construct(FIR_coeff, 4e-4, Freq_evaluate');
figure(2)
subplot(2,1,1)
plot(Freq_evaluate,gain,'o','MarkerSize',8)
xlabel('Frequency (Hz)')
ylabel('Gain')
h = gca;
h.FontSize = 14;
subplot(2,1,2)
plot(Freq_evaluate,phase,'o','MarkerSize',8)
xlabel('Frequency (Hz)')
ylabel('Phase')
h = gca;
h.FontSize = 14;

% save './data/FIR_A.mat' FIR_coeff