% Function that returns gain and phase

% Procedure: Values are given for Gain and Phase in freqency range from
% 0 to 360 HZ. These 33 point values are used to build a spline. Since the
% Gain and Phase are given for six different pertubation levels the two
% confining splines (above and below) are chosen and a linear interpolation
% is performed.

function [gain, phase] = FDFrand(omega, amplitude,n,FLAME)

% Load Freq, Gain and Phase, initialize for debugging
Freq=[];
Gain=[];
Phase=[];

if strcmp(FLAME,'A');
    load('FDF_A.mat'); clear Gain Phase
    load('FDF_Arand.mat');
    ampl_val = [0, 0.07, 0.15, 0.30, 0.41, 0.51, 0.71, 1];
end

if strcmp(FLAME,'B');
    load('FDF_B.mat'); clear Gain Phase
    load('FDF_Brand.mat');
    ampl_val = [0, 0.04, 0.09, 0.19, 0.36, 0.48, 0.60, 1];
end

Gain =G(:,:,n);
Phase=Ph(:,:,n);

% Errors
if amplitude < 0 || amplitude > 1 
    error('Amplitude out of range')
end

if omega < 0 || omega/(2*pi) > Freq(end)
    error('Omega out of range')
end

for i=0:6
    
    if amplitude >= ampl_val(i+1) && amplitude <= ampl_val(i+2);
        lower=i;
        upper=i+1;
    end
end

if lower == 0
    glowspline = spline(Freq,Gain(:,lower+1));
    plowspline = spline(Freq,Phase(:,lower+1));
else
    glowspline = spline(Freq,Gain(:,lower));
    plowspline = spline(Freq,Phase(:,lower));
end

if upper == 7
    guppspline = spline(Freq,Gain(:,upper-1));
    puppspline = spline(Freq,Phase(:,upper-1));
else
    guppspline = spline(Freq,Gain(:,upper));
    puppspline = spline(Freq,Phase(:,upper));
end

freq=omega/(2*pi);

glowval = ppval(glowspline,freq);
guppval = ppval(guppspline,freq);

plowval = ppval(plowspline,freq);
puppval = ppval(puppspline,freq);

% Interpolation between upp and low
gain  = glowval + (guppval-glowval)/(ampl_val(upper+1)-ampl_val(lower+1)) * (amplitude-ampl_val(lower+1));

phase = plowval + (puppval-plowval)/(ampl_val(upper+1)-ampl_val(lower+1)) * (amplitude-ampl_val(lower+1));

end

       

