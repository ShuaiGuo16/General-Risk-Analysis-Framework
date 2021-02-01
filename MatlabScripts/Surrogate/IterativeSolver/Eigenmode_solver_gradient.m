function [Fun,J] = Eigenmode_solver_gradient(x,FIR,R_in,R_out,alpha,GP_freq,GP_gref,freq_gradConst,gref_gradConst)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the system of nonlinear equations Eq. (2) in Sec 2 
% Surrogate-based UQ strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
%      x:  x(1)=frequency,  x(2)=growth rate 
%      FIR:  FIR coefficients
%      R_in: reflection coefficient at inlet
%      R_out: reflection coefficient at outlet
%      alpha: acoustic damping coefficient
%      GP_freq:  GP model for frequency
%      GP_gref:  GP model for growth rate
%      freq_gradConst:  constant term in taking derivatives of frequency GP model
%      freq_gradConst:  constant term in taking derivatives of growth rate GP model
% Output:
%      Fun: Fun(1)=frequency match, Eq. (3), Fun(2)=growth rate match, Eq. (4)
%      J:     Jacobian matrix (gradients)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Sept. 2019
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b
% Ref: [1] S. Guo et al, A Gaussian-Process-based framework for
% high-dimensional uncertainty quantification analysis in thermoacoustic
% instability prediction, 38th international symposium on Combustion, 2020,
% Adelaide, Australia.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Nominal equations
delta_t = 4e-4; nb = size(FIR,2);
s = x(2)+1i*x(1)*2*pi;
k = [1:nb]';
FTF = FIR*exp(-(k-1)*delta_t*s);

% Rescale the samples
input = zeros(1,5);
input(1) = (abs(FTF)-0.5)/2.5;     % Gain
input(2) = angle(FTF)/pi;    % Phase
input(3) = (R_in-0.7)/0.3;           % R_in
input(4) = (R_out-0.6)/0.4;          % R_out
input(5) = (alpha-100)/60;           % alpha
% Nominal prediction
predict_freq = uq_evalModel(GP_freq,input);
predict_gref = uq_evalModel(GP_gref,input);

% Construct equations
Fun(1) = x(1)-predict_freq;
Fun(2) = x(2)-predict_gref;

%% Jacobian matrix
% From GP model side (Freq GP model)
theta_VecSqr = 1./GP_freq.Kriging.theta.^2;
freq_r = exp((repmat(input,GP_freq.ExpDesign.NSamples,1)-GP_freq.ExpDesign.X).^2*theta_VecSqr'*-0.5);
theta_G = GP_freq.Kriging.theta(1);
theta_phi = GP_freq.Kriging.theta(2);
df1dG = (-1/theta_G^2*freq_r.*(input(1)-GP_freq.ExpDesign.X(:,1)))'*freq_gradConst;
df1dphi = (-1/theta_phi^2*freq_r.*(input(2)-GP_freq.ExpDesign.X(:,2)))'*freq_gradConst;

% From GP model side (Gref GP model)
theta_VecSqr = 1./GP_gref.Kriging.theta.^2;
gref_r = exp((repmat(input,GP_gref.ExpDesign.NSamples,1)-GP_gref.ExpDesign.X).^2*theta_VecSqr'*-0.5);
theta_G = GP_gref.Kriging.theta(1);
theta_phi = GP_gref.Kriging.theta(2);
df2dG = (-1/theta_G^2*gref_r.*(input(1)-GP_gref.ExpDesign.X(:,1)))'*gref_gradConst;
df2dphi = (-1/theta_phi^2*gref_r.*(input(2)-GP_gref.ExpDesign.X(:,2)))'*gref_gradConst;

% From FIR side (frequency derivative)
A_f = [cos(angle(FTF)),-sin(angle(FTF))*abs(FTF);sin(angle(FTF)),cos(angle(FTF))*abs(FTF)];
dFTFdf = FIR*(-1i*delta_t*(k-1).*exp(-(k-1)*delta_t*s)*2*pi);
B_f = [real(dFTFdf);imag(dFTFdf)];
df = A_f\B_f;   df(1) = df(1)/2.5;  df(2) = df(2)/pi;   % Scale the derivatives

A_g = A_f;
dFTFdg = FIR*(-delta_t*(k-1).*exp(-(k-1)*delta_t*s));
B_g = [real(dFTFdg);imag(dFTFdg)];
dg = A_g\B_g;   dg(1) = dg(1)/2.5;  dg(2) = dg(2)/pi; 

J(1,1) = 1-(df1dG*df(1)+df1dphi*df(2));
J(1,2) = -(df1dG*dg(1)+df1dphi*dg(2));
J(2,1) = -(df2dG*df(1)+df2dphi*df(2));
J(2,2) = 1-(df2dG*dg(1)+df2dphi*dg(2));

end

