function [ Metaopts ] = CreateMetaOpts_Halton( X, Y )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%   ===> Specify options for Gaussian Process model training
%            in UQLab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%   ===> X: matrix, each row represents a training sample
%   ===> Y: matrix, each row represents the four responses of the 
%                corresponding training sample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
%   ===> Metaopts: Options for Gaussian Process model training
%            in UQLab (for details please see "Kriging UserManual" of
%            UQLab in www.uqlab.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by S. Guo (TUM), Oct. 2018
% Email: guo@tfd.mw.tum.de
% Version: MATLAB R2018b
% Ref: [1] S. Guo, C. F. Silva, W. Polifke, "Efficient robust design for
% thermoacoustic instability analysis: A Gaussian process approach",
% 2019, ASME Turo Expo, Phoenix, USA, GT2019-90732
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Metaopts.Type = 'Metamodel';
Metaopts.MetaType = 'Kriging';
Metaopts.ExpDesign.Sampling = 'user';
Metaopts.ExpDesign.X = X;
Metaopts.ExpDesign.Y = Y;

Metaopts.Scaling = 0;     % impose no scaling
Metaopts.Trend.Type = 'ordinary';   % constant trend
Metaopts.Corr.Type = 'Separable';
Metaopts.Corr.Family = 'Gaussian';
Metaopts.EstimMethod = 'ML';
Metaopts.Optim.Method = 'HGA';
Metaopts.Optim.MaxIter = 500;
Metaopts.Optim.HGA.nPop = 50;

end

