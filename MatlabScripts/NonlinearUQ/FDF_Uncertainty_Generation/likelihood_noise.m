function [NegLnLike,Psi,U]=likelihood_noise(x,ModelInfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE
%    ===> Calculates the negative of the concentrated ln-likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%	x - vector of log(theta) & SigmaSqr parameters
%	ModelInfo.X - n x k matrix of sample locations
%	ModelInfo.y - n x 1 vector of observed data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%	NegLnLike - concentrated log-likelihood *-1 for minimising
%	Psi - correlation matrix
%	U - Choleski factorisation of correlation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2007 A I J Forrester
%
% This program is free software: you can redistribute it and/or modify  it
% under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or any
% later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
% General Public License for more details.
% 
% You should have received a copy of the GNU General Public License and GNU
% Lesser General Public License along with this program. If not, see
% <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=ModelInfo.X;
y=ModelInfo.y;
theta=10.^x(1:2);
SigmaSqr = x(3);
p=2;  
n=size(X,1);

one=ones(n,1);   

% Pre-allocate memory
Psi=zeros(n,n);
% Build upper half of correlation matrix
for i=1:n
	for j=i+1:n
		Psi(i,j)=exp(-sum(theta.*abs(X(i,:)-X(j,:)).^p)); % abs added (February 10)
	end
end

% Add upper and lower halves and diagonal of ones plus 
% small number to reduce ill-conditioning
Psi = (Psi+Psi'+eye(n))*SigmaSqr+ModelInfo.noise;   % Upgrade Psi to covariance matrix

% Cholesky factorisation
[U,p]=chol(Psi);

% Use penalty if ill-conditioned
if p>0
    NegLnLike=1e4;
else
    
    % Sum lns of diagonal to find ln(abs(det(Psi)))
    LnDetPsi=2*sum(log(abs(diag(U))));

    % Use back-substitution of Cholesky instead of inverse
    mu=(one'*(U\(U'\y)))/(one'*(U\(U'\one)));
    NegLnLike = LnDetPsi + (y-mu*one)'*(U\(U'\(y-mu*one)));
end
