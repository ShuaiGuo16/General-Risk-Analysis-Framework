function [indi,indj,D]= stampInner(nodesC_X, dimY, dimX, distances, areas, ind, boundary, dim)
%stampInner calculates the linear equation for the inner node (i, j)

% Calculate some help values
% Nomencature:
%
%    NW(i-1,j-1)   Nw - - (i-1,j) - -Ne     NE(i-1,j+1)
%
%                 |                 |
%
%       nW - - - - nw ------ n ------ ne - - - nE
%                 |                 |
%       |         |        |        |       |
%                 |                 |
%   W(i, j-1) - - w - - P (i,j) - - e - -  E (i,j+1)
%                 |                 |
%       |         |        |        |       |
%                 |                 |
%      sW - - - - sw ------ s ------ se - - - sE
%
%                 |                 |
%
%   SW(i+1,j-1)   Sw - - (i+1,j) - - Se      SE(i+1,j+1)
%
%

lambda=boundary.lambda; 

% Stecil 

% East 
D3=((pi.*(-distances.dX_nodesC_V(:,2:end).*(-distances.dX_nodesV_H(1:end-1,2:end)./4 + distances.dX_nodesV_H(2:end,2:end)./4 + -distances.dX_nodesV_V(:,3:end)) + -distances.dY_nodesC_V(:,2:end).*(-distances.dY_nodesV_H(1:end-1,2:end)./4 + distances.dY_nodesV_H(2:end,2:end)./4 + -distances.dY_nodesV_V(:,3:end))).*(-distances.dX_nodesC_V(:,2:end) + 2.*nodesC_X(2:end,2:end)).*(lambda(2:end-1,3:end)./2 + lambda(2:end-1,2:end-1)./2).^2)./areas.S_H(:,2:end) + (pi.*(-distances.dX_nodesC_H(1:end-1,:) + 2.*nodesC_X(1:end-1,2:end)).*((-distances.dX_nodesH_V(1:end-1,2:end).*-distances.dX_nodesC_H(1:end-1,:))./4 + (-distances.dY_nodesH_V(1:end-1,2:end).*-distances.dY_nodesC_H(1:end-1,:))./4).*(lambda(1:end-2,2:end-1)./2 + lambda(2:end-1,2:end-1)./2).^2)./areas.S_V(1:end-1,:) + (pi.*(distances.dX_nodesC_H(2:end,:) + 2.*nodesC_X(2:end,1:end-1)).*((-distances.dX_nodesH_V(2:end,2:end).*distances.dX_nodesC_H(2:end,:))./4 + (-distances.dY_nodesH_V(2:end,2:end).*distances.dY_nodesC_H(2:end,:))./4).*(lambda(2:end-1,2:end-1)./2 + lambda(3:end,2:end-1)./2).^2)./areas.S_V(2:end,:))./areas.S_C; 

% West 
D_3=((pi.*(distances.dX_nodesC_V(:,1:end-1).*(-distances.dX_nodesV_H(1:end-1,1:end-1)./4 + distances.dX_nodesV_H(2:end,1:end-1)./4 + distances.dX_nodesV_V(:,1:end-2)) + distances.dY_nodesC_V(:,1:end-1).*(-distances.dY_nodesV_H(1:end-1,1:end-1)./4 + distances.dY_nodesV_H(2:end,1:end-1)./4 + distances.dY_nodesV_V(:,1:end-2))).*(distances.dX_nodesC_V(:,1:end-1) + 2.*nodesC_X(1:end-1,1:end-1)).*(lambda(2:end-1,2:end-1)./2 + lambda(2:end-1,1:end-2)./2).^2)./areas.S_H(:,1:end-1) + (pi.*(-distances.dX_nodesC_H(1:end-1,:) + 2.*nodesC_X(1:end-1,2:end)).*((distances.dX_nodesH_V(1:end-1,1:end-1).*-distances.dX_nodesC_H(1:end-1,:))./4 + (distances.dY_nodesH_V(1:end-1,1:end-1).*-distances.dY_nodesC_H(1:end-1,:))./4).*(lambda(1:end-2,2:end-1)./2 + lambda(2:end-1,2:end-1)./2).^2)./areas.S_V(1:end-1,:) + (pi.*(distances.dX_nodesC_H(2:end,:) + 2.*nodesC_X(2:end,1:end-1)).*((distances.dX_nodesH_V(2:end,1:end-1).*distances.dX_nodesC_H(2:end,:))./4 + (distances.dY_nodesH_V(2:end,1:end-1).*distances.dY_nodesC_H(2:end,:))./4).*(lambda(2:end-1,2:end-1)./2 + lambda(3:end,2:end-1)./2).^2)./areas.S_V(2:end,:))./areas.S_C; 

% South 
D1=((pi.*(distances.dX_nodesC_H(2:end,:).*(-distances.dX_nodesH_V(2:end,2:end)./4 + distances.dX_nodesH_V(2:end,1:end-1)./4 + distances.dX_nodesH_H(3:end,:)) + distances.dY_nodesC_H(2:end,:).*(-distances.dY_nodesH_V(2:end,2:end)./4 + distances.dY_nodesH_V(2:end,1:end-1)./4 + distances.dY_nodesH_H(3:end,:))).*(distances.dX_nodesC_H(2:end,:) + 2.*nodesC_X(2:end,1:end-1)).*(lambda(2:end-1,2:end-1)./2 + lambda(3:end,2:end-1)./2).^2)./areas.S_V(2:end,:) + (pi.*(-distances.dX_nodesC_V(:,2:end) + 2.*nodesC_X(2:end,2:end)).*((distances.dX_nodesV_H(2:end,2:end).*-distances.dX_nodesC_V(:,2:end))./4 + (distances.dY_nodesV_H(2:end,2:end).*-distances.dY_nodesC_V(:,2:end))./4).*(lambda(2:end-1,3:end)./2 + lambda(2:end-1,2:end-1)./2).^2)./areas.S_H(:,2:end) + (pi.*(distances.dX_nodesC_V(:,1:end-1) + 2.*nodesC_X(1:end-1,1:end-1)).*((distances.dX_nodesV_H(2:end,1:end-1).*distances.dX_nodesC_V(:,1:end-1))./4 + (distances.dY_nodesV_H(2:end,1:end-1).*distances.dY_nodesC_V(:,1:end-1))./4).*(lambda(2:end-1,2:end-1)./2 + lambda(2:end-1,1:end-2)./2).^2)./areas.S_H(:,1:end-1))./areas.S_C; 

% North 
D_1=((pi.*(-distances.dX_nodesC_H(1:end-1,:).*(-distances.dX_nodesH_V(1:end-1,2:end)./4 + distances.dX_nodesH_V(1:end-1,1:end-1)./4 + -distances.dX_nodesH_H(1:end-2,:)) + -distances.dY_nodesC_H(1:end-1,:).*(-distances.dY_nodesH_V(1:end-1,2:end)./4 + distances.dY_nodesH_V(1:end-1,1:end-1)./4 + -distances.dY_nodesH_H(1:end-2,:))).*(-distances.dX_nodesC_H(1:end-1,:) + 2.*nodesC_X(1:end-1,2:end)).*(lambda(1:end-2,2:end-1)./2 + lambda(2:end-1,2:end-1)./2).^2)./areas.S_V(1:end-1,:) + (pi.*(-distances.dX_nodesC_V(:,2:end) + 2.*nodesC_X(2:end,2:end)).*((-distances.dX_nodesV_H(1:end-1,2:end).*-distances.dX_nodesC_V(:,2:end))./4 + (-distances.dY_nodesV_H(1:end-1,2:end).*-distances.dY_nodesC_V(:,2:end))./4).*(lambda(2:end-1,3:end)./2 + lambda(2:end-1,2:end-1)./2).^2)./areas.S_H(:,2:end) + (pi.*(distances.dX_nodesC_V(:,1:end-1) + 2.*nodesC_X(1:end-1,1:end-1)).*((-distances.dX_nodesV_H(1:end-1,1:end-1).*distances.dX_nodesC_V(:,1:end-1))./4 + (-distances.dY_nodesV_H(1:end-1,1:end-1).*distances.dY_nodesC_V(:,1:end-1))./4).*(lambda(2:end-1,2:end-1)./2 + lambda(2:end-1,1:end-2)./2).^2)./areas.S_H(:,1:end-1))./areas.S_C; 

% NW 
D_4=((pi.*(-distances.dX_nodesC_H(1:end-1,:) + 2.*nodesC_X(1:end-1,2:end)).*((distances.dX_nodesH_V(1:end-1,1:end-1).*-distances.dX_nodesC_H(1:end-1,:))./4 + (distances.dY_nodesH_V(1:end-1,1:end-1).*-distances.dY_nodesC_H(1:end-1,:))./4).*(lambda(1:end-2,2:end-1)./2 + lambda(2:end-1,2:end-1)./2).^2)./areas.S_V(1:end-1,:) + (pi.*(distances.dX_nodesC_V(:,1:end-1) + 2.*nodesC_X(1:end-1,1:end-1)).*((-distances.dX_nodesV_H(1:end-1,1:end-1).*distances.dX_nodesC_V(:,1:end-1))./4 + (-distances.dY_nodesV_H(1:end-1,1:end-1).*distances.dY_nodesC_V(:,1:end-1))./4).*(lambda(2:end-1,2:end-1)./2 + lambda(2:end-1,1:end-2)./2).^2)./areas.S_H(:,1:end-1))./areas.S_C; 

% NE 
D2=((pi.*(-distances.dX_nodesC_V(:,2:end) + 2.*nodesC_X(2:end,2:end)).*((-distances.dX_nodesV_H(1:end-1,2:end).*-distances.dX_nodesC_V(:,2:end))./4 + (-distances.dY_nodesV_H(1:end-1,2:end).*-distances.dY_nodesC_V(:,2:end))./4).*(lambda(2:end-1,3:end)./2 + lambda(2:end-1,2:end-1)./2).^2)./areas.S_H(:,2:end) + (pi.*(-distances.dX_nodesC_H(1:end-1,:) + 2.*nodesC_X(1:end-1,2:end)).*((-distances.dX_nodesH_V(1:end-1,2:end).*-distances.dX_nodesC_H(1:end-1,:))./4 + (-distances.dY_nodesH_V(1:end-1,2:end).*-distances.dY_nodesC_H(1:end-1,:))./4).*(lambda(1:end-2,2:end-1)./2 + lambda(2:end-1,2:end-1)./2).^2)./areas.S_V(1:end-1,:))./areas.S_C; 

% SW 
D_2=((pi.*(distances.dX_nodesC_H(2:end,:) + 2.*nodesC_X(2:end,1:end-1)).*((distances.dX_nodesH_V(2:end,1:end-1).*distances.dX_nodesC_H(2:end,:))./4 + (distances.dY_nodesH_V(2:end,1:end-1).*distances.dY_nodesC_H(2:end,:))./4).*(lambda(2:end-1,2:end-1)./2 + lambda(3:end,2:end-1)./2).^2)./areas.S_V(2:end,:) + (pi.*(distances.dX_nodesC_V(:,1:end-1) + 2.*nodesC_X(1:end-1,1:end-1)).*((distances.dX_nodesV_H(2:end,1:end-1).*distances.dX_nodesC_V(:,1:end-1))./4 + (distances.dY_nodesV_H(2:end,1:end-1).*distances.dY_nodesC_V(:,1:end-1))./4).*(lambda(2:end-1,2:end-1)./2 + lambda(2:end-1,1:end-2)./2).^2)./areas.S_H(:,1:end-1))./areas.S_C; 

% SE 
D4=((pi.*(-distances.dX_nodesC_V(:,2:end) + 2.*nodesC_X(2:end,2:end)).*((distances.dX_nodesV_H(2:end,2:end).*-distances.dX_nodesC_V(:,2:end))./4 + (distances.dY_nodesV_H(2:end,2:end).*-distances.dY_nodesC_V(:,2:end))./4).*(lambda(2:end-1,3:end)./2 + lambda(2:end-1,2:end-1)./2).^2)./areas.S_H(:,2:end) + (pi.*(distances.dX_nodesC_H(2:end,:) + 2.*nodesC_X(2:end,1:end-1)).*((-distances.dX_nodesH_V(2:end,2:end).*distances.dX_nodesC_H(2:end,:))./4 + (-distances.dY_nodesH_V(2:end,2:end).*distances.dY_nodesC_H(2:end,:))./4).*(lambda(2:end-1,2:end-1)./2 + lambda(3:end,2:end-1)./2).^2)./areas.S_V(2:end,:))./areas.S_C; 

% P 
D0=((pi.*(-distances.dX_nodesC_V(:,2:end).*(distances.dX_nodesV_V(:,2:end-1) + -distances.dX_nodesV_H(1:end-1,2:end)./4 + distances.dX_nodesV_H(2:end,2:end)./4) + -distances.dY_nodesC_V(:,2:end).*(distances.dY_nodesV_V(:,2:end-1) + -distances.dY_nodesV_H(1:end-1,2:end)./4 + distances.dY_nodesV_H(2:end,2:end)./4)).*(-distances.dX_nodesC_V(:,2:end) + 2.*nodesC_X(2:end,2:end)).*(lambda(2:end-1,3:end)./2 + lambda(2:end-1,2:end-1)./2).^2)./areas.S_H(:,2:end) + (pi.*(-distances.dX_nodesC_H(1:end-1,:).*(distances.dX_nodesH_H(2:end-1,:) + -distances.dX_nodesH_V(1:end-1,2:end)./4 + distances.dX_nodesH_V(1:end-1,1:end-1)./4) + -distances.dY_nodesC_H(1:end-1,:).*(distances.dY_nodesH_H(2:end-1,:) + -distances.dY_nodesH_V(1:end-1,2:end)./4 + distances.dY_nodesH_V(1:end-1,1:end-1)./4)).*(-distances.dX_nodesC_H(1:end-1,:) + 2.*nodesC_X(1:end-1,2:end)).*(lambda(1:end-2,2:end-1)./2 + lambda(2:end-1,2:end-1)./2).^2)./areas.S_V(1:end-1,:) + (pi.*(distances.dX_nodesC_H(2:end,:).*(-distances.dX_nodesH_H(2:end-1,:) + -distances.dX_nodesH_V(2:end,2:end)./4 + distances.dX_nodesH_V(2:end,1:end-1)./4) + distances.dY_nodesC_H(2:end,:).*(-distances.dY_nodesH_H(2:end-1,:) + -distances.dY_nodesH_V(2:end,2:end)./4 + distances.dY_nodesH_V(2:end,1:end-1)./4)).*(distances.dX_nodesC_H(2:end,:) + 2.*nodesC_X(2:end,1:end-1)).*(lambda(2:end-1,2:end-1)./2 + lambda(3:end,2:end-1)./2).^2)./areas.S_V(2:end,:) + (pi.*(distances.dX_nodesC_V(:,1:end-1).*( -distances.dX_nodesV_V(:,2:end-1) + -distances.dX_nodesV_H(1:end-1,1:end-1)./4 + distances.dX_nodesV_H(2:end,1:end-1)./4) + distances.dY_nodesC_V(:,1:end-1).*(-distances.dY_nodesV_V(:,2:end-1) + -distances.dY_nodesV_H(1:end-1,1:end-1)./4 + distances.dY_nodesV_H(2:end,1:end-1)./4)).*(distances.dX_nodesC_V(:,1:end-1) + 2.*nodesC_X(1:end-1,1:end-1)).*(lambda(2:end-1,2:end-1)./2 + lambda(2:end-1,1:end-2)./2).^2)./areas.S_H(:,1:end-1))./areas.S_C; 


ind0 = ind(2:dimY-1,2:dimX-1);
ind1 = ind(3:dimY,2:dimX-1);
ind_1 = ind(1:dimY-2,2:dimX-1);
ind2 = ind(1:dimY-2,3:dimX);
ind_2 = ind(3:dimY,1:dimX-2);
ind3 = ind(2:dimY-1,3:dimX);
ind_3 = ind(2:dimY-1,1:dimX-2);
ind4 = ind(3:dimY,3:dimX);
ind_4 = ind(1:dimY-2,1:dimX-2);


D =    [D0(:)   ; D1(:)   ; D_1(:)   ; D2(:)   ; D_2(:)   ; D3(:)   ; D_3(:)   ; D4(:)   ; D_4(:)  ];
indi = [ind0(:) ; ind0(:) ; ind0(:)  ; ind0(:) ; ind0(:)  ; ind0(:) ; ind0(:)  ; ind0(:) ; ind0(:) ];
indj = [ind0(:) ; ind1(:) ; ind_1(:) ; ind2(:) ; ind_2(:) ; ind3(:) ; ind_3(:) ; ind4(:) ; ind_4(:)];

a= D0 + D1 + D_1 +D2 +D_2 +D3 +D_3 + D4 +D_4;

end
