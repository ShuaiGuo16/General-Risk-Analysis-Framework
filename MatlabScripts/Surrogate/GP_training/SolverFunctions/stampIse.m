function [indi,indj,D,b] = stampIse(nodes, X, Y, distances, areas, ind, boundary,dim)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

lambda=boundary.lambda; 

% Stecil 

% East 
D3=((pi.*((nodes.nodesH_X(2,2)-nodes.nodesC_X(2,2)).*((X(2,2)-X(2,3))./2 + (3.*(X(2,3)-nodes.nodesV_X(2,3)))./4 + distances.dX_nodesV_H(2:end,2:end)./4) + (nodes.nodesH_Y(2,2)-nodes.nodesC_Y(2,2)).*((Y(2,2)-Y(2,3))./2 + (3.*(Y(2,3)-nodes.nodesV_Y(2,3)))./4 + distances.dY_nodesV_H(2:end,2:end)./4)).*((nodes.nodesH_X(2,2)-nodes.nodesC_X(2,2)) + 2.*nodes.nodesC_X(2,2)).*((3.*lambda(2:end-1,3:end))./8 + (3.*lambda(2:end-1,2:end-1))./8 + lambda(3:end,2:end-1)./8 + lambda(3,3)./8).^2)./areas.S_Ise_right + (pi.*(distances.dX_nodesC_H(2:end,:) + 2.*nodes.nodesC_X(2,1)).*((-distances.dX_nodesH_V(2:end,2:end).*distances.dX_nodesC_H(2:end,:))./4 + (-distances.dY_nodesH_V(2:end,2:end).*distances.dY_nodesC_H(2:end,:))./4).*(lambda(2:end-1,2:end-1)./2 + lambda(3:end,2:end-1)./2).^2)./areas.S_V(2:end,:))./areas.S_P_Ise; 

% West 
D_3=((pi.*((nodes.nodesC_X(1,1)-nodes.nodesV_X(1,2)) + 2.*nodes.nodesV_X(1,2)).*((nodes.nodesC_X(1,1)-nodes.nodesV_X(1,2)).*((X(2,2)-nodes.nodesH_X(2,1))./4 + distances.dX_nodesH_V(1:end-1,1:end-1)./4) + (nodes.nodesC_Y(1,1)-nodes.nodesV_Y(1,2)).*((Y(2,2)-nodes.nodesH_Y(2,1))./4 + distances.dY_nodesH_V(1:end-1,1:end-1)./4)).*((3.*lambda(1:end-2,2:end-1))./8 + (3.*lambda(2:end-1,2:end-1))./8 + lambda(2:end-1,1:end-2)./8 + lambda(1,1)./8).^2)./areas.S_Ise_left + (pi.*(distances.dX_nodesC_V(:,1:end-1).*(-distances.dX_nodesV_H(1:end-1,1:end-1)./4 + distances.dX_nodesV_H(2:end,1:end-1)./4 + distances.dX_nodesV_V(:,1:end-2)) + distances.dY_nodesC_V(:,1:end-1).*(-distances.dY_nodesV_H(1:end-1,1:end-1)./4 + distances.dY_nodesV_H(2:end,1:end-1)./4 + distances.dY_nodesV_V(:,1:end-2))).*(distances.dX_nodesC_V(:,1:end-1) + 2.*nodes.nodesC_X(1,1)).*(lambda(2:end-1,2:end-1)./2 + lambda(2:end-1,1:end-2)./2).^2)./areas.S_H(:,1:end-1) + (pi.*(distances.dX_nodesC_H(2:end,:) + 2.*nodes.nodesC_X(2,1)).*((distances.dX_nodesH_V(2:end,1:end-1).*distances.dX_nodesC_H(2:end,:))./4 + (distances.dY_nodesH_V(2:end,1:end-1).*distances.dY_nodesC_H(2:end,:))./4).*(lambda(2:end-1,2:end-1)./2 + lambda(3:end,2:end-1)./2).^2)./areas.S_V(2:end,:))./areas.S_P_Ise; 

% South 
D1=((pi.*((nodes.nodesH_X(2,2)-nodes.nodesC_X(2,2)) + 2.*nodes.nodesC_X(2,2)).*((nodes.nodesH_X(2,2)-nodes.nodesC_X(2,2)).*((nodes.nodesV_X(2,2)-X(2,2))./4 + distances.dX_nodesV_H(2:end,2:end)./4) + (nodes.nodesH_Y(2,2)-nodes.nodesC_Y(2,2)).*((nodes.nodesV_Y(2,2)-Y(2,2))./4 + distances.dY_nodesV_H(2:end,2:end)./4)).*((3.*lambda(2:end-1,3:end))./8 + (3.*lambda(2:end-1,2:end-1))./8 + lambda(3:end,2:end-1)./8 + lambda(3,3)./8).^2)./areas.S_Ise_right + (pi.*(distances.dX_nodesC_H(2:end,:).*(-distances.dX_nodesH_V(2:end,2:end)./4 + distances.dX_nodesH_V(2:end,1:end-1)./4 + distances.dX_nodesH_H(3:end,:)) + distances.dY_nodesC_H(2:end,:).*(-distances.dY_nodesH_V(2:end,2:end)./4 + distances.dY_nodesH_V(2:end,1:end-1)./4 + distances.dY_nodesH_H(3:end,:))).*(distances.dX_nodesC_H(2:end,:) + 2.*nodes.nodesC_X(2,1)).*(lambda(2:end-1,2:end-1)./2 + lambda(3:end,2:end-1)./2).^2)./areas.S_V(2:end,:) + (pi.*(distances.dX_nodesC_V(:,1:end-1) + 2.*nodes.nodesC_X(1,1)).*((distances.dX_nodesV_H(2:end,1:end-1).*distances.dX_nodesC_V(:,1:end-1))./4 + (distances.dY_nodesV_H(2:end,1:end-1).*distances.dY_nodesC_V(:,1:end-1))./4).*(lambda(2:end-1,2:end-1)./2 + lambda(2:end-1,1:end-2)./2).^2)./areas.S_H(:,1:end-1))./areas.S_P_Ise; 

% North 
D_1=((pi.*((nodes.nodesC_X(1,1)-nodes.nodesV_X(1,2)).*((X(1,2)-X(2,2))./2 + (3.*(nodes.nodesH_X(1,1)-X(1,2)))./4 + distances.dX_nodesH_V(1:end-1,1:end-1)./4) + (nodes.nodesC_Y(1,1)-nodes.nodesV_Y(1,2)).*((Y(1,2)-Y(2,2))./2 + (3.*(nodes.nodesH_Y(1,1)-Y(1,2)))./4 + distances.dY_nodesH_V(1:end-1,1:end-1)./4)).*((nodes.nodesC_X(1,1)-nodes.nodesV_X(1,2)) + 2.*nodes.nodesV_X(1,2)).*((3.*lambda(1:end-2,2:end-1))./8 + (3.*lambda(2:end-1,2:end-1))./8 + lambda(2:end-1,1:end-2)./8 + lambda(1,1)./8).^2)./areas.S_Ise_left + (pi.*(distances.dX_nodesC_V(:,1:end-1) + 2.*nodes.nodesC_X(1,1)).*((-distances.dX_nodesV_H(1:end-1,1:end-1).*distances.dX_nodesC_V(:,1:end-1))./4 + (-distances.dY_nodesV_H(1:end-1,1:end-1).*distances.dY_nodesC_V(:,1:end-1))./4).*(lambda(2:end-1,2:end-1)./2 + lambda(2:end-1,1:end-2)./2).^2)./areas.S_H(:,1:end-1))./areas.S_P_Ise; 

% NW 
D_4=((pi.*((nodes.nodesC_X(1,1)-nodes.nodesV_X(1,2)) + 2.*nodes.nodesV_X(1,2)).*((nodes.nodesC_X(1,1)-nodes.nodesV_X(1,2)).*((nodes.nodesH_X(1,1)-X(1,2))./4 + distances.dX_nodesH_V(1:end-1,1:end-1)./4) + (nodes.nodesC_Y(1,1)-nodes.nodesV_Y(1,2)).*((nodes.nodesH_Y(1,1)-Y(1,2))./4 + distances.dY_nodesH_V(1:end-1,1:end-1)./4)).*((3.*lambda(1:end-2,2:end-1))./8 + (3.*lambda(2:end-1,2:end-1))./8 + lambda(2:end-1,1:end-2)./8 + lambda(1,1)./8).^2)./areas.S_Ise_left + (pi.*(distances.dX_nodesC_V(:,1:end-1) + 2.*nodes.nodesC_X(1,1)).*((-distances.dX_nodesV_H(1:end-1,1:end-1).*distances.dX_nodesC_V(:,1:end-1))./4 + (-distances.dY_nodesV_H(1:end-1,1:end-1).*distances.dY_nodesC_V(:,1:end-1))./4).*(lambda(2:end-1,2:end-1)./2 + lambda(2:end-1,1:end-2)./2).^2)./areas.S_H(:,1:end-1))./areas.S_P_Ise; 

% SW 
D_2=((pi.*(distances.dX_nodesC_H(2:end,:) + 2.*nodes.nodesC_X(2,1)).*((distances.dX_nodesH_V(2:end,1:end-1).*distances.dX_nodesC_H(2:end,:))./4 + (distances.dY_nodesH_V(2:end,1:end-1).*distances.dY_nodesC_H(2:end,:))./4).*(lambda(2:end-1,2:end-1)./2 + lambda(3:end,2:end-1)./2).^2)./areas.S_V(2:end,:) + (pi.*(distances.dX_nodesC_V(:,1:end-1) + 2.*nodes.nodesC_X(1,1)).*((distances.dX_nodesV_H(2:end,1:end-1).*distances.dX_nodesC_V(:,1:end-1))./4 + (distances.dY_nodesV_H(2:end,1:end-1).*distances.dY_nodesC_V(:,1:end-1))./4).*(lambda(2:end-1,2:end-1)./2 + lambda(2:end-1,1:end-2)./2).^2)./areas.S_H(:,1:end-1))./areas.S_P_Ise; 

% SE 
D4=((pi.*((nodes.nodesH_X(2,2)-nodes.nodesC_X(2,2)) + 2.*nodes.nodesC_X(2,2)).*((nodes.nodesH_X(2,2)-nodes.nodesC_X(2,2)).*((X(2,3)-nodes.nodesV_X(2,3))./4 + distances.dX_nodesV_H(2:end,2:end)./4) + (nodes.nodesH_Y(2,2)-nodes.nodesC_Y(2,2)).*((Y(2,3)-nodes.nodesV_Y(2,3))./4 + distances.dY_nodesV_H(2:end,2:end)./4)).*((3.*lambda(2:end-1,3:end))./8 + (3.*lambda(2:end-1,2:end-1))./8 + lambda(3:end,2:end-1)./8 + lambda(3,3)./8).^2)./areas.S_Ise_right + (pi.*(distances.dX_nodesC_H(2:end,:) + 2.*nodes.nodesC_X(2,1)).*((-distances.dX_nodesH_V(2:end,2:end).*distances.dX_nodesC_H(2:end,:))./4 + (-distances.dY_nodesH_V(2:end,2:end).*distances.dY_nodesC_H(2:end,:))./4).*(lambda(2:end-1,2:end-1)./2 + lambda(3:end,2:end-1)./2).^2)./areas.S_V(2:end,:))./areas.S_P_Ise; 

% P 
D0=((pi.*((nodes.nodesC_X(1,1)-nodes.nodesV_X(1,2)).*((X(1,2)-X(2,2))./2 + (3.*(X(2,2)-nodes.nodesH_X(2,1)))./4 + distances.dX_nodesH_V(1:end-1,1:end-1)./4) + (nodes.nodesC_Y(1,1)-nodes.nodesV_Y(1,2)).*((Y(1,2)-Y(2,2))./2 + (3.*(Y(2,2)-nodes.nodesH_Y(2,1)))./4 + distances.dY_nodesH_V(1:end-1,1:end-1)./4)).*((nodes.nodesC_X(1,1)-nodes.nodesV_X(1,2)) + 2.*nodes.nodesV_X(1,2)).*((3.*lambda(1:end-2,2:end-1))./8 + (3.*lambda(2:end-1,2:end-1))./8 + lambda(2:end-1,1:end-2)./8 + lambda(1,1)./8).^2)./areas.S_Ise_left + (pi.*((nodes.nodesH_X(2,2)-nodes.nodesC_X(2,2)).*((X(2,2)-X(2,3))./2 + (3.*(nodes.nodesV_X(2,2)-X(2,2)))./4 + distances.dX_nodesV_H(2:end,2:end)./4) + (nodes.nodesH_Y(2,2)-nodes.nodesC_Y(2,2)).*((Y(2,2)-Y(2,3))./2 + (3.*(nodes.nodesV_Y(2,2)-Y(2,2)))./4 + distances.dY_nodesV_H(2:end,2:end)./4)).*((nodes.nodesH_X(2,2)-nodes.nodesC_X(2,2)) + 2.*nodes.nodesC_X(2,2)).*((3.*lambda(2:end-1,3:end))./8 + (3.*lambda(2:end-1,2:end-1))./8 + lambda(3:end,2:end-1)./8 + lambda(3,3)./8).^2)./areas.S_Ise_right + (pi.*(distances.dX_nodesC_H(2:end,:).*(-distances.dX_nodesH_H(2:end-1,:) + -distances.dX_nodesH_V(2:end,2:end)./4 + distances.dX_nodesH_V(2:end,1:end-1)./4) + distances.dY_nodesC_H(2:end,:).*(-distances.dY_nodesH_H(2:end-1,:) + -distances.dY_nodesH_V(2:end,2:end)./4 + distances.dY_nodesH_V(2:end,1:end-1)./4)).*(distances.dX_nodesC_H(2:end,:) + 2.*nodes.nodesC_X(2,1)).*(lambda(2:end-1,2:end-1)./2 + lambda(3:end,2:end-1)./2).^2)./areas.S_V(2:end,:) + (pi.*(distances.dX_nodesC_V(:,1:end-1).*( -distances.dX_nodesV_V(:,2:end-1) + -distances.dX_nodesV_H(1:end-1,1:end-1)./4 + distances.dX_nodesV_H(2:end,1:end-1)./4) + distances.dY_nodesC_V(:,1:end-1).*(-distances.dY_nodesV_V(:,2:end-1) + -distances.dY_nodesV_H(1:end-1,1:end-1)./4 + distances.dY_nodesV_H(2:end,1:end-1)./4)).*(distances.dX_nodesC_V(:,1:end-1) + 2.*nodes.nodesC_X(1,1)).*(lambda(2:end-1,2:end-1)./2 + lambda(2:end-1,1:end-2)./2).^2)./areas.S_H(:,1:end-1))./areas.S_P_Ise; 


ind0 = ind(2,2);
ind1 = ind(3,2);
ind_1 = ind(1,2);
ind_2 = ind(3,1);
ind3 = ind(2,3);
ind_3 = ind(2,1);
ind4 = ind(3,3);
ind_4 = ind(1,1);


D =    [D0(:)   ; D1(:)   ; D_1(:)   ; D_2(:)   ; D3(:)   ; D_3(:)   ; D4(:)   ; D_4(:)  ];
indi = [ind0(:) ; ind0(:) ; ind0(:)  ; ind0(:)  ; ind0(:) ; ind0(:)  ; ind0(:) ; ind0(:) ];
indj = [ind0(:) ; ind1(:) ; ind_1(:) ; ind_2(:) ; ind3(:) ; ind_3(:) ; ind4(:) ; ind_4(:)];

a= D0 + D1 + D_1 +D_2 +D3 +D_3 +D4 +D_4;

end
