function [indi,indj,D,b] = stampW(nodes, dimY, dimX, distances, distances_west, areas, ind, boundary,dim)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

switch boundary.west
    case 'Dirichlet'
        indi = reshape(ind(2:dimY-1,1),(dimY-2),1);
        indj = reshape(ind(2:dimY-1,1),(dimY-2),1);
        D = ones(dimY-2,1);
        b = boundary.TD.west*ones(dimY-2,1);
        
    case {'Neumann'}
        
        % Stecil
        lambda=boundary.lambda;
        
        % Right 
        D1=((pi.*(distances_west.dX_nodesC_H(2:end).*(distances_west.dX_nodesH_V(2:end)./2 + (3.*distances_west.dX_nodesH_H(3:end))./4 + -distances.dX_nodesH_V(2:end,1)./4) + distances_west.dY_nodesC_H(2:end).*(distances_west.dY_nodesH_V(2:end)./2 + (3.*distances_west.dY_nodesH_H(3:end))./4 + -distances.dY_nodesH_V(2:end,1)./4)).*(distances_west.dX_nodesC_H(2:end) + 2.*nodes.nodesV_X(2:end,1)).*((3.*lambda(2:end-1,1))./8 + lambda(2:end-1,2)./8 + (3.*lambda(3:end,1))./8 + lambda(3:end,2)./8).^2)./areas.S_W(2:end) + (pi.*(-distances.dX_nodesC_V(:,1) + 2.*nodes.nodesC_X(2:end,1)).*((distances.dX_nodesV_H(2:end,1).*-distances.dX_nodesC_V(:,1))./4 + (distances.dY_nodesV_H(2:end,1).*-distances.dY_nodesC_V(:,1))./4).*(lambda(2:end-1,1)./2 + lambda(2:end-1,2)./2).^2)./areas.S_H(:,1))./areas.S_P_west; 

        % DomainRight 
        D4=((pi.*(distances_west.dX_nodesC_H(2:end) + 2.*nodes.nodesV_X(2:end,1)).*(distances_west.dX_nodesC_H(2:end).*(distances_west.dX_nodesH_H(3:end)./4 + -distances.dX_nodesH_V(2:end,1)./4) + distances_west.dY_nodesC_H(2:end).*(distances_west.dY_nodesH_H(3:end)./4 + -distances.dY_nodesH_V(2:end,1)./4)).*((3.*lambda(2:end-1,1))./8 + lambda(2:end-1,2)./8 + (3.*lambda(3:end,1))./8 + lambda(3:end,2)./8).^2)./areas.S_W(2:end) + (pi.*(-distances.dX_nodesC_V(:,1) + 2.*nodes.nodesC_X(2:end,1)).*((distances.dX_nodesV_H(2:end,1).*-distances.dX_nodesC_V(:,1))./4 + (distances.dY_nodesV_H(2:end,1).*-distances.dY_nodesC_V(:,1))./4).*(lambda(2:end-1,1)./2 + lambda(2:end-1,2)./2).^2)./areas.S_H(:,1))./areas.S_P_west; 

        % Domain 
        D3=((pi.*(-distances_west.dX_nodesC_H(1:end-1) + 2.*nodes.nodesC_X(1:end-1,1)).*(-distances_west.dX_nodesC_H(1:end-1).*(distances_west.dX_nodesH_H(2:end-1)./4 + -distances.dX_nodesH_V(1:end-1,1)./4) + -distances_west.dY_nodesC_H(1:end-1).*(distances_west.dY_nodesH_H(2:end-1)./4 + -distances.dY_nodesH_V(1:end-1,1)./4)).*((3.*lambda(2:end-1,1))./8 + lambda(2:end-1,2)./8 + (3.*lambda(1:end-2,1))./8 + lambda(1:end-2,2)./8).^2)./areas.S_W(1:end-1) + (pi.*(distances_west.dX_nodesC_H(2:end) + 2.*nodes.nodesV_X(2:end,1)).*(distances_west.dX_nodesC_H(2:end).*(-distances_west.dX_nodesH_H(2:end-1)./4 + -distances.dX_nodesH_V(2:end,1)./4) + distances_west.dY_nodesC_H(2:end).*(-distances_west.dY_nodesH_H(2:end-1)./4 + -distances.dY_nodesH_V(2:end,1)./4)).*((3.*lambda(2:end-1,1))./8 + lambda(2:end-1,2)./8 + (3.*lambda(3:end,1))./8 + lambda(3:end,2)./8).^2)./areas.S_W(2:end) + (pi.*(-distances.dX_nodesC_V(:,1).*(-distances.dX_nodesV_H(1:end-1,1)./4 + distances.dX_nodesV_H(2:end,1)./4 + -distances.dX_nodesV_V(:,2)) + -distances.dY_nodesC_V(:,1).*(-distances.dY_nodesV_H(1:end-1,1)./4 + distances.dY_nodesV_H(2:end,1)./4 + -distances.dY_nodesV_V(:,2))).*(-distances.dX_nodesC_V(:,1) + 2.*nodes.nodesC_X(2:end,1)).*(lambda(2:end-1,1)./2 + lambda(2:end-1,2)./2).^2)./areas.S_H(:,1))./areas.S_P_west; 

        % DomainLeft 
        D2=((pi.*(-distances_west.dX_nodesC_H(1:end-1) + 2.*nodes.nodesC_X(1:end-1,1)).*(-distances_west.dX_nodesC_H(1:end-1).*(-distances_west.dX_nodesH_H(1:end-2)./4 + -distances.dX_nodesH_V(1:end-1,1)./4) + -distances_west.dY_nodesC_H(1:end-1).*(-distances_west.dY_nodesH_H(1:end-2)./4 + -distances.dY_nodesH_V(1:end-1,1)./4)).*((3.*lambda(2:end-1,1))./8 + lambda(2:end-1,2)./8 + (3.*lambda(1:end-2,1))./8 + lambda(1:end-2,2)./8).^2)./areas.S_W(1:end-1) + (pi.*(-distances.dX_nodesC_V(:,1) + 2.*nodes.nodesC_X(2:end,1)).*((-distances.dX_nodesV_H(1:end-1,1).*-distances.dX_nodesC_V(:,1))./4 + (-distances.dY_nodesV_H(1:end-1,1).*-distances.dY_nodesC_V(:,1))./4).*(lambda(2:end-1,1)./2 + lambda(2:end-1,2)./2).^2)./areas.S_H(:,1))./areas.S_P_west; 

        % Left 
        D_1=((pi.*(-distances_west.dX_nodesC_H(1:end-1).*(distances_west.dX_nodesH_V(1:end-1)./2 + (3.*-distances_west.dX_nodesH_H(1:end-2))./4 + -distances.dX_nodesH_V(1:end-1,1)./4) + -distances_west.dY_nodesC_H(1:end-1).*(distances_west.dY_nodesH_V(1:end-1)./2 + (3.*-distances_west.dY_nodesH_H(1:end-2))./4 + -distances.dY_nodesH_V(1:end-1,1)./4)).*(-distances_west.dX_nodesC_H(1:end-1) + 2.*nodes.nodesC_X(1:end-1,1)).*((3.*lambda(2:end-1,1))./8 + lambda(2:end-1,2)./8 + (3.*lambda(1:end-2,1))./8 + lambda(1:end-2,2)./8).^2)./areas.S_W(1:end-1) + (pi.*(-distances.dX_nodesC_V(:,1) + 2.*nodes.nodesC_X(2:end,1)).*((-distances.dX_nodesV_H(1:end-1,1).*-distances.dX_nodesC_V(:,1))./4 + (-distances.dY_nodesV_H(1:end-1,1).*-distances.dY_nodesC_V(:,1))./4).*(lambda(2:end-1,1)./2 + lambda(2:end-1,2)./2).^2)./areas.S_H(:,1))./areas.S_P_west; 

        % Boundary 
        D0=((pi.*(-distances_west.dX_nodesC_H(1:end-1).*(distances_west.dX_nodesH_V(1:end-1)./2 + (3.*distances_west.dX_nodesH_H(2:end-1))./4 + -distances.dX_nodesH_V(1:end-1,1)./4) + -distances_west.dY_nodesC_H(1:end-1).*(distances_west.dY_nodesH_V(1:end-1)./2 + (3.*distances_west.dY_nodesH_H(2:end-1))./4 + -distances.dY_nodesH_V(1:end-1,1)./4)).*(-distances_west.dX_nodesC_H(1:end-1) + 2.*nodes.nodesC_X(1:end-1,1)).*((3.*lambda(2:end-1,1))./8 + lambda(2:end-1,2)./8 + (3.*lambda(1:end-2,1))./8 + lambda(1:end-2,2)./8).^2)./areas.S_W(1:end-1) + (pi.*(distances_west.dX_nodesC_H(2:end).*(distances_west.dX_nodesH_V(2:end)./2 + (3.*-distances_west.dX_nodesH_H(2:end-1))./4 + -distances.dX_nodesH_V(2:end,1)./4) + distances_west.dY_nodesC_H(2:end).*(distances_west.dY_nodesH_V(2:end)./2 + (3.*-distances_west.dY_nodesH_H(2:end-1))./4 + -distances.dY_nodesH_V(2:end,1)./4)).*(distances_west.dX_nodesC_H(2:end) + 2.*nodes.nodesV_X(2:end,1)).*((3.*lambda(2:end-1,1))./8 + lambda(2:end-1,2)./8 + (3.*lambda(3:end,1))./8 + lambda(3:end,2)./8).^2)./areas.S_W(2:end) + (pi.*(-distances.dX_nodesC_V(:,1).*(distances.dX_nodesV_V(:,1) + -distances.dX_nodesV_H(1:end-1,1)./4 + distances.dX_nodesV_H(2:end,1)./4) + -distances.dY_nodesC_V(:,1).*(distances.dY_nodesV_V(:,1) + -distances.dY_nodesV_H(1:end-1,1)./4 + distances.dY_nodesV_H(2:end,1)./4)).*(-distances.dX_nodesC_V(:,1) + 2.*nodes.nodesC_X(2:end,1)).*(lambda(2:end-1,1)./2 + lambda(2:end-1,2)./2).^2)./areas.S_H(:,1))./areas.S_P_west; 

        
        % Define boundary lengths
        boundlength=(distances.dX_nodesV_V(:,1).^2+distances.dY_nodesV_V(:,1).^2).^0.5;
        b = boundary.qdot.west*boundlength;
        
        ind0 = ind(2:dimY-1,1);
        ind1 = ind(3:dimY,1);
        ind_1 = ind(1:dimY-2,1);
        ind2 = ind(1:dimY-2,2);
        ind3 = ind(2:dimY-1,2);
        ind4 = ind(3:dimY,2);
        
        D =    [D0(:)   ; D1(:)   ; D_1(:)   ; D2(:)   ; D3(:)   ;   D4(:) ];
        indi = [ind0(:) ; ind0(:) ; ind0(:)  ; ind0(:) ; ind0(:) ; ind0(:) ];
        indj = [ind0(:) ; ind1(:) ; ind_1(:) ; ind2(:) ; ind3(:) ; ind4(:) ];
        
        %a  =  D0(:) + D1(:) + D_1(:) + D2(:) +  D3(:) +  D4(:);
end

end





%
% % East
%         D3 = ...
%             (distances_west.dY_nodesC_H(2:dimY-1,1).*(-distances.dY_nodesH_V(2:dimY-1,1)-distances_west.dY_nodesH_H(2:dimY-1,1)) +...
%             distances_west.dX_nodesC_H(2:dimY-1,1).*(-distances.dX_nodesH_V(2:dimY-1,1)-distances_west.dX_nodesH_H(2:dimY-1,1)))./(4*areas.S_W(2:dimY-1,1)) +...
%             (-distances.dY_nodesC_V(1:dimY-2,1).*(-distances.dY_nodesV_V(1:dimY-2,2)+distances.dY_nodesV_H(2:dimY-1,1)/4-distances.dY_nodesV_H(1:dimY-2,1)/4) +...
%             -distances.dX_nodesC_V(1:dimY-2,1).*(-distances.dX_nodesV_V(1:dimY-2,2)+distances.dX_nodesV_H(2:dimY-1,1)/4-distances.dX_nodesV_H(1:dimY-2,1)/4))./areas.S_H(1:dimY-2,1) +...
%             (-distances_west.dY_nodesC_H(1:dimY-2,1).*(-distances.dY_nodesH_V(1:dimY-2,1)+distances_west.dY_nodesH_H(2:dimY-1,1)) +...
%             -distances_west.dX_nodesC_H(1:dimY-2,1).*(-distances.dX_nodesH_V(1:dimY-2,1)+distances_west.dX_nodesH_H(2:dimY-1,1)))./(4*areas.S_W(1:dimY-2,1));
%
%         D3 = D3./areas.S_P_west;
%
%         % South
%         D1 = ...
%             (distances_west.dY_nodesC_H(2:dimY-1,1).*(3*distances_west.dY_nodesH_H(3:dimY,1)-distances.dY_nodesH_V(2:dimY-1,1)+2*distances_west.dY_nodesH_V(2:dimY-1,1)) +...
%             distances_west.dX_nodesC_H(2:dimY-1,1).*(3*distances_west.dX_nodesH_H(3:dimY,1)-distances.dX_nodesH_V(2:dimY-1,1)+2*distances_west.dX_nodesH_V(2:dimY-1,1)))./(4*areas.S_W(2:dimY-1,1)) +...
%             (-distances.dY_nodesC_V(1:dimY-2,1).*distances.dY_nodesV_H(2:dimY-1,1) +...
%             -distances.dX_nodesC_V(1:dimY-2,1).*distances.dX_nodesV_H(2:dimY-1,1))./(4*areas.S_H(1:dimY-2,1));
%
%         D1 = D1./areas.S_P_west;
%
%         % North
%         D_1 = ...
%             (-distances.dY_nodesC_V(1:dimY-2,1).*(-distances.dY_nodesV_H(1:dimY-2,1)) +...
%             -distances.dX_nodesC_V(1:dimY-2,1).*(-distances.dX_nodesV_H(1:dimY-2,1)))./(4*areas.S_H(1:dimY-2,1)) +...
%             (-distances_west.dY_nodesC_H(1:dimY-2,1).*(-3*distances_west.dY_nodesH_H(1:dimY-2,1)-distances.dY_nodesH_V(1:dimY-2,1)+2*distances_west.dY_nodesH_V(1:dimY-2,1)) +...
%             -distances_west.dX_nodesC_H(1:dimY-2,1).*(-3*distances_west.dX_nodesH_H(1:dimY-2,1)-distances.dX_nodesH_V(1:dimY-2,1)+2*distances_west.dX_nodesH_V(1:dimY-2,1)))./(4*areas.S_W(1:dimY-2,1));
%
%         D_1 = D_1./areas.S_P_west;
%
%         % NE
%         D2 = ...
%             (-distances.dY_nodesC_V(1:dimY-2,1).*(-distances.dY_nodesV_H(1:dimY-2,1)) +...
%             -distances.dX_nodesC_V(1:dimY-2,1).*(-distances.dX_nodesV_H(1:dimY-2,1)))./(4*areas.S_H(1:dimY-2,1)) +...
%             (-distances_west.dY_nodesC_H(1:dimY-2,1).*(-distances.dY_nodesH_V(1:dimY-2,1)-distances_west.dY_nodesH_H(1:dimY-2,1)) +...
%             -distances_west.dX_nodesC_H(1:dimY-2,1).*(-distances.dX_nodesH_V(1:dimY-2,1)-distances_west.dX_nodesH_H(1:dimY-2,1)))./(4*areas.S_W(1:dimY-2,1));
%
%         D2 = D2./areas.S_P_west;
%
%         % SE
%         D4 = ...
%             (distances_west.dY_nodesC_H(2:dimY-1,1).*(-distances.dY_nodesH_V(2:dimY-1,1)+distances_west.dY_nodesH_H(3:dimY,1)) +...
%             distances_west.dX_nodesC_H(2:dimY-1,1).*(-distances.dX_nodesH_V(2:dimY-1,1)+distances_west.dX_nodesH_H(3:dimY,1)))./(4*areas.S_W(2:dimY-1,1)) +...
%             (-distances.dY_nodesC_V(1:dimY-2,1).*distances.dY_nodesV_H(2:dimY-1,1) +...
%             -distances.dX_nodesC_V(1:dimY-2,1).*distances.dX_nodesV_H(2:dimY-1,1))./(4*areas.S_H(1:dimY-2,1));
%
%         D4 = D4./areas.S_P_west;
%
%         % P
%         D0 = ...
%             (distances_west.dY_nodesC_H(2:dimY-1,1).*(-3*distances_west.dY_nodesH_H(2:dimY-1,1)-distances.dY_nodesH_V(2:dimY-1,1)+2*distances_west.dY_nodesH_V(2:dimY-1,1)) +...
%             distances_west.dX_nodesC_H(2:dimY-1,1).*(-3*distances_west.dX_nodesH_H(2:dimY-1,1)-distances.dX_nodesH_V(2:dimY-1,1)+2*distances_west.dX_nodesH_V(2:dimY-1,1)))./(4*areas.S_W(2:dimY-1,1)) +...
%             (-distances.dY_nodesC_V(1:dimY-2,1).*(distances.dY_nodesV_V(1:dimY-2,1)+distances.dY_nodesV_H(2:dimY-1,1)/4-distances.dY_nodesV_H(1:dimY-2,1)/4) +...
%             -distances.dX_nodesC_V(1:dimY-2,1).*(distances.dX_nodesV_V(1:dimY-2,1)+distances.dX_nodesV_H(2:dimY-1,1)/4-distances.dX_nodesV_H(1:dimY-2,1)/4))./areas.S_H(1:dimY-2,1) +...
%             (-distances_west.dY_nodesC_H(1:dimY-2,1).*(3*distances_west.dY_nodesH_H(2:dimY-1,1)-distances.dY_nodesH_V(1:dimY-2,1)+2*distances_west.dY_nodesH_V(1:dimY-2,1)) +...
%             -distances_west.dX_nodesC_H(1:dimY-2,1).*(3*distances_west.dX_nodesH_H(2:dimY-1,1)-distances.dX_nodesH_V(1:dimY-2,1)+2*distances_west.dX_nodesH_V(1:dimY-2,1)))./(4*areas.S_W(1:dimY-2,1));
%
%         D0 = D0./areas.S_P_west;



