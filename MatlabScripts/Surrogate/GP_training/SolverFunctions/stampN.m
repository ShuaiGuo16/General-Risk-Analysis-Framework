function [indi,indj,D,b,Dx,indx] = stampN(nodes, dimY, dimX, distances, distances_north, areas, ind, boundary,dim)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Dx = [];
indx = [];

if strcmp(boundary.north,'Dirichlet')

        D = ones(dimX-2,1);
        indi = reshape(ind(1,2:dimX-1),(dimX-2),1);
        indj = reshape(ind(1,2:dimX-1),(dimX-2),1);
        b = boundary.TD.north*ones(dimX-2,1);
        
else
        % Stecil
        lambda=boundary.lambda;

        % Right 
        D_3=((pi.*(distances_north.dX_nodesC_V(1:end-1).*(-distances_north.dX_nodesV_H(1:end-1)./2 + (3.*distances_north.dX_nodesV_V(1:end-2))./4 + distances.dX_nodesV_H(1,1:end-1)./4) + distances_north.dY_nodesC_V(1:end-1).*(-distances_north.dY_nodesV_H(1:end-1)./2 + (3.*distances_north.dY_nodesV_V(1:end-2))./4 + distances.dY_nodesV_H(1,1:end-1)./4)).*(distances_north.dX_nodesC_V(1:end-1) + 2.*nodes.nodesH_X(1,1:end-1)).*((3.*lambda(1,2:end-1))./8 + lambda(2,2:end-1)./8 + (3.*lambda(1,1:end-2))./8 + lambda(2,1:end-2)./8).^2)./areas.S_N(1:end-1) + (pi.*(distances.dX_nodesC_H(1,:) + 2.*nodes.nodesC_X(1,1:end-1)).*((distances.dX_nodesH_V(1,1:end-1).*distances.dX_nodesC_H(1,:))./4 + (distances.dY_nodesH_V(1,1:end-1).*distances.dY_nodesC_H(1,:))./4).*(lambda(1,2:end-1)./2 + lambda(2,2:end-1)./2).^2)./areas.S_V(1,:))./areas.S_P_north; 

        % DomainRight 
        D_2=((pi.*(distances_north.dX_nodesC_V(1:end-1) + 2.*nodes.nodesH_X(1,1:end-1)).*(distances_north.dX_nodesC_V(1:end-1).*(distances_north.dX_nodesV_V(1:end-2)./4 + distances.dX_nodesV_H(1,1:end-1)./4) + distances_north.dY_nodesC_V(1:end-1).*(distances_north.dY_nodesV_V(1:end-2)./4 + distances.dY_nodesV_H(1,1:end-1)./4)).*((3.*lambda(1,2:end-1))./8 + lambda(2,2:end-1)./8 + (3.*lambda(1,1:end-2))./8 + lambda(2,1:end-2)./8).^2)./areas.S_N(1:end-1) + (pi.*(distances.dX_nodesC_H(1,:) + 2.*nodes.nodesC_X(1,1:end-1)).*((distances.dX_nodesH_V(1,1:end-1).*distances.dX_nodesC_H(1,:))./4 + (distances.dY_nodesH_V(1,1:end-1).*distances.dY_nodesC_H(1,:))./4).*(lambda(1,2:end-1)./2 + lambda(2,2:end-1)./2).^2)./areas.S_V(1,:))./areas.S_P_north; 

        % Domain 
        D1=((pi.*(-distances_north.dX_nodesC_V(2:end) + 2.*nodes.nodesC_X(1,2:end)).*(-distances_north.dX_nodesC_V(2:end).*(distances_north.dX_nodesV_V(2:end-1)./4 + distances.dX_nodesV_H(1,2:end)./4) + -distances_north.dY_nodesC_V(2:end).*(distances_north.dY_nodesV_V(2:end-1)./4 + distances.dY_nodesV_H(2,2:end)./4)).*((3.*lambda(1,2:end-1))./8 + lambda(2,2:end-1)./8 + (3.*lambda(1,3:end))./8 + lambda(2,3:end)./8).^2)./areas.S_N(2:end) + (pi.*(distances_north.dX_nodesC_V(1:end-1) + 2.*nodes.nodesH_X(1,1:end-1)).*(distances_north.dX_nodesC_V(1:end-1).*(-distances_north.dX_nodesV_V(2:end-1)./4 + distances.dX_nodesV_H(1,1:end-1)./4) + distances_north.dY_nodesC_V(1:end-1).*(-distances_north.dY_nodesV_V(2:end-1)./4 + distances.dY_nodesV_H(1,1:end-1)./4)).*((3.*lambda(1,2:end-1))./8 + lambda(2,2:end-1)./8 + (3.*lambda(1,1:end-2))./8 + lambda(2,1:end-2)./8).^2)./areas.S_N(1:end-1) + (pi.*(distances.dX_nodesC_H(1,:).*(-distances.dX_nodesH_V(1,2:end)./4 + distances.dX_nodesH_V(1,1:end-1)./4 + distances.dX_nodesH_H(2,:)) + distances.dY_nodesC_H(1,:).*(-distances.dY_nodesH_V(1,2:end)./4 + distances.dY_nodesH_V(1,1:end-1)./4 + distances.dY_nodesH_H(2,:))).*(distances.dX_nodesC_H(1,:) + 2.*nodes.nodesC_X(1,1:end-1)).*(lambda(1,2:end-1)./2 + lambda(2,2:end-1)./2).^2)./areas.S_V(1,:))./areas.S_P_north; 

        % DomainLeft 
        D4=((pi.*(-distances_north.dX_nodesC_V(2:end) + 2.*nodes.nodesC_X(1,2:end)).*(-distances_north.dX_nodesC_V(2:end).*(-distances_north.dX_nodesV_V(3:end)./4 + distances.dX_nodesV_H(1,2:end)./4) + -distances_north.dY_nodesC_V(2:end).*(-distances_north.dY_nodesV_V(3:end)./4 + distances.dY_nodesV_H(2,2:end)./4)).*((3.*lambda(1,2:end-1))./8 + lambda(2,2:end-1)./8 + (3.*lambda(1,3:end))./8 + lambda(2,3:end)./8).^2)./areas.S_N(2:end) + (pi.*(distances.dX_nodesC_H(1,:) + 2.*nodes.nodesC_X(1,1:end-1)).*((-distances.dX_nodesH_V(1,2:end).*distances.dX_nodesC_H(1,:))./4 + (-distances.dY_nodesH_V(1,2:end).*distances.dY_nodesC_H(1,:))./4).*(lambda(1,2:end-1)./2 + lambda(2,2:end-1)./2).^2)./areas.S_V(1,:))./areas.S_P_north; 

        % Left 
        D3=((pi.*(-distances_north.dX_nodesC_V(2:end).*(-distances_north.dX_nodesV_H(2:end)./2 + (3.*-distances_north.dX_nodesV_V(3:end))./4 + distances.dX_nodesV_H(1,2:end)./4) + -distances_north.dY_nodesC_V(2:end).*(-distances_north.dY_nodesV_H(2:end)./2 + (3.*-distances_north.dY_nodesV_V(3:end))./4 + distances.dY_nodesV_H(2,2:end)./4)).*(-distances_north.dX_nodesC_V(2:end) + 2.*nodes.nodesC_X(1,2:end)).*((3.*lambda(1,2:end-1))./8 + lambda(2,2:end-1)./8 + (3.*lambda(1,3:end))./8 + lambda(2,3:end)./8).^2)./areas.S_N(2:end) + (pi.*(distances.dX_nodesC_H(1,:) + 2.*nodes.nodesC_X(1,1:end-1)).*((-distances.dX_nodesH_V(1,2:end).*distances.dX_nodesC_H(1,:))./4 + (-distances.dY_nodesH_V(1,2:end).*distances.dY_nodesC_H(1,:))./4).*(lambda(1,2:end-1)./2 + lambda(2,2:end-1)./2).^2)./areas.S_V(1,:))./areas.S_P_north; 

        % Boundary 
        D0=((pi.*(-distances_north.dX_nodesC_V(2:end).*(-distances_north.dX_nodesV_H(2:end)./2 + (3.*distances_north.dX_nodesV_V(2:end-1))./4 + distances.dX_nodesV_H(1,2:end)./4) + -distances_north.dY_nodesC_V(2:end).*(-distances_north.dY_nodesV_H(2:end)./2 + (3.*distances_north.dY_nodesV_V(2:end-1))./4 + distances.dY_nodesV_H(2,2:end)./4)).*(-distances_north.dX_nodesC_V(2:end) + 2.*nodes.nodesC_X(1,2:end)).*((3.*lambda(1,2:end-1))./8 + lambda(2,2:end-1)./8 + (3.*lambda(1,3:end))./8 + lambda(2,3:end)./8).^2)./areas.S_N(2:end) + (pi.*(distances_north.dX_nodesC_V(1:end-1).*(-distances_north.dX_nodesV_H(1:end-1)./2 + (3.*-distances_north.dX_nodesV_V(2:end-1))./4 + distances.dX_nodesV_H(1,1:end-1)./4) + distances_north.dY_nodesC_V(1:end-1).*(-distances_north.dY_nodesV_H(1:end-1)./2 + (3.*-distances_north.dY_nodesV_V(2:end-1))./4 + distances.dY_nodesV_H(1,1:end-1)./4)).*(distances_north.dX_nodesC_V(1:end-1) + 2.*nodes.nodesH_X(1,1:end-1)).*((3.*lambda(1,2:end-1))./8 + lambda(2,2:end-1)./8 + (3.*lambda(1,1:end-2))./8 + lambda(2,1:end-2)./8).^2)./areas.S_N(1:end-1) + (pi.*(distances.dX_nodesC_H(1,:).*(-distances.dX_nodesH_H(1,:) + -distances.dX_nodesH_V(1,2:end)./4 + distances.dX_nodesH_V(1,1:end-1)./4) + distances.dY_nodesC_H(1,:).*(-distances.dY_nodesH_H(1,:) + -distances.dY_nodesH_V(1,2:end)./4 + distances.dY_nodesH_V(1,1:end-1)./4)).*(distances.dX_nodesC_H(1,:) + 2.*nodes.nodesC_X(1,1:end-1)).*(lambda(1,2:end-1)./2 + lambda(2,2:end-1)./2).^2)./areas.S_V(1,:))./areas.S_P_north; 

        % Define boundary lengths
        boundlength=pi*(2*nodes.nodesH_X(1,1:end-1)+distances.dX_nodesH_H(1,:)).*(distances.dX_nodesH_H(1,:).^2+distances.dY_nodesH_H(1,:).^2).^0.5;
        
        if strcmp(boundary.north,'Neumann')
            b = boundary.qdot.north*boundlength';
        elseif strcmp(boundary.north,'Robin')
            b = - boundary.alpha*boundary.Tinf*boundlength';
            
            if boundary.splitA
                Dx = - lambda(1,2:end-1).^2.*boundlength./areas.S_P_north;
                indx = ind(1,2:dimX-1)';
            else 
                D0 = D0 - boundary.alpha*lambda(1,2:end-1).^2.*boundlength./areas.S_P_north;
            end
            
        else
            error('unknown BC type')
        end
        
        ind0 = ind(1,2:dimX-1);
        ind1 = ind(2,2:dimX-1);
        ind_2 = ind(2,1:dimX-2);
        ind3 = ind(1,3:dimX);
        ind_3 = ind(1,1:dimX-2);
        ind4 = ind(2,3:dimX);
        
        D =    [D0(:)   ; D1(:)   ; D_2(:)   ; D3(:)   ; D_3(:)   ; D4(:)   ];
        indi = [ind0(:) ; ind0(:) ; ind0(:)  ; ind0(:) ; ind0(:)  ; ind0(:) ];
        indj = [ind0(:) ; ind1(:) ; ind_2(:) ; ind3(:) ; ind_3(:) ; ind4(:) ];
        
        
        
        
        %a  =  D0(:) + D1(:) + D_2(:) + D3(:) + D_3(:) + D4(:);
        
end
end



%
%
% % East
%         D3 = ...
%             (distances.dY_nodesC_H(1,1:dimX-2).*(-distances.dY_nodesH_V(1,2:dimX-1)) +...
%             distances.dX_nodesC_H(1,1:dimX-2).*(-distances.dX_nodesH_V(1,2:dimX-1)))./(4*areas.S_V(1,1:dimX-2)) +...
%             (-distances_north.dY_nodesC_V(1,2:dimX-1).*(-3*distances_north.dY_nodesV_V(1,3:dimX)+distances.dY_nodesV_H(1,2:dimX-1)-2*distances_north.dY_nodesV_H(1,2:dimX-1)) +...
%             -distances_north.dX_nodesC_V(1,2:dimX-1).*(-3*distances_north.dX_nodesV_V(1,3:dimX)+distances.dX_nodesV_H(1,2:dimX-1)-2*distances_north.dX_nodesV_H(1,2:dimX-1)))./(4*areas.S_N(1,2:dimX-1));
%
%         D3 = D3./areas.S_P_north;
%
%         % West
%         D_3 = ...
%             (distances.dY_nodesC_H(1,1:dimX-2).*(distances.dY_nodesH_V(1,1:dimX-2)) +...
%             distances.dX_nodesC_H(1,1:dimX-2).*(distances.dX_nodesH_V(1,1:dimX-2)))./(4*areas.S_V(1,1:dimX-2)) +...
%             (distances_north.dY_nodesC_V(1,1:dimX-2).*(3*distances_north.dY_nodesV_V(1,1:dimX-2)-2*distances_north.dY_nodesV_H(1,1:dimX-2)+distances.dY_nodesV_H(1,1:dimX-2)) +...
%             distances_north.dX_nodesC_V(1,1:dimX-2).*(3*distances_north.dX_nodesV_V(1,1:dimX-2)-2*distances_north.dX_nodesV_H(1,1:dimX-2)+distances.dX_nodesV_H(1,1:dimX-2)))./(4*areas.S_N(1,1:dimX-2));
%
%         D_3 = D_3./areas.S_P_north;
%
%         % South
%         D1 = ...
%             (distances.dY_nodesC_H(1,1:dimX-2).*(distances.dY_nodesH_H(2,1:dimX-2)-distances.dY_nodesH_V(1,2:dimX-1)/4+distances.dY_nodesH_V(1,1:dimX-2)/4) +...
%             distances.dX_nodesC_H(1,1:dimX-2).*(distances.dX_nodesH_H(2,1:dimX-2)-distances.dX_nodesH_V(1,2:dimX-1)/4+distances.dX_nodesH_V(1,1:dimX-2)/4))./areas.S_V(1,1:dimX-2) +...
%             (-distances_north.dY_nodesC_V(1,2:dimX-1).*(distances.dY_nodesV_H(1,2:dimX-1)+distances_north.dY_nodesV_V(1,2:dimX-1)) +...
%             -distances_north.dX_nodesC_V(1,2:dimX-1).*(distances.dX_nodesV_H(1,2:dimX-1)+distances_north.dX_nodesV_V(1,2:dimX-1)))./(4*areas.S_N(1,2:dimX-1)) +...
%             (distances_north.dY_nodesC_V(1,1:dimX-2).*(distances.dY_nodesV_H(1,1:dimX-2)-distances_north.dY_nodesV_V(1,2:dimX-1)) +...
%             distances_north.dX_nodesC_V(1,1:dimX-2).*(distances.dX_nodesV_H(1,1:dimX-2)-distances_north.dX_nodesV_V(1,2:dimX-1)))./(4*areas.S_N(1,1:dimX-2));
%
%         D1 = D1./areas.S_P_north;
%
%         % SW
%         D_2 = ...
%             (distances.dY_nodesC_H(1,1:dimX-2).*distances.dY_nodesH_V(1,1:dimX-2) +...
%             distances.dX_nodesC_H(1,1:dimX-2).*distances.dX_nodesH_V(1,1:dimX-2))./(4*areas.S_V(1,1:dimX-2)) +...
%             (distances_north.dY_nodesC_V(1,1:dimX-2).*(distances.dY_nodesV_H(1,1:dimX-2)+distances_north.dY_nodesV_V(1,1:dimX-2)) +...
%             distances_north.dX_nodesC_V(1,1:dimX-2).*(distances.dX_nodesV_H(1,1:dimX-2)+distances_north.dX_nodesV_V(1,1:dimX-2)))./(4*areas.S_N(1,1:dimX-2));
%
%         D_2 = D_2./areas.S_P_north;
%
%         % SE
%         D4 = ...
%             (distances.dY_nodesC_H(1,1:dimX-2).*(-distances.dY_nodesH_V(1,2:dimX-1)) +...
%             distances.dX_nodesC_H(1,1:dimX-2).*(-distances.dX_nodesH_V(1,2:dimX-1)))./(4*areas.S_V(1,1:dimX-2)) +...
%             (-distances_north.dY_nodesC_V(1,2:dimX-1).*(distances.dY_nodesV_H(1,2:dimX-1)-distances_north.dY_nodesV_V(1,3:dimX)) +...
%             -distances_north.dX_nodesC_V(1,2:dimX-1).*(distances.dX_nodesV_H(1,2:dimX-1)-distances_north.dX_nodesV_V(1,3:dimX)))./(4*areas.S_N(1,2:dimX-1));
%
%         D4 = D4./areas.S_P_north;
%
%         % P
%         D0 = ...
%             (distances.dY_nodesC_H(1,1:dimX-2).*(-distances.dY_nodesH_H(1,1:dimX-2)-distances.dY_nodesH_V(1,2:dimX-1)/4+distances.dY_nodesH_V(1,1:dimX-2)/4) +...
%             distances.dX_nodesC_H(1,1:dimX-2).*(-distances.dX_nodesH_H(1,1:dimX-2)-distances.dX_nodesH_V(1,2:dimX-1)/4+distances.dX_nodesH_V(1,1:dimX-2)/4))./areas.S_V(1,1:dimX-2) +...
%             (-distances_north.dY_nodesC_V(1,2:dimX-1).*(3*distances_north.dY_nodesV_V(1,2:dimX-1)+distances.dY_nodesV_H(1,2:dimX-1)-2*distances_north.dY_nodesV_H(1,2:dimX-1)) +...
%             -distances_north.dX_nodesC_V(1,2:dimX-1).*(3*distances_north.dX_nodesV_V(1,2:dimX-1)+distances.dX_nodesV_H(1,2:dimX-1)-2*distances_north.dX_nodesV_H(1,2:dimX-1)))./(4*areas.S_N(1,2:dimX-1)) +...
%             (distances_north.dY_nodesC_V(1,1:dimX-2).*(-3*distances_north.dY_nodesV_V(1,2:dimX-1)-2*distances_north.dY_nodesV_H(1,1:dimX-2)+distances.dY_nodesV_H(1,1:dimX-2)) +...
%             distances_north.dX_nodesC_V(1,1:dimX-2).*(-3*distances_north.dX_nodesV_V(1,2:dimX-1)-2*distances_north.dX_nodesV_H(1,1:dimX-2)+distances.dX_nodesV_H(1,1:dimX-2)))./(4*areas.S_N(1,1:dimX-2));
%
%         D0 = D0./areas.S_P_north;
