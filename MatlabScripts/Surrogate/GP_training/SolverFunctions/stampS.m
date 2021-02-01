function [indi,indj,D,b,Dx,indx] = stampS(nodes, dimY, dimX, distances, distances_south, areas, ind, boundary,dim)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Dx = [];
indx = [];

if strcmp(boundary.south,'Dirichlet')
        
        D = ones(dimX-2,1);
        indi = reshape(ind(dimY,2:dimX-1),(dimX-2),1);
        indj = reshape(ind(dimY,2:dimX-1),(dimX-2),1);
        b = boundary.TD.south*ones(dimX-2,1);        
        
else
        
        % Stecil
        lambda=boundary.lambda;
        
        % Right 
        D3=((pi.*(-distances_south.dX_nodesC_V(2:end).*(distances_south.dX_nodesV_H(2:end)./2 + (3.*-distances_south.dX_nodesV_V(3:end))./4 + -distances.dX_nodesV_H(end,2:end)./4) + -distances_south.dY_nodesC_V(2:end).*(distances_south.dY_nodesV_H(2:end)./2 + (3.*-distances_south.dY_nodesV_V(3:end))./4 + -distances.dY_nodesV_H(end,2:end)./4)).*(-distances_south.dX_nodesC_V(2:end) + 2.*nodes.nodesH_X(end,2:end)).*((3.*lambda(end,2:end-1))./8 + lambda(end-1,2:end-1)./8 + (3.*lambda(end,3:end))./8 + lambda(end-1,3:end)./8).^2)./areas.S_S(2:end) + (pi.*(-distances.dX_nodesC_H(end,:) + 2.*nodes.nodesC_X(end,2:end)).*((-distances.dX_nodesH_V(end,2:end).*-distances.dX_nodesC_H(end,:))./4 + (-distances.dY_nodesH_V(end,2:end).*-distances.dY_nodesC_H(end,:))./4).*(lambda(end,2:end-1)./2 + lambda(end-1,2:end-1)./2).^2)./areas.S_V(end,:))./areas.S_P_south; 

        % DomainRight 
        D2=((pi.*(-distances_south.dX_nodesC_V(2:end) + 2.*nodes.nodesH_X(end,2:end)).*(-distances_south.dX_nodesC_V(2:end).*(-distances_south.dX_nodesV_V(3:end)./4 + -distances.dX_nodesV_H(end,2:end)./4) + -distances_south.dY_nodesC_V(2:end).*(-distances_south.dY_nodesV_V(3:end)./4 + -distances.dY_nodesV_H(end,2:end)./4)).*((3.*lambda(end,2:end-1))./8 + lambda(end-1,2:end-1)./8 + (3.*lambda(end,3:end))./8 + lambda(end-1,3:end)./8).^2)./areas.S_S(2:end) + (pi.*(-distances.dX_nodesC_H(end,:) + 2.*nodes.nodesC_X(end,2:end)).*((-distances.dX_nodesH_V(end,2:end).*-distances.dX_nodesC_H(end,:))./4 + (-distances.dY_nodesH_V(end,2:end).*-distances.dY_nodesC_H(end,:))./4).*(lambda(end,2:end-1)./2 + lambda(end-1,2:end-1)./2).^2)./areas.S_V(end,:))./areas.S_P_south; 

        % Domain 
        D_1=((pi.*(distances_south.dX_nodesC_V(1:end-1) + 2.*nodes.nodesC_X(end,1:end-1)).*(distances_south.dX_nodesC_V(1:end-1).*(-distances_south.dX_nodesV_V(2:end-1)./4 + -distances.dX_nodesV_H(end,1:end-1)./4) + distances_south.dY_nodesC_V(1:end-1).*(-distances_south.dY_nodesV_V(2:end-1)./4 + -distances.dY_nodesV_H(end,1:end-1)./4)).*((3.*lambda(end,2:end-1))./8 + lambda(end-1,2:end-1)./8 + (3.*lambda(end,1:end-2))./8 + lambda(end-1,1:end-2)./8).^2)./areas.S_S(1:end-1) + (pi.*(-distances_south.dX_nodesC_V(2:end) + 2.*nodes.nodesH_X(end,2:end)).*(-distances_south.dX_nodesC_V(2:end).*(distances_south.dX_nodesV_V(2:end-1)./4 + -distances.dX_nodesV_H(end,2:end)./4) + -distances_south.dY_nodesC_V(2:end).*(distances_south.dY_nodesV_V(2:end-1)./4 + -distances.dY_nodesV_H(end,2:end)./4)).*((3.*lambda(end,2:end-1))./8 + lambda(end-1,2:end-1)./8 + (3.*lambda(end,3:end))./8 + lambda(end-1,3:end)./8).^2)./areas.S_S(2:end) + (pi.*(-distances.dX_nodesC_H(end,:).*(distances.dX_nodesH_V(end,1:end-1)./4 + -distances.dX_nodesH_V(end,2:end)./4 + -distances.dX_nodesH_H(end-1,:)) + -distances.dY_nodesC_H(end,:).*(distances.dY_nodesH_V(end,1:end-1)./4 + -distances.dY_nodesH_V(end,2:end)./4 + -distances.dY_nodesH_H(end-1,:))).*(-distances.dX_nodesC_H(end,:) + 2.*nodes.nodesC_X(end,2:end)).*(lambda(end,2:end-1)./2 + lambda(end-1,2:end-1)./2).^2)./areas.S_V(end,:))./areas.S_P_south; 

        % DomainLeft 
        D_4=((pi.*(distances_south.dX_nodesC_V(1:end-1) + 2.*nodes.nodesC_X(end,1:end-1)).*(distances_south.dX_nodesC_V(1:end-1).*(distances_south.dX_nodesV_V(1:end-2)./4 + -distances.dX_nodesV_H(end,1:end-1)./4) + distances_south.dY_nodesC_V(1:end-1).*(distances_south.dY_nodesV_V(1:end-2)./4 + -distances.dY_nodesV_H(end,1:end-1)./4)).*((3.*lambda(end,2:end-1))./8 + lambda(end-1,2:end-1)./8 + (3.*lambda(end,1:end-2))./8 + lambda(end-1,1:end-2)./8).^2)./areas.S_S(1:end-1) + (pi.*(-distances.dX_nodesC_H(end,:) + 2.*nodes.nodesC_X(end,2:end)).*((distances.dX_nodesH_V(end,1:end-1).*-distances.dX_nodesC_H(end,:))./4 + (distances.dY_nodesH_V(end,1:end-1).*-distances.dY_nodesC_H(end,:))./4).*(lambda(end,2:end-1)./2 + lambda(end-1,2:end-1)./2).^2)./areas.S_V(end,:))./areas.S_P_south; 

        % Left 
        D_3=((pi.*(distances_south.dX_nodesC_V(1:end-1).*(distances_south.dX_nodesV_H(1:end-1)./2 + (3.*distances_south.dX_nodesV_V(1:end-2))./4 + -distances.dX_nodesV_H(end,1:end-1)./4) + distances_south.dY_nodesC_V(1:end-1).*(distances_south.dY_nodesV_H(1:end-1)./2 + (3.*distances_south.dY_nodesV_V(1:end-2))./4 + -distances.dY_nodesV_H(end,1:end-1)./4)).*(distances_south.dX_nodesC_V(1:end-1) + 2.*nodes.nodesC_X(end,1:end-1)).*((3.*lambda(end,2:end-1))./8 + lambda(end-1,2:end-1)./8 + (3.*lambda(end,1:end-2))./8 + lambda(end-1,1:end-2)./8).^2)./areas.S_S(1:end-1) + (pi.*(-distances.dX_nodesC_H(end,:) + 2.*nodes.nodesC_X(end,2:end)).*((distances.dX_nodesH_V(end,1:end-1).*-distances.dX_nodesC_H(end,:))./4 + (distances.dY_nodesH_V(end,1:end-1).*-distances.dY_nodesC_H(end,:))./4).*(lambda(end,2:end-1)./2 + lambda(end-1,2:end-1)./2).^2)./areas.S_V(end,:))./areas.S_P_south; 

        % Boundary 
        D0=((pi.*(distances_south.dX_nodesC_V(1:end-1).*(distances_south.dX_nodesV_H(1:end-1)./2 + (3.*-distances_south.dX_nodesV_V(2:end-1))./4 + -distances.dX_nodesV_H(end,1:end-1)./4) + distances_south.dY_nodesC_V(1:end-1).*(distances_south.dY_nodesV_H(1:end-1)./2 + (3.*-distances_south.dY_nodesV_V(2:end-1))./4 + -distances.dY_nodesV_H(end,1:end-1)./4)).*(distances_south.dX_nodesC_V(1:end-1) + 2.*nodes.nodesC_X(end,1:end-1)).*((3.*lambda(end,2:end-1))./8 + lambda(end-1,2:end-1)./8 + (3.*lambda(end,1:end-2))./8 + lambda(end-1,1:end-2)./8).^2)./areas.S_S(1:end-1) + (pi.*(-distances_south.dX_nodesC_V(2:end).*(distances_south.dX_nodesV_H(2:end)./2 + (3.*distances_south.dX_nodesV_V(2:end-1))./4 + -distances.dX_nodesV_H(end,2:end)./4) + -distances_south.dY_nodesC_V(2:end).*(distances_south.dY_nodesV_H(2:end)./2 + (3.*distances_south.dY_nodesV_V(2:end-1))./4 + -distances.dY_nodesV_H(end,2:end)./4)).*(-distances_south.dX_nodesC_V(2:end) + 2.*nodes.nodesH_X(end,2:end)).*((3.*lambda(end,2:end-1))./8 + lambda(end-1,2:end-1)./8 + (3.*lambda(end,3:end))./8 + lambda(end-1,3:end)./8).^2)./areas.S_S(2:end) + (pi.*(-distances.dX_nodesC_H(end,:).*(distances.dX_nodesH_H(end,:) + distances.dX_nodesH_V(end,1:end-1)./4 + -distances.dX_nodesH_V(end,2:end)./4) + -distances.dY_nodesC_H(end,:).*(distances.dY_nodesH_H(end,:) + distances.dY_nodesH_V(end,1:end-1)./4 + -distances.dY_nodesH_V(end,2:end)./4)).*(-distances.dX_nodesC_H(end,:) + 2.*nodes.nodesC_X(end,2:end)).*(lambda(end,2:end-1)./2 + lambda(end-1,2:end-1)./2).^2)./areas.S_V(end,:))./areas.S_P_south; 


        % Define boundary lengths
        boundlength=pi*(2*nodes.nodesH_X(end,1:end-1)+distances.dX_nodesH_H(end,:)).*(distances.dX_nodesH_H(end,:).^2+distances.dY_nodesH_H(end,:).^2).^0.5;

        if strcmp(boundary.south,'Neumann')
            b = boundary.qdot.south*boundlength';
        elseif strcmp(boundary.south,'Robin')
            b = -boundary.alpha*boundary.Tinf*boundlength';
            
            if boundary.splitA
                Dx = - lambda(end,2:end-1).^2.*boundlength./areas.S_P_south;
                indx = ind(dimY,2:dimX-1)';
            else
                 D0 = D0 - boundary.alpha*lambda(end,2:end-1).^2.*boundlength./areas.S_P_south;
            end
        else
            error('unknown BC type')
        end
        
        ind0 = ind(dimY,2:dimX-1);
        ind_1 = ind(dimY-1,2:dimX-1);
        ind2 = ind(dimY-1,3:dimX);
        ind3 = ind(dimY,3:dimX);
        ind_3 = ind(dimY,1:dimX-2);
        ind_4 = ind(dimY-1,1:dimX-2);
        
        D =    [D0(:)   ; D_1(:)   ; D2(:)   ; D3(:)   ; D_3(:)   ; D_4(:)  ];
        indi = [ind0(:) ; ind0(:)  ; ind0(:) ; ind0(:) ; ind0(:)  ; ind0(:) ];
        indj = [ind0(:) ; ind_1(:) ; ind2(:) ; ind3(:) ; ind_3(:) ; ind_4(:)];
        
        a  =  D0(:) +  D_1(:) + D2(:) +  D3(:) + D_3(:) +  D_4(:);
        
end
end





% % East
%         D3 = ...
%             (-distances_south.dY_nodesC_V(1,2:dimX-1).*(-3*distances_south.dY_nodesV_V(1,3:dimX)+2*distances_south.dY_nodesV_H(1,2:dimX-1)-distances.dY_nodesV_H(dimY-1,2:dimX-1)) +...
%             -distances_south.dX_nodesC_V(1,2:dimX-1).*(-3*distances_south.dX_nodesV_V(1,3:dimX)+2*distances_south.dX_nodesV_H(1,2:dimX-1)-distances.dX_nodesV_H(dimY-1,2:dimX-1)))./(4*areas.S_S(1,2:dimX-1)) +...
%             (-distances.dY_nodesC_H(dimY-1,1:dimX-2).*(-distances.dY_nodesH_V(dimY-1,2:dimX-1)) +...
%             -distances.dX_nodesC_H(dimY-1,1:dimX-2).*(-distances.dX_nodesH_V(dimY-1,2:dimX-1)))./(4*areas.S_V(dimY-1,1:dimX-2));
%         
%         D3 = D3./areas.S_P_south;
%         
%         % West
%         D_3 = ...
%             (-distances.dY_nodesC_H(dimY-1,1:dimX-2).*(distances.dY_nodesH_V(dimY-1,1:dimX-2)) +...
%             -distances.dX_nodesC_H(dimY-1,1:dimX-2).*(distances.dX_nodesH_V(dimY-1,1:dimX-2)))./(4*areas.S_V(dimY-1,1:dimX-2)) +...
%             (distances_south.dY_nodesC_V(1,1:dimX-2).*(3*distances_south.dY_nodesV_V(1,1:dimX-2)+2*distances_south.dY_nodesV_H(1,1:dimX-2)-distances.dY_nodesV_H(dimY-1,1:dimX-2)) +...
%             distances_south.dX_nodesC_V(1,1:dimX-2).*(3*distances_south.dX_nodesV_V(1,1:dimX-2)+2*distances_south.dX_nodesV_H(1,1:dimX-2)-distances.dX_nodesV_H(dimY-1,1:dimX-2)))./(4*areas.S_S(1,1:dimX-2));
%         
%         D_3 = D_3./areas.S_P_south;
%         
%         % North
%         D_1 = ...
%             (-distances_south.dY_nodesC_V(1,2:dimX-1).*(-distances.dY_nodesV_H(dimY-1,2:dimX-1)+distances_south.dY_nodesV_V(1,2:dimX-1)) +...
%             -distances_south.dX_nodesC_V(1,2:dimX-1).*(-distances.dX_nodesV_H(dimY-1,2:dimX-1)+distances_south.dX_nodesV_V(1,2:dimX-1)))./(4*areas.S_S(1,2:dimX-1)) +...
%             (-distances.dY_nodesC_H(dimY-1,1:dimX-2).*(-distances.dY_nodesH_H(dimY-1,1:dimX-2)-distances.dY_nodesH_V(dimY-1,2:dimX-1)/4+distances.dY_nodesH_V(dimY-1,1:dimX-2)/4) +...
%             -distances.dX_nodesC_H(dimY-1,1:dimX-2).*(-distances.dX_nodesH_H(dimY-1,1:dimX-2)-distances.dX_nodesH_V(dimY-1,2:dimX-1)/4+distances.dX_nodesH_V(dimY-1,1:dimX-2)/4))./areas.S_V(dimY-1,1:dimX-2) +...
%             (distances_south.dY_nodesC_V(1,1:dimX-2).*(-distances.dY_nodesV_H(dimY-1,1:dimX-2)-distances_south.dY_nodesV_V(1,2:dimX-1)) +...
%             distances_south.dX_nodesC_V(1,1:dimX-2).*(-distances.dX_nodesV_H(dimY-1,1:dimX-2)-distances_south.dX_nodesV_V(1,2:dimX-1)))./(4*areas.S_S(1,1:dimX-2));
%         
%         D_1 = D_1./areas.S_P_south;
%         
%         % NW
%         D_4 = ...
%             (-distances.dY_nodesC_H(dimY-1,1:dimX-2).*distances.dY_nodesH_V(dimY-1,1:dimX-2) +...
%             -distances.dX_nodesC_H(dimY-1,1:dimX-2).*distances.dX_nodesH_V(dimY-1,1:dimX-2))./(4*areas.S_V(dimY-1,1:dimX-2)) +...
%             (distances_south.dY_nodesC_V(1,1:dimX-2).*(-distances.dY_nodesV_H(dimY-1,1:dimX-2)+distances_south.dY_nodesV_V(1,1:dimX-2)) +...
%             distances_south.dX_nodesC_V(1,1:dimX-2).*(-distances.dX_nodesV_H(dimY-1,1:dimX-2)+distances_south.dX_nodesV_V(1,1:dimX-2)))./(4*areas.S_S(1,1:dimX-2));
%         
%         D_4 = D_4./areas.S_P_south;
%         
%         % NE
%         D2 = ...
%             (-distances_south.dY_nodesC_V(1,2:dimX-1).*(-distances.dY_nodesV_H(dimY-1,2:dimX-1)-distances_south.dY_nodesV_V(1,3:dimX)) +...
%             -distances_south.dX_nodesC_V(1,2:dimX-1).*(-distances.dX_nodesV_H(dimY-1,2:dimX-1)-distances_south.dX_nodesV_V(1,3:dimX)))./(4*areas.S_S(1,2:dimX-1)) +...
%             (-distances.dY_nodesC_H(dimY-1,1:dimX-2).*(-distances.dY_nodesH_V(dimY-1,2:dimX-1)) +...
%             -distances.dX_nodesC_H(dimY-1,1:dimX-2).*(-distances.dX_nodesH_V(dimY-1,2:dimX-1)))./(4*areas.S_V(dimY-1,1:dimX-2));
%         
%         D2 = D2./areas.S_P_south;
%         
%         % P
%         D0 = ...
%             (-distances_south.dY_nodesC_V(1,2:dimX-1).*(3*distances_south.dY_nodesV_V(1,2:dimX-1)+2*distances_south.dY_nodesV_H(1,2:dimX-1)-distances.dY_nodesV_H(dimY-1,2:dimX-1)) +...
%             -distances_south.dX_nodesC_V(1,2:dimX-1).*(3*distances_south.dX_nodesV_V(1,2:dimX-1)+2*distances_south.dX_nodesV_H(1,2:dimX-1)-distances.dX_nodesV_H(dimY-1,2:dimX-1)))./(4*areas.S_S(1,2:dimX-1)) +...
%             (-distances.dY_nodesC_H(dimY-1,1:dimX-2).*(distances.dY_nodesH_H(dimY,1:dimX-2)-distances.dY_nodesH_V(dimY-1,2:dimX-1)/4+distances.dY_nodesH_V(dimY-1,1:dimX-2)/4) +...
%             -distances.dX_nodesC_H(dimY-1,1:dimX-2).*(distances.dX_nodesH_H(dimY,1:dimX-2)-distances.dX_nodesH_V(dimY-1,2:dimX-1)/4+distances.dX_nodesH_V(dimY-1,1:dimX-2)/4))./areas.S_V(dimY-1,1:dimX-2) +...
%             (distances_south.dY_nodesC_V(1,1:dimX-2).*(-3*distances_south.dY_nodesV_V(1,2:dimX-1)-distances.dY_nodesV_H(dimY-1,1:dimX-2)+2*distances_south.dY_nodesV_H(1,1:dimX-2)) +...
%             distances_south.dX_nodesC_V(1,1:dimX-2).*(-3*distances_south.dX_nodesV_V(1,2:dimX-1)-distances.dX_nodesV_H(dimY-1,1:dimX-2)+2*distances_south.dX_nodesV_H(1,1:dimX-2)))./(4*areas.S_S(1,1:dimX-2));
%         
%         D0 = D0./areas.S_P_south;


