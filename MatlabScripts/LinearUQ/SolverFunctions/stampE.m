function [indi,indj,D,b] = stampE(nodes, dimY, dimX, distances, distances_east, areas, ind, boundary,dim)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

switch boundary.east
    case 'Dirichlet'
        
        D = ones(dimY-2,1);
        indi = reshape(ind(2:dimY-1,dimX),(dimY-2),1);
        indj = reshape(ind(2:dimY-1,dimX),(dimY-2),1);
        b = boundary.TD.east*ones(dimY-2,1);
        
    case {'Neumann'}
        
        %Stecil
        lambda=boundary.lambda;
        
        % Right 
        D_1=((pi.*(-distances_east.dX_nodesC_H(1:end-1).*(-distances_east.dX_nodesH_V(1:end-1)./2 + (3.*-distances_east.dX_nodesH_H(1:end-2))./4 + distances.dX_nodesH_V(1:end-1,end)./4) + -distances_east.dY_nodesC_H(1:end-1).*(-distances_east.dY_nodesH_V(1:end-1)./2 + (3.*-distances_east.dY_nodesH_H(1:end-2))./4 + distances.dY_nodesH_V(1:end-1,end)./4)).*(-distances_east.dX_nodesC_H(1:end-1) + 2.*nodes.nodesV_X(1:end-1,end)).*((3.*lambda(2:end-1,end))./8 + lambda(2:end-1,end-1)./8 + (3.*lambda(1:end-2,end))./8 + lambda(1:end-2,end-1)./8).^2)./areas.S_E(1:end-1) + (pi.*(distances.dX_nodesC_V(:,end) + 2.*nodes.nodesC_X(1:end-1,end)).*((-distances.dX_nodesV_H(1:end-1,end).*distances.dX_nodesC_V(:,end))./4 + (-distances.dY_nodesV_H(1:end-1,end).*distances.dY_nodesC_V(:,end))./4).*(lambda(2:end-1,end)./2 + lambda(2:end-1,end-1)./2).^2)./areas.S_H(:,end))./areas.S_P_east; 

        % DomainRight 
        D_4=((pi.*(-distances_east.dX_nodesC_H(1:end-1) + 2.*nodes.nodesV_X(1:end-1,end)).*(-distances_east.dX_nodesC_H(1:end-1).*(-distances_east.dX_nodesH_H(1:end-2)./4 + distances.dX_nodesH_V(1:end-1,end)./4) + -distances_east.dY_nodesC_H(1:end-1).*(-distances_east.dY_nodesH_H(1:end-2)./4 + distances.dY_nodesH_V(1:end-1,end)./4)).*((3.*lambda(2:end-1,end))./8 + lambda(2:end-1,end-1)./8 + (3.*lambda(1:end-2,end))./8 + lambda(1:end-2,end-1)./8).^2)./areas.S_E(1:end-1) + (pi.*(distances.dX_nodesC_V(:,end) + 2.*nodes.nodesC_X(1:end-1,end)).*((-distances.dX_nodesV_H(1:end-1,end).*distances.dX_nodesC_V(:,end))./4 + (-distances.dY_nodesV_H(1:end-1,end).*distances.dY_nodesC_V(:,end))./4).*(lambda(2:end-1,end)./2 + lambda(2:end-1,end-1)./2).^2)./areas.S_H(:,end))./areas.S_P_east; 

        % Domain 
        D_3=((pi.*(distances_east.dX_nodesC_H(2:end) + 2.*nodes.nodesC_X(2:end,end)).*(distances_east.dX_nodesC_H(2:end).*(-distances_east.dX_nodesH_H(2:end-1)./4 + distances.dX_nodesH_V(2:end,end)./4) + distances_east.dY_nodesC_H(2:end).*(-distances_east.dY_nodesH_H(2:end-1)./4 + distances.dY_nodesH_V(2:end,end)./4)).*((3.*lambda(2:end-1,end))./8 + lambda(2:end-1,end-1)./8 + (3.*lambda(3:end,end))./8 + lambda(3:end,end-1)./8).^2)./areas.S_E(2:end) + (pi.*(-distances_east.dX_nodesC_H(1:end-1) + 2.*nodes.nodesV_X(1:end-1,end)).*(-distances_east.dX_nodesC_H(1:end-1).*(distances_east.dX_nodesH_H(2:end-1)./4 + distances.dX_nodesH_V(1:end-1,end)./4) + -distances_east.dY_nodesC_H(1:end-1).*(distances_east.dY_nodesH_H(2:end-1)./4 + distances.dY_nodesH_V(1:end-1,end)./4)).*((3.*lambda(2:end-1,end))./8 + lambda(2:end-1,end-1)./8 + (3.*lambda(1:end-2,end))./8 + lambda(1:end-2,end-1)./8).^2)./areas.S_E(1:end-1) + (pi.*(distances.dX_nodesC_V(:,end).*(distances.dX_nodesV_H(2:end,end)./4 + -distances.dX_nodesV_H(1:end-1,end)./4 + distances.dX_nodesV_V(:,end-1)) + distances.dY_nodesC_V(:,end).*(distances.dY_nodesV_H(2:end,end)./4 + -distances.dY_nodesV_H(1:end-1,end)./4 + distances.dY_nodesV_V(:,end-1))).*(distances.dX_nodesC_V(:,end) + 2.*nodes.nodesC_X(1:end-1,end)).*(lambda(2:end-1,end)./2 + lambda(2:end-1,end-1)./2).^2)./areas.S_H(:,end))./areas.S_P_east; 

        % DomainLeft 
        D_2=((pi.*(distances_east.dX_nodesC_H(2:end) + 2.*nodes.nodesC_X(2:end,end)).*(distances_east.dX_nodesC_H(2:end).*(distances_east.dX_nodesH_H(3:end)./4 + distances.dX_nodesH_V(2:end,end)./4) + distances_east.dY_nodesC_H(2:end).*(distances_east.dY_nodesH_H(3:end)./4 + distances.dY_nodesH_V(2:end,end)./4)).*((3.*lambda(2:end-1,end))./8 + lambda(2:end-1,end-1)./8 + (3.*lambda(3:end,end))./8 + lambda(3:end,end-1)./8).^2)./areas.S_E(2:end) + (pi.*(distances.dX_nodesC_V(:,end) + 2.*nodes.nodesC_X(1:end-1,end)).*((distances.dX_nodesV_H(2:end,end).*distances.dX_nodesC_V(:,end))./4 + (distances.dY_nodesV_H(2:end,end).*distances.dY_nodesC_V(:,end))./4).*(lambda(2:end-1,end)./2 + lambda(2:end-1,end-1)./2).^2)./areas.S_H(:,end))./areas.S_P_east; 

        % Left 
        D1=((pi.*(distances_east.dX_nodesC_H(2:end).*(-distances_east.dX_nodesH_V(2:end)./2 + (3.*distances_east.dX_nodesH_H(3:end))./4 + distances.dX_nodesH_V(2:end,end)./4) + distances_east.dY_nodesC_H(2:end).*(-distances_east.dY_nodesH_V(2:end)./2 + (3.*distances_east.dY_nodesH_H(3:end))./4 + distances.dY_nodesH_V(2:end,end)./4)).*(distances_east.dX_nodesC_H(2:end) + 2.*nodes.nodesC_X(2:end,end)).*((3.*lambda(2:end-1,end))./8 + lambda(2:end-1,end-1)./8 + (3.*lambda(3:end,end))./8 + lambda(3:end,end-1)./8).^2)./areas.S_E(2:end) + (pi.*(distances.dX_nodesC_V(:,end) + 2.*nodes.nodesC_X(1:end-1,end)).*((distances.dX_nodesV_H(2:end,end).*distances.dX_nodesC_V(:,end))./4 + (distances.dY_nodesV_H(2:end,end).*distances.dY_nodesC_V(:,end))./4).*(lambda(2:end-1,end)./2 + lambda(2:end-1,end-1)./2).^2)./areas.S_H(:,end))./areas.S_P_east; 

        % Boundary 
        D0=((pi.*(distances_east.dX_nodesC_H(2:end).*(-distances_east.dX_nodesH_V(2:end)./2 + (3.*-distances_east.dX_nodesH_H(2:end-1))./4 + distances.dX_nodesH_V(2:end,end)./4) + distances_east.dY_nodesC_H(2:end).*(-distances_east.dY_nodesH_V(2:end)./2 + (3.*-distances_east.dY_nodesH_H(2:end-1))./4 + distances.dY_nodesH_V(2:end,end)./4)).*(distances_east.dX_nodesC_H(2:end) + 2.*nodes.nodesC_X(2:end,end)).*((3.*lambda(2:end-1,end))./8 + lambda(2:end-1,end-1)./8 + (3.*lambda(3:end,end))./8 + lambda(3:end,end-1)./8).^2)./areas.S_E(2:end) + (pi.*(-distances_east.dX_nodesC_H(1:end-1).*(-distances_east.dX_nodesH_V(1:end-1)./2 + (3.*distances_east.dX_nodesH_H(2:end-1))./4 + distances.dX_nodesH_V(1:end-1,end)./4) + -distances_east.dY_nodesC_H(1:end-1).*(-distances_east.dY_nodesH_V(1:end-1)./2 + (3.*distances_east.dY_nodesH_H(2:end-1))./4 + distances.dY_nodesH_V(1:end-1,end)./4)).*(-distances_east.dX_nodesC_H(1:end-1) + 2.*nodes.nodesV_X(1:end-1,end)).*((3.*lambda(2:end-1,end))./8 + lambda(2:end-1,end-1)./8 + (3.*lambda(1:end-2,end))./8 + lambda(1:end-2,end-1)./8).^2)./areas.S_E(1:end-1) + (pi.*(distances.dX_nodesC_V(:,end).*(-distances.dX_nodesV_V(:,end) + distances.dX_nodesV_H(2:end,end)./4 + -distances.dX_nodesV_H(1:end-1,end)./4) + distances.dY_nodesC_V(:,end).*(-distances.dY_nodesV_V(:,end) + distances.dY_nodesV_H(2:end,end)./4 + -distances.dY_nodesV_H(1:end-1,end)./4)).*(distances.dX_nodesC_V(:,end) + 2.*nodes.nodesC_X(1:end-1,end)).*(lambda(2:end-1,end)./2 + lambda(2:end-1,end-1)./2).^2)./areas.S_H(:,end))./areas.S_P_east; 


        
        % Define boundary lengths
        boundlength=(distances.dX_nodesV_V(:,end).^2+distances.dY_nodesV_V(:,end).^2).^0.5;
        b = boundary.qdot.east*boundlength;
        
        ind0 = ind(2:dimY-1,dimX);
        ind1 = ind(3:dimY,dimX);
        ind_1 = ind(1:dimY-2,dimX);
        ind_2 = ind(3:dimY,dimX-1);
        ind_3 = ind(2:dimY-1,dimX-1);
        ind_4 = ind(1:dimY-2,dimX-1);
        
        D =    [D0(:)   ; D1(:)   ; D_1(:)   ; D_2(:)   ; D_3(:)   ; D_4(:)  ];
        indi = [ind0(:) ; ind0(:) ; ind0(:)  ; ind0(:)  ; ind0(:)  ; ind0(:) ];
        indj = [ind0(:) ; ind1(:) ; ind_1(:) ; ind_2(:) ; ind_3(:) ; ind_4(:)];
        
        %a  =  D0(:) + D1(:) + D_1(:) +D_2(:) + D_3(:) + D_4(:)
        
end
end




%
%
% % South
%         D1= ...
%             (distances_east.dX_nodesC_H(2:end).*(-distances_east.dX_nodesH_V(2:end)/2 + (3*distances_east.dX_nodesH_H(3:end))/4 + distances.dX_nodesH_V(2:end,end)/4))./areas.S_E(2:end) + ...
%             (distances_east.dY_nodesC_H(2:end).*(-distances_east.dY_nodesH_V(2:end)/2 + (3*distances_east.dY_nodesH_H(3:end))/4 + distances.dY_nodesH_V(2:end,end)/4))./areas.S_E(2:end) + ...
%             (distances.dX_nodesV_H(2:end,end).*distances.dX_nodesC_V(:,end))./(4*areas.S_H(:,end)) + ...
%             (distances.dY_nodesV_H(2:end,end).*distances.dY_nodesC_V(:,end))./(4*areas.S_H(:,end));
%
%         D1 = D1./areas.S_P_east;
%
%         % Southwest
%         D_2= ...
%             (distances_east.dX_nodesC_H(2:end).*(distances_east.dX_nodesH_H(3:end)/4 + distances.dX_nodesH_V(2:end,end)/4))./areas.S_E(2:end) + ...
%             (distances_east.dY_nodesC_H(2:end).*(distances_east.dY_nodesH_H(3:end)/4 + distances.dY_nodesH_V(2:end,end)/4))./areas.S_E(2:end) + ...
%             (distances.dX_nodesV_H(2:end,end).*distances.dX_nodesC_V(:,end))./(4*areas.S_H(:,end)) + ...
%             (distances.dY_nodesV_H(2:end,end).*distances.dY_nodesC_V(:,end))./(4*areas.S_H(:,end));
%
%         D_2 = D_2./areas.S_P_east;
%
%         % West
%         D_3= ...
%             (distances_east.dX_nodesC_H(2:end).*(-distances_east.dX_nodesH_H(2:end-1)/4 + distances.dX_nodesH_V(2:end,end)/4))./areas.S_E(2:end) + ...
%             (-distances_east.dX_nodesC_H(1:end-1).*(distances_east.dX_nodesH_H(2:end-1)/4 + distances.dX_nodesH_V(1:end-1,end)/4))./areas.S_E(1:end-1) + ...
%             (distances_east.dY_nodesC_H(2:end).*(-distances_east.dY_nodesH_H(2:end-1)/4 + distances.dY_nodesH_V(2:end,end)/4))./areas.S_E(2:end) + ...
%             (-distances_east.dY_nodesC_H(1:end-1).*(distances_east.dY_nodesH_H(2:end-1)/4 + distances.dY_nodesH_V(1:end-1,end)/4))./areas.S_E(1:end-1) + ...
%             (distances.dX_nodesC_V(:,end).*(distances.dX_nodesV_H(2:end,end)/4 - distances.dX_nodesV_H(1:end-1,end)/4 + distances.dX_nodesV_V(:,end-1)))./areas.S_H(:,end) + ...
%             (distances.dY_nodesC_V(:,end).*(distances.dY_nodesV_H(2:end,end)/4 - distances.dY_nodesV_H(1:end-1,end)/4 + distances.dY_nodesV_V(:,end-1)))./areas.S_H(:,end);
%
%         D_3 = D_3./areas.S_P_east;
%
%         % Northwest
%         D_4= ...
%             (-distances_east.dX_nodesC_H(1:end-1).*(-distances_east.dX_nodesH_H(1:end-2)/4 + distances.dX_nodesH_V(1:end-1,end)/4))./areas.S_E(1:end-1) + ...
%             (-distances_east.dY_nodesC_H(1:end-1).*(-distances_east.dY_nodesH_H(1:end-2)/4 + distances.dY_nodesH_V(1:end-1,end)/4))./areas.S_E(1:end-1) + ...
%             (-distances.dX_nodesV_H(1:end-1,end).*distances.dX_nodesC_V(:,end))./(4*areas.S_H(:,end)) + ...
%             (-distances.dY_nodesV_H(1:end-1,end).*distances.dY_nodesC_V(:,end))./(4*areas.S_H(:,end));
%
%         D_4 = D_4./areas.S_P_east;
%
%         % North
%         D_1= ...
%             (distances_east.dX_nodesC_H(1:end-1).*(distances_east.dX_nodesH_V(1:end-1)/2 + (3*distances_east.dX_nodesH_H(1:end-2))/4 - distances.dX_nodesH_V(1:end-1,end)/4))./areas.S_E(1:end-1) + ...
%             (distances_east.dY_nodesC_H(1:end-1).*(distances_east.dY_nodesH_V(1:end-1)/2 + (3*distances_east.dY_nodesH_H(1:end-2))/4 - distances.dY_nodesH_V(1:end-1,end)/4))./areas.S_E(1:end-1) + ...
%             (-distances.dX_nodesV_H(1:end-1,end).*distances.dX_nodesC_V(:,end))./(4*areas.S_H(:,end)) + ...
%             (-distances.dY_nodesV_H(1:end-1,end).*distances.dY_nodesC_V(:,end))./(4*areas.S_H(:,end));
%
%         D_1 = D_1./areas.S_P_east;
%
%         % P
%         D0= ...
%             (distances_east.dX_nodesC_H(2:end).*(-distances_east.dX_nodesH_V(2:end)/2 + (3*-distances_east.dX_nodesH_H(2:end-1))/4 + distances.dX_nodesH_V(2:end,end)/4))./areas.S_E(2:end) + ...
%             (-distances_east.dX_nodesC_H(1:end-1).*(-distances_east.dX_nodesH_V(1:end-1)/2 + (3*distances_east.dX_nodesH_H(2:end-1))/4 + distances.dX_nodesH_V(1:end-1,end)/4))./areas.S_E(1:end-1) + ...
%             (distances_east.dY_nodesC_H(2:end).*(-distances_east.dY_nodesH_V(2:end)/2 - (3*distances_east.dY_nodesH_H(2:end-1))/4 + distances.dY_nodesH_V(2:end,end)/4))./areas.S_E(2:end) + ...
%             (-distances_east.dY_nodesC_H(1:end-1).*(-distances_east.dY_nodesH_V(1:end-1)/2 + (3*distances_east.dY_nodesH_H(2:end-1))/4 + distances.dY_nodesH_V(1:end-1,end)/4))./areas.S_E(1:end-1) + ...
%             (distances.dX_nodesC_V(:,end).*(-distances.dX_nodesV_V(:,end) + distances.dX_nodesV_H(2:end,end)/4 - distances.dX_nodesV_H(1:end-1,end)/4))./areas.S_H(:,end) + ...
%             (distances.dY_nodesC_V(:,end).*(-distances.dY_nodesV_V(:,end) + distances.dY_nodesV_H(2:end,end)/4 - distances.dY_nodesV_H(1:end-1,end)/4))./areas.S_H(:,end);
%
%         D0 = D0./areas.S_P_east;

