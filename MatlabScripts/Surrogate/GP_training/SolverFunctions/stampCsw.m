function [indi,indj,D,b,Dx,indx] = stampCsw(nodes, dimY, dimX, distances, distances_south, distances_west, areas, ind, boundary, dim)

boundleft='west';
boundright='south';

Dx = [];
indx = [];

% Decide for Dirichlet or flow boundary
lefttype=strcat('boundary.',boundleft);
righttype=strcat('boundary.',boundright);

if strcmp(eval(lefttype),'Dirichlet') || strcmp(eval(righttype),'Dirichlet')
    
    if strcmp(eval(lefttype),'Dirichlet')
        TempDiri=eval(strcat('boundary.TD.',boundleft));
    else
        TempDiri=eval(strcat('boundary.TD.',boundright));
    end
    
    indi = ind(dimY,1);
    indj = ind(dimY,1);
    
    D=1;
    b = TempDiri;
    disp('Edge Dirichlet')
    
else
    
    % X-Position
    rx_m = nodes.nodesC_X(end,1);
    rx_r = nodes.nodesH_X(end,1);
    
    % X Define distances
    dx_r_dL = -distances.dX_nodesH_V(end,1);
    dx_dL_L = -distances_west.dX_nodesH_H(end-1,1);
    dx_L_P  =  distances_west.dX_nodesH_V(end,1);
    dx_P_r  =  distances_west.dX_nodesH_H(end,1);
    
    dx_R_dR = -distances_south.dX_nodesV_V(1,2);
    dx_dR_l = -distances.dX_nodesV_H(end,1);
    dx_l_P  =  distances_south.dX_nodesV_V(1,1);
    dx_P_R  =  distances_south.dX_nodesV_H(1,1);
    
    dx_m_l  = -distances_west.dX_nodesC_H(end,1);
    dx_r_m  = -distances_south.dX_nodesC_V(1,1);
    
    % Y Define distances
    dy_r_dL = -distances.dY_nodesH_V(end,1);
    dy_dL_L = -distances_west.dY_nodesH_H(end-1,1);
    dy_L_P  =  distances_west.dY_nodesH_V(end,1);
    dy_P_r  =  distances_west.dY_nodesH_H(end,1);
    
    dy_R_dR = -distances_south.dY_nodesV_V(1,2);
    dy_dR_l = -distances.dY_nodesV_H(end,1);
    dy_l_P  =  distances_south.dY_nodesV_V(1,1);
    dy_P_R  =  distances_south.dY_nodesV_H(1,1);
    
    dy_m_l  = -distances_west.dY_nodesC_H(end,1);
    dy_r_m  = -distances_south.dY_nodesC_V(1,1);
    
    % Define boundary lengths
    dBl=pi*(2*nodes.nodesV_X(end,1)+dx_l_P)*(dx_l_P^2+dy_l_P^2)^0.5;
    dBr=pi*(2*nodes.nodesH_X(end,1)-dx_P_r)*(dx_P_r^2+dy_P_r^2)^0.5;
    
    % Define Areas
    S_gamma = areas.S_S(1,1);
    S_DELTA = areas.S_W(end,1);
    S_P = areas.S_P_sw;
    
    
    % Lambda
    lambda_R = boundary.lambda(end,2);
    lambda_D = boundary.lambda(end-1,2);
    lambda_L = boundary.lambda(end-1,1);
    lambda_P = boundary.lambda(end,1);
    
    % Stecil
    
	% Right 
	D_R=((pi*(dx_r_m*(dx_P_R/2 + (3*dx_R_dR)/4 + dx_dR_l/4) + dy_r_m*(dy_P_R/2 + (3*dy_R_dR)/4 + dy_dR_l/4))*(dx_m_l + 2*rx_r)*(lambda_D/8 + lambda_L/8 + (3*lambda_P)/8 + (3*lambda_R)/8)^2)/S_gamma + (pi*(dx_m_l + 2*rx_m)*(dx_m_l*(dx_P_r/4 + dx_r_dL/4) + dy_m_l*(dy_P_r/4 + dy_r_dL/4))*(lambda_D/8 + (3*lambda_L)/8 + (3*lambda_P)/8 + lambda_R/8)^2)/S_DELTA)/S_P; 

	% Domain 
	D_D=((pi*(dx_m_l + 2*rx_m)*(dx_m_l*(dx_dL_L/4 + dx_r_dL/4) + dy_m_l*(dy_dL_L/4 + dy_r_dL/4))*(lambda_D/8 + (3*lambda_L)/8 + (3*lambda_P)/8 + lambda_R/8)^2)/S_DELTA + (pi*(dx_m_l + 2*rx_r)*(dx_r_m*(dx_R_dR/4 + dx_dR_l/4) + dy_r_m*(dy_R_dR/4 + dy_dR_l/4))*(lambda_D/8 + lambda_L/8 + (3*lambda_P)/8 + (3*lambda_R)/8)^2)/S_gamma)/S_P; 

	% Left 
	D_L=((pi*(dx_m_l*(dx_L_P/2 + (3*dx_dL_L)/4 + dx_r_dL/4) + dy_m_l*(dy_L_P/2 + (3*dy_dL_L)/4 + dy_r_dL/4))*(dx_m_l + 2*rx_m)*(lambda_D/8 + (3*lambda_L)/8 + (3*lambda_P)/8 + lambda_R/8)^2)/S_DELTA + (pi*(dx_m_l + 2*rx_r)*(dx_r_m*(dx_l_P/4 + dx_dR_l/4) + dy_r_m*(dy_l_P/4 + dy_dR_l/4))*(lambda_D/8 + lambda_L/8 + (3*lambda_P)/8 + (3*lambda_R)/8)^2)/S_gamma)/S_P; 

	% Point 
	D_P=((pi*(dx_m_l*(dx_L_P/2 + (3*dx_P_r)/4 + dx_r_dL/4) + dy_m_l*(dy_L_P/2 + (3*dy_P_r)/4 + dy_r_dL/4))*(dx_m_l + 2*rx_m)*(lambda_D/8 + (3*lambda_L)/8 + (3*lambda_P)/8 + lambda_R/8)^2)/S_DELTA + (pi*(dx_r_m*(dx_P_R/2 + (3*dx_l_P)/4 + dx_dR_l/4) + dy_r_m*(dy_P_R/2 + (3*dy_l_P)/4 + dy_dR_l/4))*(dx_m_l + 2*rx_r)*(lambda_D/8 + lambda_L/8 + (3*lambda_P)/8 + (3*lambda_R)/8)^2)/S_gamma)/S_P; 


    
    % RHS
    qdotL = eval(strcat('boundary.qdot.',boundleft));
    qdotR = eval(strcat('boundary.qdot.',boundright));
    
    if strcmp(eval(lefttype),'Neumann') && strcmp(eval(righttype),'Neumann')
        b = +qdotL*dBl + qdotR*dBr;
%         disp('Edge Neumann - Neumann')
        
    elseif strcmp(eval(lefttype),'Neumann') && strcmp(eval(righttype),'Robin')
        b = +qdotL*dBl - boundary.alpha*boundary.Tinf*dBr;
        
        if boundary.splitA
            Dx = lambda_P^2*dBr/S_P;
            indx = ind(dimY,1);
        else
            D_P = D_P - boundary.alpha*lambda_P^2*dBr/S_P;
        end
        
%         disp('Edge Neumann - Robin')
        
    elseif strcmp(eval(lefttype),'Robin') && strcmp(eval(righttype),'Neumann')
        b = -boundary.alpha*boundary.Tinf*dBl + qdotR*dBr;
        
        if boundary.splitA
            Dx = - lambda_P^2*dBl/S_P;
            indx = ind(dimY,1);
        else
            D_P = D_P - boundary.alpha*lambda_P^2*dBl/S_P;
        end
        
%         disp('Edge Robin - Neumann')
        
    elseif strcmp(eval(lefttype),'Robin') && strcmp(eval(righttype),'Robin')
        b = -boundary.alpha*boundary.Tinf*(dBl+dBr);
        
        if boundary.splitA
            Dx = - lambda_P^2*(dBl+dBr)/S_P;
            indx = ind(dimY,1);
        else
            D_P = D_P - boundary.alpha*lambda_P^2*(dBl+dBr)/S_P;
        end
        
%         disp('Edge Robin - Robin')
    
    else
        error('unknown BC type')
    end
    
    ind0 = ind(dimY,1);
    ind_1 = ind(dimY-1,1);
    ind2 = ind(dimY-1,2);
    ind3 = ind(dimY,2);
    
    D =    [ D_P   ; D_R   ; D_L    ; D_D  ];
    indi = [ ind0  ; ind0  ; ind0   ; ind0 ];
    indj = [ ind0  ; ind3  ; ind_1  ; ind2 ];
    
end
end





