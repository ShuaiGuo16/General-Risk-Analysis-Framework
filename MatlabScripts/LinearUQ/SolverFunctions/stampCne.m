function [indi,indj,D,b,Dx,indx] = stampCne(nodes, dimY, dimX, distances, distances_north, distances_east, areas, ind, boundary, dim)

boundleft='east';
boundright='north';

Dx = [];
indx = [];

% Decide for Dirichlet or flow boundary
lefttype=strcat('boundary.',boundleft);
righttype=strcat('boundary.',boundright);

if  strcmp(eval(lefttype),'Dirichlet') || strcmp(eval(righttype),'Dirichlet')
    
    if strcmp(eval(lefttype),'Dirichlet')
        TempDiri=eval(strcat('boundary.TD.',boundleft));
    else
        TempDiri=eval(strcat('boundary.TD.',boundright));
    end
    
    indi = ind(1,dimX);
    indj = ind(1,dimX);
    
    D=1;
    b = TempDiri;
    disp('Edge Dirichlet')
    
else
    
    % X-Position
    rx_m = nodes.nodesC_X(1,end);
    rx_r = nodes.nodesH_X(1,end);
    
    % X Define distances
    dx_r_dL =  distances.dX_nodesH_V(1,end);
    dx_dL_L =  distances_east.dX_nodesH_H(2,1);
    dx_L_P  = -distances_east.dX_nodesH_V(1,1);
    dx_P_r  = -distances_east.dX_nodesH_H(1,1);
    
    dx_R_dR =  distances_north.dX_nodesV_V(1,end-1);
    dx_dR_l =  distances.dX_nodesV_H(1,end);
    dx_l_P  = -distances_north.dX_nodesV_V(1,end);
    dx_P_R  = -distances_north.dX_nodesV_H(1,end);
    
    dx_m_l  =  distances_east.dX_nodesC_H(1,1);
    dx_r_m  =  distances_north.dX_nodesC_V(1,end);
    
    % Y Define distances
    dy_r_dL =  distances.dY_nodesH_V(1,end);
    dy_dL_L =  distances_east.dY_nodesH_H(2,1);
    dy_L_P  = -distances_east.dY_nodesH_V(1,1);
    dy_P_r  = -distances_east.dY_nodesH_H(1,1);
    
    dy_R_dR =  distances_north.dY_nodesV_V(1,end-1);
    dy_dR_l =  distances.dY_nodesV_H(1,end);
    dy_l_P  = -distances_north.dY_nodesV_V(1,end);
    dy_P_R  = -distances_north.dY_nodesV_H(1,end);
    
    dy_m_l  =  distances_east.dY_nodesC_H(1,1);
    dy_r_m  =  distances_north.dY_nodesC_V(1,end);
    
    % Define boundary lengths
    dBl=pi*(2*nodes.nodesV_X(1,end)+dx_l_P)*(dx_l_P^2+dy_l_P^2)^0.5;
    dBr=pi*(2*nodes.nodesH_X(1,end)-dx_P_r)*(dx_P_r^2+dy_P_r^2)^0.5;
    
    % Define Areas
    S_gamma = areas.S_N(1,end);
    S_DELTA = areas.S_E(1,1);
    S_P = areas.S_P_ne;
    
    % Lambda
    lambda_R = boundary.lambda(1,end);
    lambda_D = boundary.lambda(2,end-1);
    lambda_L = boundary.lambda(2,end);
    lambda_P = boundary.lambda(1,end-1);
    
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
            Dx = - lambda_P^2*dBr/S_P;
            indx = ind(1,dimX);
        else
        D_P = D_P - boundary.alpha*lambda_P^2*dBr/S_P;
        end
        
%         disp('Edge Neumann - Robin')
        
    elseif strcmp(eval(lefttype),'Robin') && strcmp(eval(righttype),'Neumann')
        b = -boundary.alpha*boundary.Tinf*dBl + qdotR*dBr;
        
        if boundary.splitA
            Dx = - lambda_P^2*dBl/S_P;
            indx = ind(1,dimX);
        else
        D_P = D_P - boundary.alpha*lambda_P^2*dBl/S_P;
        end
        
%         disp('Edge Robin - Neumann')
        
    elseif strcmp(eval(lefttype),'Robin') && strcmp(eval(righttype),'Robin')
        b = -boundary.alpha*boundary.Tinf*(dBl+dBr);
        
        if boundary.splitA
            Dx = - lambda_P^2*(dBl+dBr)/S_P;
            indx = ind(1,dimX);
        else
        D_P = D_P - boundary.alpha*lambda_P^2*(dBl+dBr)/S_P;
        end
        
%         disp('Edge Robin - Robin')
    else
        error('unknown BC type')
    end
    
    ind0 = ind(1,dimX);
    ind1 = ind(2,dimX);
    ind_2 = ind(2,dimX-1);
    ind_3 = ind(1,dimX-1);
    
    D =    [ D_P   ; D_R    ; D_L   ; D_D   ];
    indi = [ ind0  ; ind0   ; ind0  ; ind0  ];
    indj = [ ind0  ; ind_3  ; ind1  ; ind_2 ];
    
end
end