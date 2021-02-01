%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the parameters of the simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Geometrical parameters

% Name for GIF
setting = 'pipe';

% Chose Configuration form 1 to 12

corr  = 1;  % on:1   off:0

pattern = [1,2,3,4,1,2,3,4,1,2,3,4; ...
           1,1,1,1,2,2,2,2,3,3,3,3];

if shape == 1
    upstream  = [0.096, 0.160, 0.224];
    flametube = [0.1, 0.15, 0.2, 0.4];      
    l.d1  = 0.07;                        %the top rectangle
    l.l1  = flametube(pattern(1,conf)) + 0.4*l.d1*corr; 
    l.d2 = 0.022;                        % the middle rectangle
    l.l2 = 0.078;
    l.d3  = 0.030;                       %the bottom rectangle
    l.d3u = 0.065; 
    l.l3o = 0.0675; 
    l.l3u = upstream(pattern(2,conf));  
    l.l3  = l.l3o + l.l3u;
else
    upstream  = [0.1248, 0.1888, 0.2528];
    flametube = [0.1, 0.15, 0.2, 0.4];
    l.d1  = 0.07;                        %the top rectangle
    l.l1  = flametube(pattern(1,conf)) + 0.4*l.d1*corr; 
    l.d2 = 0.02117;
    l.l2 = 0.1167;
    l.d3  = 0.065;                       %the bottom rectangle
    l.d3u = 0.065; 
    l.l3o = 0.0248; 
    l.l3u = upstream(pattern(2,conf))- 0.0248 ; 
    l.l3  = l.l3o + l.l3u;
end

% thining
%thin  = 0.01;
l.d2  = thin*l.d2;
l.d1  = thin*l.d1;
l.d3  = thin*l.d3;
l.d3u = thin*l.d3u;

% Geometrie of bottom cavitiy
% Vector for diameter

%l.d3vec = [l.d3,l.d3u*ones(1, round(l.l3u/0.005))];
%l.l3vec = [0, linspace(-l.l3o,-l.l3,length(l.d3vec)-1)];
%l.pp=spline(l.l3vec, [0, l.d3vec/2, 0]);
%yy = linspace(0, -l.l3);
%xx = ppval(l.pp, yy);
%l.pp = spline(yy, xx);

l.d3vec = [l.d3,l.d3u];

l.l3vec = [0, -l.l3o];

l.pp=spline(l.l3vec, [0, l.d3vec/2, slope]);
yy = linspace(0, -l.l3o);
xx = ppval(l.pp, yy);

xx2 = [xx(1:end-1),0.5*l.d3u*ones(1, round(l.l3u/0.005))];
yy2 = [yy(1:end-1),linspace(-l.l3o,-l.l3,length(xx2) - length(xx)+1) ];

yy3 = linspace(0, -l.l3);

xx3 = spline(yy2,xx2,yy3);


l.pp = spline(yy3, xx3);



% Approx size of dx and dy

size_dx = 0.001*thin;
size_dy = 0.003;

% Number of nodes for each direction in the middle rectangle

dim.dimX2 = round(0.5*l.d2/size_dx); % 0.5 because of radius/diameter
dim.dimY2 = round(l.l2/size_dy);

%sizes of the cells in the middle rectangle
x2 = l.d2/(dim.dimX2-1); 
y2 = l.l2/(dim.dimY2-1);

%sizes of the cells in the other rectangles 
%chosen closest to those of the cells in the middle rectangle
dim.dimX1bis = round((l.d1-l.d2)/x2);
dim.dimX1 = dim.dimX1bis + dim.dimX2; % number of nodes in the top rectangle
dim.dimY1 = round(l.l1/y2) + 1;
dim.dimX3bis = round((l.d3-l.d2)/x2);
dim.dimX3 = dim.dimX3bis + dim.dimX2; % number of nodes in the bottom rectangle
dim.dimY3 = round(l.l3/y2) + 1;

% Total number of nodes
dim.total = dim.dimX1*dim.dimY1 + dim.dimX2*(dim.dimY2-2) + dim.dimX3*dim.dimY3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lambda Values
c_sound1=sqrt(1.4*287*300);
c_sound2=sqrt(1.4*287*1600);
boundary(1).lambda=c_sound2*ones(dim.dimY1,dim.dimX1);
boundary(2).lambda=c_sound1*ones(dim.dimY2,dim.dimX2);
boundary(3).lambda=c_sound1*ones(dim.dimY3,dim.dimX3);

boundary(12).lambda = [boundary(1).lambda(end-1:end, 1:dim.dimX2); ...
                       boundary(2).lambda( 2 , :)];

boundary(23).lambda = [boundary(2).lambda(end-1, :); ...
                       boundary(3).lambda( 1:2 , 1:dim.dimX2)];
                   
boundary(11).lambda = [boundary(1).lambda(end-1:end, dim.dimX2-1:dim.dimX2+1); ...
                       boundary(2).lambda(2, dim.dimX2-1:end), boundary(2).lambda(2, end)];
                   
boundary(33).lambda = [boundary(2).lambda(end-1, dim.dimX2-1:end), boundary(2).lambda(end-1, end); ...
                       boundary(3).lambda(1:2, dim.dimX2-1:dim.dimX2+1)];
                   
% Parameter for Conjugated Heat Transfer
% Reflex = MagR*exp(1i*ArgR); 
% Z=(1+Reflex)/(1-Reflex);
% 
% boundary(1).alpha = 1i*omega/(Z*c_sound2);
boundary(1).alpha = 0;
boundary(1).splitA = 1; % Split Matrix A in constant and alpha dependent part
boundary(1).Tinf = 0;
boundary(2).alpha = 0;
boundary(2).Tinf = 0;
boundary(3).alpha = 0;
boundary(3).splitA = 1;
boundary(3).Tinf = 0;
                  
% Boundary conditions %Pay attention to correct spelling, there is no
% complete warning for wrong spelling so far. Condition will be Robin then.
% Upper body
boundary(1).south = 'Neumann';
boundary(1).north = 'Robin';
%boundary(1).north = 'Dirichlet';
boundary(1).east = 'Neumann';
boundary(1).west = 'Neumann';
boundary(1).isInner = 0;
% Middle body
boundary(2).south = 'Neumann';
boundary(2).north = 'Neumann';
boundary(2).east = 'Neumann';
boundary(2).west = 'Neumann';
boundary(2).isInner = 0;
% Lower body
boundary(3).south = 'Robin';
%boundary(3).south = 'Neumann';
boundary(3).north = 'Neumann';
boundary(3).east = 'Neumann';
boundary(3).west = 'Neumann';
boundary(3).isInner = 0;

% Interfaces
boundary(12).isInner = 1;
boundary(12).west = 'Neumann';
boundary(23).isInner = 1;
boundary(23).west = 'Neumann';

% Corners (Dirichlet and Robin to be added)
boundary(11).isInner = 2;
boundary(33).isInner = 3;

% Values for  BC
% Upper body
boundary(1).TD.north = 0;
boundary(1).TD.south= 0;
boundary(1).TD.west=0;
boundary(1).TD.east=0;
boundary(1).qdot.north = 0;
boundary(1).qdot.south = 0;
boundary(1).qdot.west = 0;
boundary(1).qdot.east = 0;
% Middle body
boundary(2).TD.north = 0;
boundary(2).TD.south= 10;
boundary(2).TD.west= 0;
boundary(2).TD.east=0;
boundary(2).qdot.north = 0;
boundary(2).qdot.south = 0;
boundary(2).qdot.west = 0;
boundary(2).qdot.east = 0;
% Lower body
boundary(3).TD.north = 0;
boundary(3).TD.south= 0;
boundary(3).TD.west= 0;
boundary(3).TD.east=0;
boundary(3).qdot.north = 0;
boundary(3).qdot.south = 0;
boundary(3).qdot.west = 0;
boundary(3).qdot.east = 0;
% Interfaces
boundary(12).TD.west= 0;
boundary(12).qdot.west = 0;
boundary(23).TD.west= 0;
boundary(23).qdot.west = 0;

bound = boundary;



