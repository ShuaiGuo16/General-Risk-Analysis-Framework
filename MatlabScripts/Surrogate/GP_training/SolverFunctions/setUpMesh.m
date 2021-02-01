function [ mesh, index ] = setUpMesh(l,dim)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%File setUpMesh
%Calculates all the parameters for setting up the mesh 
%and sets it up.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dimX1 = dim.dimX1;
dimY1 = dim.dimY1;
dimX2 = dim.dimX2;
dimY2 = dim.dimY2;
dimX3 = dim.dimX3;
dimY3 = dim.dimY3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create the mesh for the top rectangle
%Matrices containing the position of the nodes for each rectangle
mesh.X1 = zeros(dimY1,dimX1);
mesh.X1(:,(1:dimX2)) = repmat(linspace(0,l.d2/2,dimX2),dimY1,1);
mesh.X1(:,(dimX2:dimX1)) = repmat(linspace(l.d2/2,l.d1/2,dimX1-dimX2+1),dimY1,1);

mesh.Y1 = repmat(linspace(l.l1+l.l2,l.l2,dimY1)',1,dimX1);

mesh.X2 = repmat(linspace(0,l.d2/2,dimX2),dimY2,1);
mesh.Y2 = repmat(linspace(l.l2,0,dimY2)',1,dimX2);


mesh.X3 = zeros(dimY3,dimX3);
mesh.X3(1,(1:dimX2)) = linspace(0,l.d2/2,dimX2);
mesh.X3(1,(dimX2:dimX3)) = linspace(l.d2/2,l.d3/2,dimX3-dimX2+1);

% Make Spline
mesh.X3 = ppval(l.pp,linspace(0,-l.l3,dimY3)')*mesh.X3(1,:)/(l.d3/2);
mesh.Y3 = repmat(linspace(0,-l.l3,dimY3)',1,dimX3);

% Interface meshes
mesh.X12=[mesh.X1(end-1:end, 1 : dimX2); ...
          mesh.X2(  2, : )];
          
mesh.Y12=[mesh.Y1(end-1:end, 1 : dimX2); ...
          mesh.Y2(  2, : )];
          
mesh.X23=[mesh.X2(end-1, : ); ...
          mesh.X3(  1:2, 1 : dimX2)]; 
          
mesh.Y23=[mesh.Y2(end-1, : ); ...
          mesh.Y3(  1:2, 1 : dimX2)];
      
% Corner meshes
mesh.X11=[mesh.X1(end-1:end, dimX2-1:dimX2+1); ...
          mesh.X2(2,end-1:end), mesh.X2(2,end)];
mesh.Y11=[mesh.Y1(end-1:end, dimX2-1:dimX2+1); ...
          mesh.Y2(2,end-1:end), mesh.Y2(2,end)];
mesh.X33=[mesh.X2(end-1,end-1:end), mesh.X2(end-1,end); ...
          mesh.X3(1:2,dimX2-1:dimX2+1)];
mesh.Y33=[mesh.Y2(end-1,end-1:end), mesh.Y2(end-1,end); ...
          mesh.Y3(1:2,dimX2-1:dimX2+1)];

% mesh.X3 = zeros(dimY3,dimX3);
% mesh.X3(:,(1:(dimX3-dimX2+2)/2)) = repmat(linspace(-l.d3/2,-l.d2/2,(dimX3-dimX2+2)/2),dimY3,1);
% mesh.X3(:,((dimX3-dimX2+2)/2:(dimX3+dimX2)/2)) = repmat(linspace(-l.d2/2,l.d2/2,dimX2),dimY3,1);
% mesh.X3(:,((dimX3+dimX2)/2:dimX3)) = repmat(linspace(l.d2/2,l.d3/2,(dimX3-dimX2+2)/2),dimY3,1);
% mesh.Y3 = repmat(linspace(0,-l.l3,dimY3)',1,dimX3);

%{
figure(8)
plot(mesh.X1,mesh.Y1,'k.',mesh.X2,mesh.Y2,'k.',mesh.X3,mesh.Y3,'k.')
hold on
plot([mesh.X1(1,1),mesh.X1(end,1)],[mesh.Y1(1,1),mesh.Y1(end,1)],'k-')
plot([mesh.X1(end,1),mesh.X1(end,end)],[mesh.Y1(end,1),mesh.Y1(end,end)],'k-')
plot([mesh.X1(end,end),mesh.X1(1,end)],[mesh.Y1(end,end),mesh.Y1(1,end)],'k-')
plot([mesh.X1(1,end),mesh.X1(1,1)],[mesh.Y1(1,end),mesh.Y1(1,1)],'k-')
plot([mesh.X2(1,1),mesh.X2(end,1)],[mesh.Y2(1,1),mesh.Y2(end,1)],'k-')
plot([mesh.X2(1,end),mesh.X2(end,end)],[mesh.Y2(1,end),mesh.Y2(end,end)],'k-')
plot([mesh.X2(end,1),mesh.X2(end,end)],[mesh.Y2(end,1),mesh.Y2(end,end)],'k-')
plot([mesh.X3(:,1),mesh.X3(:,1)],[mesh.Y3(:,1),mesh.Y3(:,1)],'k-')
plot([mesh.X3(end,1),mesh.X3(end,end)],[mesh.Y3(end,1),mesh.Y3(end,end)],'k-')
plot([mesh.X3(:,end),mesh.X3(:,end)],[mesh.Y3(:,end),mesh.Y3(:,end)],'k-')
plot([mesh.X3(1,end),mesh.X3(1,1)],[mesh.Y3(1,end),mesh.Y3(1,1)],'k-')
axis off
axis equal;
drawnow
%}

%%
%Number the nodes
index1=@(i, j) i + (j-1)*dimY1;
index2=@(i, j) (i-1) + (j-1)*(dimY2-2) + dimX1*dimY1;
index3=@(i, j) i + (j-1)*dimY3 + dimX1*dimY1 + dimX2*(dimY2-2);

index.block1 = index1(repmat((1:dimY1)',1,dimX1),repmat(1:dimX1,dimY1,1));

index.block3 = index3(repmat((1:dimY3)',1,dimX3),repmat(1:dimX3,dimY3,1));

index.block2 = zeros(dimY2,dimX2);
index.block2(1,:) = index.block1(end,1:dimX2);
index.block2(2:end-1,:) = index2(repmat((2:dimY2-1)',1,dimX2),repmat(1:dimX2,dimY2-2,1));
index.block2(end,:) = index.block3(1,1:dimX2);

% Number of Interfaces
index.block12 = [index.block1(end-1:end, 1 : dimX2); ...
                 index.block2( 2, : )];
       
index.block23 = [index.block2(end-1, : ); ...
                 index.block3( 1:2, 1 : dimX2)];
             
% Number of Corners
index.block11 = [index.block1(end-1:end, dimX2-1:dimX2+1); ...
                 index.block2(2, dimX2-1:end), index.block2(2, end)];
             
index.block33 = [index.block2(end-1, dimX2-1:end), index.block2(end-1, end); ...
                 index.block3(1:2, dimX2-1:dimX2+1)];




 end