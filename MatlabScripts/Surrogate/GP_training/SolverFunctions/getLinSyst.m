% This function assembles the discretization matrix from the submatrices
% for each tube. The code is set up for three tubes but can be exteded to
% more complex geometries. A seperate equation is used for the inner corner
% nodes, see A11 and A33.

function [A,B,Axn,Axs]  = getLinSyst(mesh, index, bound, dim)

boundary = bound;

[A1,B1,Axn,~] = solveFVM(mesh.X1, mesh.Y1, index.block1, boundary(1), dim);
[A2,B2] = solveFVM(mesh.X2, mesh.Y2, index.block2, boundary(2), dim);
[A3,B3,~,Axs] = solveFVM(mesh.X3, mesh.Y3, index.block3, boundary(3), dim);

A=A1+A2+A3; % Assemble with inaccurate interface

% Correct the interfaces between block 1&2 and 2&3 and inner Corners

% Delete "wrong" Matrix entries
A(index.block2(  1,:),:)=0;
A(index.block2(end,:),:)=0; 

% New entries for interfaces
A12 = solveFVM(mesh.X12, mesh.Y12, index.block12, boundary(12), dim);
A23 = solveFVM(mesh.X23, mesh.Y23, index.block23, boundary(23), dim);

% New entries for inner egdes
A11 = solveFVM(mesh.X11, mesh.Y11, index.block11, boundary(11), dim);
A33 = solveFVM(mesh.X33, mesh.Y33, index.block33, boundary(33), dim);

% Assemble with correct interfaces
A=A+A12+A23+A11+A33;
B=B1+B2+B3;

end


