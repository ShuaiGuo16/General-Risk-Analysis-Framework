% Flame describung function
function H = sflame(mesh, index, dim, l, flame, thin)

% Number of rows where flame is applied
n = 3;

% Coefficients of 1st derivate downwind scheme of first order
dy2 = mesh.Y2(1,1)-mesh.Y2(2,1);
coeff = [1.5/dy2; -2/dy2; 0.5/dy2];

dy1= mesh.Y1(end-1,1)-mesh.Y1(end,1);
flame.qref = flame.Qref/(0.25*(((l.d1-0.5*mesh.X1(1,2))/thin)^2-((0.5*mesh.X1(1,2))/thin)^2)*pi*dy1*n);
scale = (flame.gamma-1)*flame.qref/(flame.rhoref*flame.uref)*flame.FTF;

Hvalues = repmat(scale*coeff,(dim.dimX1-2)*n,1);
posi = repmat(reshape(index.block1(end-n:end-1,2:end-1),1,(dim.dimX1-2)*n),3,1);
posj = repmat(index.block2(2:4,1),(dim.dimX1-2)*n,1); % 1 at symmetry line

H = sparse(posi,posj,Hvalues,dim.total,dim.total);

end