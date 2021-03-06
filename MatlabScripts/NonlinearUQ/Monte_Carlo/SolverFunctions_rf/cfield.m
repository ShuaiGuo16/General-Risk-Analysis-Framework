


c_sound1=sqrt(1.4*287*300);
c_sound2=sqrt(1.4*287*1600);
c1=600;

yf=0.1167;  

y1=mesh.Y1(:,1); y2=mesh.Y2(:,1); y3=mesh.Y3(:,1);

VAL1=ones(dim.dimY1,dim.dimX1);
VAL2=ones(dim.dimY2,dim.dimX2);
VAL3=ones(dim.dimY3,dim.dimX3);

val1=0.5*(c_sound2-c_sound1)*(1+tanh((c1*(y1-yf))))+c_sound1;
val2=0.5*(c_sound2-c_sound1)*(1+tanh((c1*(y2-yf))))+c_sound1;
val3=0.5*(c_sound2-c_sound1)*(1+tanh((c1*(y3-yf))))+c_sound1;

for nx=1:dim.dimX1; VAL1(:,nx)=VAL1(:,nx).*val1; end
for nx=1:dim.dimX2; VAL2(:,nx)=VAL2(:,nx).*val2; end
for nx=1:dim.dimX3; VAL3(:,nx)=VAL3(:,nx).*val3; end



%%


%{
    
    figure(3)
    plot(mesh.Y1(:,1),VAL1(:,1),'o')
    hold on
    plot(mesh.Y2(:,1),VAL2(:,1),'o')
    hold on
    plot(mesh.Y3(:,1),VAL3(:,1),'o')
    
%}

% Lambda Values

boundary(1).lambda=VAL1;
boundary(2).lambda=VAL2;
boundary(3).lambda=VAL3;

boundary(12).lambda = [boundary(1).lambda(end-1:end, 1:dim.dimX2); ...
    boundary(2).lambda( 2 , :)];

boundary(23).lambda = [boundary(2).lambda(end-1, :); ...
    boundary(3).lambda( 1:2 , 1:dim.dimX2)];

boundary(11).lambda = [boundary(1).lambda(end-1:end, dim.dimX2-1:dim.dimX2+1); ...
    boundary(2).lambda(2, dim.dimX2-1:end), boundary(2).lambda(2, end)];

boundary(33).lambda = [boundary(2).lambda(end-1, dim.dimX2-1:end), boundary(2).lambda(end-1, end); ...
    boundary(3).lambda(1:2, dim.dimX2-1:dim.dimX2+1)];
